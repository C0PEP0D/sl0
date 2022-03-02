#ifndef SL0_CHAIN_DYNAMIC_H
#define SL0_CHAIN_DYNAMIC_H
#pragma once

// std includes
#include <memory> // shared_ptr
#include <vector>
#include <numeric> // iota
#include <execution>
// module includes
#include "sl0/group/dynamic.h"
#include "p0l/interpolation.h"
#include "m0sh/non_uniform.h"
#include "m0sh/structured_sub.h"
// test
#include <Eigen/Dense> // TEST

namespace sl0 {

template<template<int> typename TypeVector, unsigned int DIM, template<typename...> class TypeView, template<typename...> class TypeRef, typename TypeMemberStep>
class StepChainDynamic : public StepGroupDynamicHomogeneous<TypeVector, DIM, TypeView, TypeMemberStep> {
	public:
		using Type = StepGroupDynamicHomogeneous<TypeVector, DIM, TypeView, TypeMemberStep>;
		using typename Type::TypeSpaceVector;
		using typename Type::TypeStateVectorDynamic;
	public:
		using Line = Eigen::Hyperplane<double, 2>; // TEST
	public:
		template<typename ...Args>
		using TypeContainer = std::vector<Args...>;
		using TypeMesh = m0sh::NonUniform<TypeVector<1>, TypeRef, TypeContainer>;
		using TypeMeshSub = m0sh::StructuredSub<TypeVector<1>, TypeRef, TypeContainer>;
	public:
		StepChainDynamic(std::shared_ptr<TypeMemberStep> p_sMemberStep, const double& p_dl, const unsigned int& p_interpolationOrder) : Type(p_sMemberStep), dl(p_dl), interpolationOrder(p_interpolationOrder), interpolationData(TypeMemberStep::StateSize) {
		}

		void update(std::vector<double>& state, const double& t) override {
			// interpolation data
			for(std::size_t i = 0; i < TypeMemberStep::StateSize; i++) {
				interpolationData[i].resize(Type::size());
			}
			for(std::size_t i = 0; i < Type::size(); i++) {
				const double* pMemberState = cMemberState(state.data(), i);
				for(std::size_t j = 0; j < TypeMemberStep::StateSize; j++) {
					interpolationData[j][i] = pMemberState[j];
				}
			}
			// mesh
			double l = 0.0;
			TypeContainer<double> gridPoints(Type::size());
			gridPoints[0] = 0.0;
			for(std::size_t i = 1; i < Type::size(); i++){ 
				const auto x1 = sMemberStep->cX(cMemberState(state.data(), i));
				const auto x0 = sMemberStep->cX(cMemberState(state.data(), i - 1));
				l += (x1 - x0).norm();
				gridPoints[i] = l;
			}
			for(std::size_t i = 1; i < Type::size(); i++){ 
				gridPoints[i] /= l;
			}
			sMesh = std::make_shared<TypeMesh>(gridPoints, false);
			// update size
			const double newLength = length(state.data());
			const std::size_t newSize = std::ceil(newLength/dl);
			const double newDs = 1.0 / (newSize - 1);
			// manage data
			const int difference = newSize - Type::size();
			if(difference < 0) {
				for(int i = -1; i >= difference; i--) {
					Type::removeMember(state);
				}
			} else if(difference > 0) {
				for(int i = 0; i < difference; i++) {
					// create new member
					std::shared_ptr<TypeMemberStep> sNewMemberStep = std::make_shared<TypeMemberStep>(*sMemberStep);
					// add member
					Type::addMember(state);
				}
			}
			// interpolation
			std::for_each(std::execution::par_unseq, memberIndexs.cbegin(), memberIndexs.cend(), [this, &state, newDs](const unsigned int& memberIndex){ 
				const TypeVector<1> s(memberIndex * newDs);
				double* pMemberState = memberState(state.data(), memberIndex);
				for(std::size_t i = 0; i < TypeMemberStep::StateSize; i++) {
					pMemberState[i] = p0l::lagrangeMeshPoint<TypeMesh, TypeContainer, double, TypeVector<1>, TypeRef, TypeMeshSub>(sMesh, interpolationData[i], s, interpolationOrder + 1, false);
				}
			});
			// update
			Type::update(state, t);
		}
	public:
		TypeVector<DIM> cX(const double* pState, const double& p_s) const {
			const TypeVector<1> s(p_s);
			// coordinates
			TypeContainer<TypeContainer<double>> tmpCoordinates(DIM, TypeContainer<double>(Type::size())); // TODO: do this better
			for(std::size_t i = 0; i < Type::size(); i++) {
				const TypeVector<DIM> position = sMemberStep->cX(cMemberState(pState, i));
				for(std::size_t j = 0; j < DIM; j++) {
					tmpCoordinates[j][i] = position[j];
				}
			}
			// mesh
			double l = 0.0;
			TypeContainer<double> gridPoints(Type::size());
			gridPoints[0] = 0.0;
			for(std::size_t i = 1; i < Type::size(); i++){ 
				const TypeSpaceVector x1 = sMemberStep->cX(cMemberState(pState, i));
				const TypeSpaceVector x0 = sMemberStep->cX(cMemberState(pState, i - 1));
				l += (x1 - x0).norm();
				gridPoints[i] = l;
			}
			for(std::size_t i = 1; i < Type::size(); i++){ 
				gridPoints[i] /= l;
			}
			std::shared_ptr<TypeMesh> sTmpMesh = std::make_shared<TypeMesh>(gridPoints, false);
			// interpolate for each coordinates
			TypeVector<DIM> x;
			for(std::size_t j = 0; j < DIM; j++) {
				x[j] = p0l::lagrangeMeshPoint<TypeMesh, TypeContainer, double, TypeVector<1>, TypeRef, TypeMeshSub>(sTmpMesh, tmpCoordinates[j], s, interpolationOrder + 1, false);
			}
			return x;
		}
	public:
		double length(const double* pState) const {
			double l = 0.0;
			for(std::size_t i = 1; i < Type::size(); i++){ 
				const TypeSpaceVector x1 = sMemberStep->cX(cMemberState(pState, i));
				const TypeSpaceVector x0 = sMemberStep->cX(cMemberState(pState, i-1));
				l += (x1 - x0).norm();
			}
			return l;
		}
	public:
		TypeSpaceVector closest(const double* pState, const TypeSpaceVector& x) const {
			std::vector<std::pair<TypeSpaceVector, double>> xSegments(Type::size() - 1);
			for(std::size_t i = 1; i < Type::size(); i++){ 
				const TypeSpaceVector x1 = sMemberStep->cX(cMemberState(pState, i));
				const TypeSpaceVector x0 = sMemberStep->cX(cMemberState(pState, i-1));
				// segment
				const TypeSpaceVector segment = x1 - x0;
				const TypeSpaceVector lSegment = segment.norm();
				const TypeSpaceVector dSegment = (x1 - x0) / lSegment;
				const double sSegment = (x - x0).dot(dSegment);
				// closest point on segment
				TypeSpaceVector xSegment;
				if(sSegment > 0.0) {
					if (sSegment < lSegment) {
						xSegment = x0 + dSegment * sSegment;
					} else {
						xSegment = x1;
					}
				} else {
					xSegment = x0;
				}
				xSegments[i-1].first = xSegment;
				xSegments[i-1].second = (x - xSegment).norm();
			}
			std::sort(xSegments.begin(), xSegments.end(), [](std::pair<TypeSpaceVector, double> a, std::pair<TypeSpaceVector, double> b) {
				return a.second < b.second;
			});
			return xSegments.front().first;
		}
	public:
		std::vector<TypeSpaceVector> intersections(const double* pState) const {
			std::vector<TypeSpaceVector> points;
			for(unsigned int i = 0; i < Type::size() - 1; i++){
				const TypeSpaceVector iX0 = sMemberStep->cX(cMemberState(pState, i));
				const TypeSpaceVector iX1 = sMemberStep->cX(cMemberState(pState, i + 1));
				const TypeSpaceVector iSegment = iX1 - iX0;
				const TypeSpaceVector iLength = iSegment.norm();
				const TypeSpaceVector iDir = iSegment / iLength;
				Line iLine = Line::Through(TypeVector<2>(iX0(0), iX0(1)), TypeVector<2>(iX1(0), iX1(1)));
				for(unsigned int j = i + 1; j < Type::size() - 1; j++) {
					const TypeSpaceVector jX0 = sMemberStep->cX(cMemberState(pState, j));
					const TypeSpaceVector jX1 = sMemberStep->cX(cMemberState(pState, j + 1));
					const TypeSpaceVector jSegment = jX1 - jX0;
					const TypeSpaceVector jLength = jSegment.norm();
					const TypeSpaceVector jDir = jSegment / jLength;
					Line jLine = Line::Through(TypeVector<2>({jX0(0), jX0(1)}), TypeVector<2>({jX1(0), jX1(1)}));
					// Dirty intersection code
					const TypeVector<2> point = iLine.intersection(jLine);
					if ((point - iX0).dot(iDir) > 0.0 && (point - iX0).dot(iDir) < iLength && (point - jX0).dot(jDir) > 0.0 && (point - jX0).dot(jDir) < jLength) {
						points.emplace_back({point(0), point(1)});
					}
				}
			}
			return points;
		}
	public:
		// parameters
		double dl;
		unsigned int interpolationOrder;
		// shape description
		std::shared_ptr<TypeMesh> sMesh;
		std::vector<std::vector<double>> interpolationData;
	public:
		using Type::memberIndexs;
	public:
		using Type::sMemberStep;
	public:
		using Type::cMemberState;
		using Type::memberState;
	public:
		using Type::sMemberSteps;
};

template<template<int> typename TypeVector, unsigned int DIM, template<typename...> class TypeView, template<typename...> class TypeRef, typename TypeMemberStep, typename TypeSolver>
class ChainDynamic : public ObjectDynamic<TypeSolver, StepChainDynamic<TypeVector, DIM, TypeView, TypeRef, TypeMemberStep>> {
	public:
		using TypeStep = StepChainDynamic<TypeVector, DIM, TypeView, TypeRef, TypeMemberStep>;
	public:
		template<typename ...Args>
		using TypeContainer = std::vector<Args...>;
	public:
		ChainDynamic(const std::shared_ptr<TypeMemberStep>& p_sMemberStep, const double& p_dl, const unsigned int& p_interpolationOrder) : ObjectDynamic<TypeSolver, TypeStep>(std::make_shared<TypeStep>(p_sMemberStep, p_dl, p_interpolationOrder)) {
		}
	public:
		// Inherited
		using ObjectDynamic<TypeSolver, TypeStep>::sSolver;
		using ObjectDynamic<TypeSolver, TypeStep>::sStep;
		using ObjectDynamic<TypeSolver, TypeStep>::state;
		using ObjectDynamic<TypeSolver, TypeStep>::t;
};

}

#endif
