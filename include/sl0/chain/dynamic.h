#ifndef SL0_CHAIN_DYNAMIC_H
#define SL0_CHAIN_DYNAMIC_H
#pragma once

// std includes
#include <memory> // shared_ptr
#include <vector>
#include <numeric> // iota
#include <execution>
#include <tuple>
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
		StepChainDynamic(std::shared_ptr<TypeMemberStep> p_sMemberStep, const double& p_dl, const unsigned int& p_interpolationOrder, const bool p_closed) : Type(p_sMemberStep), dl(p_dl), interpolationOrder(p_interpolationOrder), interpolationData(TypeMemberStep::StateSize), closed(p_closed) {
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
			double newLength = 0.0;
			meshGrid.resize(Type::size() + int(closed));
			meshGrid[0] = 0.0;
			for(std::size_t i = 1; i < Type::size(); i++){ 
				const auto x1 = sMemberStep->cX(cMemberState(state.data(), i));
				const auto x0 = sMemberStep->cX(cMemberState(state.data(), i - 1));
				newLength += (x1 - x0).norm();
				meshGrid[i] = newLength;
			}
			if(closed) {
				const TypeSpaceVector x1 = sMemberStep->cX(cMemberState(state.data(), 0));
				const TypeSpaceVector x0 = sMemberStep->cX(cMemberState(state.data(), Type::size() - 1));
				newLength += (x1 - x0).norm();
				meshGrid[Type::size()] = newLength;
			}
			for(std::size_t i = 1; i < meshGrid.size(); i++){ 
				meshGrid[i] /= newLength;
			}
			sMesh = std::make_shared<TypeMesh>(meshGrid, closed);
			// update size
			const std::size_t newSize = std::ceil(newLength/dl) + int(!closed);
			const double newDs = 1.0 / (newSize - int(!closed));
			// manage data
			const int difference = newSize - Type::size();
			if(difference < 0) {
				for(int i = -1; i >= difference; i--) {
					Type::removeMember(state);
				}
			} else if(difference > 0) {
				for(int i = 0; i < difference; i++) {
					Type::addMember(state);
				}
			}
			// interpolation
			std::for_each(std::execution::par_unseq, memberIndexs.cbegin(), memberIndexs.cend(), [this, &state, newDs](const unsigned int& memberIndex){ 
				const TypeVector<1> s(memberIndex * newDs);
				double* pMemberState = memberState(state.data(), memberIndex);
				for(std::size_t i = 0; i < TypeMemberStep::StateSize; i++) {
					pMemberState[i] = p0l::lagrangeMeshPoint<TypeMesh, TypeContainer, double, TypeVector<1>, TypeRef, TypeMeshSub>(sMesh, interpolationData[i], s, interpolationOrder + 1, closed); 
				}
			});
			// update
			Type::update(state, t);
		}
	public:
		TypeStateVectorDynamic cState(const double* pState, const double& p_s) const {
			const TypeVector<1> s(p_s);
			// interpolate for each coordinates
			TypeStateVectorDynamic state(sMemberStep->stateSize());
			for(std::size_t j = 0; j < sMemberStep->stateSize(); j++) {
				state[j] = p0l::lagrangeMeshPoint<TypeMesh, TypeContainer, double, TypeVector<1>, TypeRef, TypeMeshSub>(sMesh, interpolationData[j], s, interpolationOrder + 1, closed);
			}
			return state;
		}
	public:
		TypeVector<DIM> cX(const double* pState, const double& p_s) const {
			const TypeStateVectorDynamic state = cState(pState, p_s);
			return sMemberStep->cX(state.data());
		}
	public:
		double cS(const unsigned int memberIndex) const {
			return double(memberIndex) / (Type::size() - int(!closed)); // TODO: deal with closed everywhere
		}
	public:
		double length(const double* pState) const {
			double l = 0.0;
			for(std::size_t i = 1; i < Type::size(); i++){ 
				const TypeSpaceVector x1 = sMemberStep->cX(cMemberState(pState, i));
				const TypeSpaceVector x0 = sMemberStep->cX(cMemberState(pState, i-1));
				l += (x1 - x0).norm();
			}
			if(closed) {
				const TypeSpaceVector x1 = sMemberStep->cX(cMemberState(pState, 0));
				const TypeSpaceVector x0 = sMemberStep->cX(cMemberState(pState, Type::size() - 1));
				l += (x1 - x0).norm();
			}
			return l;
		}
	public:
		struct ClosestPointData {
			public:
				ClosestPointData() {
				}
			public:
				TypeSpaceVector x;
				double distance;
				double s;
		};
		
		ClosestPointData closest(const double* pState, const TypeSpaceVector& x) const {
			std::vector<ClosestPointData> closestPointsData(Type::size() - 1);
			for(std::size_t i = 1; i < Type::size(); i++){ 
				const TypeSpaceVector x1 = sMemberStep->cX(cMemberState(pState, i));
				const TypeSpaceVector x0 = sMemberStep->cX(cMemberState(pState, i-1));
				// segment
				const TypeSpaceVector segment = x1 - x0;
				const TypeSpaceVector lSegment = segment.norm();
				const TypeSpaceVector dSegment = (x1 - x0) / lSegment;
				const double sSegment = (x - x0).dot(dSegment);
				// closest point on segment
				double s;
				TypeSpaceVector xSegment;
				if(sSegment > 0.0) {
					if (sSegment < lSegment) {
						xSegment = x0 + dSegment * sSegment;
						s = cS(i - 1) + sSegment/lSegment * (cS(i) - cS(i - 1));
					} else {
						xSegment = x1;
						s = cS(i);
					}
				} else {
					xSegment = x0;
					s = cS(i - 1);
				}
				closestPointsData[i - 1].x = xSegment;
				closestPointsData[i - 1].distance = (x - xSegment).norm();
				closestPointsData[i - 1].s = s;
			}
			std::sort(closestPointsData.begin(), closestPointsData.end(), [](ClosestPointData a, ClosestPointData b) {
				return a.distance < b.distance;
			});
			return closestPointsData.front();
		}
	public:
		struct IntersectionData {
			public:
				IntersectionData() {
				}
			public:
				TypeSpaceVector point;
				unsigned int i;
				unsigned int j;
				double iS;
				double jS;
		};

		IntersectionData firstIntersection(const double* pState) const {
			IntersectionData firstIntersectionData;
			firstIntersectionData.point = TypeSpaceVector::Zero();
			firstIntersectionData.i = 0;
			firstIntersectionData.j = 0;
			firstIntersectionData.iS = 0.0;
			firstIntersectionData.jS = 0.0;
			// search for intersections
			if(Type::size() > 3) {
				for(unsigned int i = 0; i < Type::size() - 1; i++){
					const TypeSpaceVector iX0 = sMemberStep->cX(cMemberState(pState, i));
					const TypeSpaceVector iX1 = sMemberStep->cX(cMemberState(pState, i + 1));
					const TypeSpaceVector iSegment = iX1 - iX0;
					const double iLength = iSegment.norm();
					const TypeSpaceVector iDir = iSegment / iLength;
					Line iLine = Line::Through(TypeVector<2>(iX0(0), iX0(1)), TypeVector<2>(iX1(0), iX1(1)));
					for(unsigned int j = i + 2; j < Type::size() - 1; j++) {
						const TypeSpaceVector jX0 = sMemberStep->cX(cMemberState(pState, j));
						const TypeSpaceVector jX1 = sMemberStep->cX(cMemberState(pState, j + 1));
						const TypeSpaceVector jSegment = jX1 - jX0;
						const double jLength = jSegment.norm();
						const TypeSpaceVector jDir = jSegment / jLength;
						Line jLine = Line::Through(TypeVector<2>({jX0(0), jX0(1)}), TypeVector<2>({jX1(0), jX1(1)}));
						const TypeVector<2> point2d = iLine.intersection(jLine);
						const TypeSpaceVector point = TypeSpaceVector({point2d(0), point2d(1), 0.25 * (iX0(2) + iX1(2) + jX0(2) + jX1(2))});
						if ((point - iX0).dot(iDir) > 0.0 && (point - iX0).dot(iDir) < iLength && (point - jX0).dot(jDir) > 0.0 && (point - jX0).dot(jDir) < jLength) {
							firstIntersectionData.point = point;
							firstIntersectionData.i = i;
							firstIntersectionData.j = j;
							firstIntersectionData.iS = cS(i) + (point - iX0).dot(iDir)/iLength * (cS(i + 1) - cS(i));
							firstIntersectionData.jS = cS(j) + (point - jX0).dot(jDir)/jLength * (cS(j + 1) - cS(j));
							return firstIntersectionData;
						}
					}
					if(closed && i > 0 && i < Type::size() - 2){
						const TypeSpaceVector jX0 = sMemberStep->cX(cMemberState(pState, Type::size() - 1));
						const TypeSpaceVector jX1 = sMemberStep->cX(cMemberState(pState, 0));
						const TypeSpaceVector jSegment = jX1 - jX0;
						const double jLength = jSegment.norm();
						const TypeSpaceVector jDir = jSegment / jLength;
						Line jLine = Line::Through(TypeVector<2>({jX0(0), jX0(1)}), TypeVector<2>({jX1(0), jX1(1)}));
						const TypeVector<2> point2d = iLine.intersection(jLine);
						const TypeSpaceVector point = TypeSpaceVector({point2d(0), point2d(1), 0.25 * (iX0(2) + iX1(2) + jX0(2) + jX1(2))});
						if ((point - iX0).dot(iDir) > 0.0 && (point - iX0).dot(iDir) < iLength && (point - jX0).dot(jDir) > 0.0 && (point - jX0).dot(jDir) < jLength) {
							firstIntersectionData.point = point;
							firstIntersectionData.i = i;
							firstIntersectionData.j = Type::size() - 1;
							firstIntersectionData.iS = cS(i) + (point - iX0).dot(iDir)/iLength * (cS(i + 1) - cS(i));
							firstIntersectionData.jS = cS(Type::size() - 1) + (point - jX0).dot(jDir)/jLength * (1.0 - cS(Type::size() - 1));
							return firstIntersectionData;
						}
					}
				}
			}
			return firstIntersectionData;
		}

		std::vector<IntersectionData> intersections(const double* pState) const {
			std::vector<IntersectionData> intersectionsData;
			for(unsigned int i = 0; i < Type::size() - 1; i++){
				const TypeSpaceVector iX0 = sMemberStep->cX(cMemberState(pState, i));
				const TypeSpaceVector iX1 = sMemberStep->cX(cMemberState(pState, i + 1));
				const TypeSpaceVector iSegment = iX1 - iX0;
				const double iLength = iSegment.norm();
				const TypeSpaceVector iDir = iSegment / iLength;
				Line iLine = Line::Through(TypeVector<2>(iX0(0), iX0(1)), TypeVector<2>(iX1(0), iX1(1)));
				for(unsigned int j = i + 2; j < Type::size() - 1; j++) {
					const TypeSpaceVector jX0 = sMemberStep->cX(cMemberState(pState, j));
					const TypeSpaceVector jX1 = sMemberStep->cX(cMemberState(pState, j + 1));
					const TypeSpaceVector jSegment = jX1 - jX0;
					const double jLength = jSegment.norm();
					const TypeSpaceVector jDir = jSegment / jLength;
					Line jLine = Line::Through(TypeVector<2>({jX0(0), jX0(1)}), TypeVector<2>({jX1(0), jX1(1)}));
					const TypeVector<2> point2d = iLine.intersection(jLine);
					const TypeSpaceVector point = TypeSpaceVector({point2d(0), point2d(1), 0.25 * (iX0(2) + iX1(2) + jX0(2) + jX1(2))});
					if ((point - iX0).dot(iDir) > 0.0 && (point - iX0).dot(iDir) < iLength && (point - jX0).dot(jDir) > 0.0 && (point - jX0).dot(jDir) < jLength) {
						IntersectionData newIntersectionData;
						newIntersectionData.point = point;
						newIntersectionData.i = i;
						newIntersectionData.j = j;
						newIntersectionData.iS = cS(i) + (point - iX0).dot(iDir)/iLength * (cS(i + 1) - cS(i));
						newIntersectionData.jS = cS(j) + (point - jX0).dot(jDir)/jLength * (cS(j + 1) - cS(j));
						intersectionsData.push_back(newIntersectionData);
					}
				}
				if(closed && i > 0 && i < Type::size() - 2){
					const TypeSpaceVector jX0 = sMemberStep->cX(cMemberState(pState, Type::size() - 1));
					const TypeSpaceVector jX1 = sMemberStep->cX(cMemberState(pState, 0));
					const TypeSpaceVector jSegment = jX1 - jX0;
					const double jLength = jSegment.norm();
					const TypeSpaceVector jDir = jSegment / jLength;
					Line jLine = Line::Through(TypeVector<2>({jX0(0), jX0(1)}), TypeVector<2>({jX1(0), jX1(1)}));
					const TypeVector<2> point2d = iLine.intersection(jLine);
					const TypeSpaceVector point = TypeSpaceVector({point2d(0), point2d(1), 0.25 * (iX0(2) + iX1(2) + jX0(2) + jX1(2))});
					if ((point - iX0).dot(iDir) > 0.0 && (point - iX0).dot(iDir) < iLength && (point - jX0).dot(jDir) > 0.0 && (point - jX0).dot(jDir) < jLength) {
						IntersectionData newIntersectionData;
						newIntersectionData.point = point;
						newIntersectionData.i = i;
						newIntersectionData.j = Type::size() - 1;
						newIntersectionData.iS = cS(i) + (point - iX0).dot(iDir)/iLength * (cS(i + 1) - cS(i));
						newIntersectionData.jS = cS(Type::size() - 1) + (point - jX0).dot(jDir)/jLength * (1.0 - cS(Type::size() - 1));
						intersectionsData.push_back(newIntersectionData);
					}
				}
			}
			return intersectionsData;
		}
	public:
		void registerState(const std::vector<double>& state) override {
			Type::registerState(state);
			// update
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
			meshGrid.resize(Type::size() + int(closed));
			meshGrid[0] = 0.0;
			for(std::size_t i = 1; i < Type::size(); i++){ 
				const auto x1 = sMemberStep->cX(cMemberState(state.data(), i));
				const auto x0 = sMemberStep->cX(cMemberState(state.data(), i - 1));
				l += (x1 - x0).norm();
				meshGrid[i] = l;
			}
			if(closed) {
				const TypeSpaceVector x1 = sMemberStep->cX(cMemberState(state.data(), 0));
				const TypeSpaceVector x0 = sMemberStep->cX(cMemberState(state.data(), Type::size() - 1));
				l += (x1 - x0).norm();
				meshGrid[Type::size()] = l;
			}
			for(std::size_t i = 1; i < meshGrid.size(); i++){ 
				meshGrid[i] /= l;
			}
			sMesh = std::make_shared<TypeMesh>(meshGrid, closed);
		}
	public:
		// parameters
		double dl;
		unsigned int interpolationOrder;
		bool closed;
		// shape description updated at each time step
		std::vector<double> meshGrid;
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
