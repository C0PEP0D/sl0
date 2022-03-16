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
			sMembers.resize(Type::size());
			sMembers[0] = 0.0;
			for(std::size_t i = 1; i < Type::size(); i++){ 
				const auto x1 = sMemberStep->cX(cMemberState(state.data(), i));
				const auto x0 = sMemberStep->cX(cMemberState(state.data(), i - 1));
				l += (x1 - x0).norm();
				sMembers[i] = l;
			}
			for(std::size_t i = 1; i < Type::size(); i++){ 
				sMembers[i] /= l;
			}
			sMesh = std::make_shared<TypeMesh>(sMembers, false);
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
		TypeStateVectorDynamic cState(const double* pState, const double& p_s) const {
			const TypeVector<1> s(p_s);
			// interpolate for each coordinates
			TypeStateVectorDynamic state(sMemberStep->stateSize());
			for(std::size_t j = 0; j < sMemberStep->stateSize(); j++) {
				state[j] = p0l::lagrangeMeshPoint<TypeMesh, TypeContainer, double, TypeVector<1>, TypeRef, TypeMeshSub>(sMesh, interpolationData[j], s, interpolationOrder + 1, false);
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
			return memberIndex * 1.0 / (Type::size() - 1);
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
			sMembers.resize(Type::size());
			sMembers[0] = 0.0;
			for(std::size_t i = 1; i < Type::size(); i++){ 
				const auto x1 = sMemberStep->cX(cMemberState(state.data(), i));
				const auto x0 = sMemberStep->cX(cMemberState(state.data(), i - 1));
				l += (x1 - x0).norm();
				sMembers[i] = l;
			}
			for(std::size_t i = 1; i < Type::size(); i++){ 
				sMembers[i] /= l;
			}
			sMesh = std::make_shared<TypeMesh>(sMembers, false);
		}
	public:
		// parameters
		double dl;
		unsigned int interpolationOrder;
		// shape description updated at each time step
		std::vector<double> sMembers;
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
