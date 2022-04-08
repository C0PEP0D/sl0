#ifndef SL0_GROUP_DYNAMIC_H
#define SL0_GROUP_DYNAMIC_H
#pragma once

// std includes
#include <memory> // shared_ptr
#include <vector>
#include <numeric> // iota
#include <execution>
// module includes
#include "sl0/object.h"

// TODO: add possibility to change execution policy, seq, par, par_unseq

namespace sl0 {

template<template<int> typename TypeVector, unsigned int DIM>
class StepGroupDynamicBase : public StepObject<TypeVector, DIM> {
	public:
		using Type = StepObject<TypeVector, DIM>;
		using typename Type::TypeSpaceVector;
		using typename Type::TypeStateVectorDynamic;
	public:
		StepGroupDynamicBase() {
		}
	public:
		// member management
		virtual void removeMember(std::vector<double>& state, const unsigned int& memberIndex) = 0;
		virtual void removeMember(std::vector<double>& state) = 0;
		virtual void updateMemberSize(std::vector<double>& state, const unsigned int& memberIndex, const unsigned int& stateSize) = 0;
	public:
		// getters
		virtual std::size_t size() const = 0;
		virtual const double* cMemberState(const double* pState, const std::size_t& memberIndex) const = 0;
		virtual double* memberState(double* pState, const std::size_t& memberIndex) const = 0;
};

template<template<int> typename TypeVector, unsigned int DIM, template<typename...> class TypeView, typename _TypeMemberStep>
class StepGroupDynamic : public StepGroupDynamicBase<TypeVector, DIM> {
	public:
		using Type = StepGroupDynamicBase<TypeVector, DIM>;
		using typename Type::TypeSpaceVector;
		using typename Type::TypeStateVectorDynamic;
		using TypeMemberStep = _TypeMemberStep;
	public:
		StepGroupDynamic() {
		}
	public:
		// Main
		TypeStateVectorDynamic operator()(const double* pState, const double& t) const override {
			// Init dState
			TypeStateVectorDynamic dState;
			dState.resize(stateSize());
			// Set dState
			std::for_each(std::execution::par_unseq, memberIndexs.cbegin(), memberIndexs.cend(), [this, pState, t, &dState](const unsigned int& memberIndex){ 
				TypeView<TypeStateVectorDynamic> dMemberState(memberState(dState.data(), memberIndex), memberStateSizes[memberIndex]); 
				dMemberState = (*sMemberSteps[memberIndex])(cMemberState(pState, memberIndex), t); 
			});
			// Return result
			return dState;
		}

		void update(std::vector<double>& state, const double& t) override {
			std::for_each(std::execution::par_unseq, memberIndexs.cbegin(), memberIndexs.cend(), [this, &state, t](const unsigned int& memberIndex){
				sMemberSteps[memberIndex]->update(memberState(state.data(), memberIndex), t);
			});
		}

		unsigned int stateSize() const override {
			if (not memberStateIndexs.empty()) {
				return memberStateIndexs.back() + memberStateSizes.back();
			} else {
				return 0;
			}
		}
	public:
		// member management
		void addMember(std::vector<double>& state, std::shared_ptr<TypeMemberStep> sMemberStep, const unsigned int& stateSize) {
			// state resize
			state.resize(state.size() + stateSize);
			// register member
			if (not memberStateIndexs.empty()) {
				memberIndexs.push_back(memberIndexs.back() + 1);
				memberStateIndexs.push_back(memberStateIndexs.back() + memberStateSizes.back());
			} else {
				memberIndexs.push_back(0);
				memberStateIndexs.push_back(0);
			}
			sMemberSteps.push_back(sMemberStep);
			memberStateSizes.push_back(stateSize);
		}

		void removeMember(std::vector<double>& state, const unsigned int& memberIndex) override {
			// state resize
			std::copy(state.begin() + memberIndex + memberStateSizes[memberIndex], state.end(), state.begin() + memberIndex);
			state.resize(state.size() - memberStateSizes[memberIndex]);
			// unregister member
			memberIndexs.pop_back();
			sMemberSteps.erase(sMemberSteps.begin() + memberIndex);
			for(unsigned int i = memberIndex + 1; i < sMemberSteps.size(); i++) {
				memberStateIndexs[i] -= memberStateSizes[memberIndex];
			}
			memberStateIndexs.erase(memberStateIndexs.begin() + memberIndex);
			memberStateSizes.erase(memberStateSizes.begin() + memberIndex);
		}


		void removeMember(std::vector<double>& state) override {
			// state resize
			state.resize(state.size() - memberStateSizes.back());
			// unregister member
			memberIndexs.pop_back();
			sMemberSteps.pop_back();
			memberStateIndexs.pop_back();
			memberStateSizes.pop_back();
		}

		void updateMemberSize(std::vector<double>& state, const unsigned int& memberIndex, const unsigned int& stateSize) override {
			const int difference = stateSize - memberStateSizes[memberIndex];
			if(difference > 0) {
				state.resize(state.size() + difference);
				std::copy(state.begin() + memberIndex + memberStateSizes[memberIndex], state.end() - difference, state.begin() + stateSize);
			} else {
				std::copy(state.begin() + memberIndex + memberStateSizes[memberIndex], state.end(), state.begin() + stateSize);
				state.resize(state.size() + difference);
			}
			// register size change
			memberStateSizes[memberIndex] += difference;
			for(unsigned int i = memberIndex + 1; i < sMemberSteps.size(); i++) {
				memberStateIndexs[i] += difference;
			}
		}
	public:
		// getters
		std::size_t size() const override {
			return sMemberSteps.size();
		}

		const double* cMemberState(const double* pState, const std::size_t& memberIndex) const override {
			return pState + memberStateIndexs[memberIndex];
		}

		double* memberState(double* pState, const std::size_t& memberIndex) const override {
			return pState + memberStateIndexs[memberIndex];
		}
	public:
		std::vector<TypeSpaceVector> positions(const double* pState) const override {
			std::vector<TypeSpaceVector> result;
			//result.reserve(size() * sMemberSteps[0]->positions(cMemberState(state, 0)).size());
			for(unsigned int i = 0; i < size(); i++) {
				const std::vector<TypeSpaceVector> memberPositions = sMemberSteps[i]->positions(cMemberState(pState, i));
				result.insert(result.end(), memberPositions.begin(), memberPositions.end());
			}
			return result;
		}
	public:
		std::vector<unsigned int> memberIndexs;
		std::vector<std::shared_ptr<TypeMemberStep>> sMemberSteps;
		std::vector<unsigned int> memberStateSizes;
		std::vector<unsigned int> memberStateIndexs;
};

template<template<int> typename TypeVector, unsigned int DIM, template<typename...> class TypeView, typename TypeMemberStep, typename TypeSolver>
class GroupDynamic : public ObjectDynamic<TypeSolver, StepGroupDynamic<TypeVector, DIM, TypeView, TypeMemberStep>> {
	public:
		using TypeStep = StepGroupDynamic<TypeVector, DIM, TypeView, TypeMemberStep>;
	public:
		GroupDynamic() : ObjectDynamic<TypeSolver, TypeStep>(std::make_shared<TypeStep>()) {
		}
	public:
		// Inherited properties
		using ObjectDynamic<TypeSolver, TypeStep>::sSolver;
		using ObjectDynamic<TypeSolver, TypeStep>::sStep;
		using ObjectDynamic<TypeSolver, TypeStep>::state;
		using ObjectDynamic<TypeSolver, TypeStep>::t;
};

template<template<int> typename TypeVector, unsigned int DIM>
class StepGroupDynamicHomogeneousBase : public StepGroupDynamicBase<TypeVector, DIM> {
	public:
		StepGroupDynamicHomogeneousBase() {
		}
	public:
		virtual void addMember(std::vector<double>& state) = 0;
	public:
		virtual void registerState(const std::vector<double>& state) = 0;
};

template<template<int> typename TypeVector, unsigned int DIM, template<typename...> class TypeView, typename _TypeMemberStep>
class StepGroupDynamicHomogeneous : public StepGroupDynamic<TypeVector, DIM, TypeView, _TypeMemberStep>, public StepGroupDynamicHomogeneousBase<TypeVector, DIM> {
	public:
		using TypeGroupDynamic = StepGroupDynamic<TypeVector, DIM, TypeView, _TypeMemberStep>;
		using typename TypeGroupDynamic::TypeSpaceVector;
		using typename TypeGroupDynamic::TypeStateVectorDynamic;
		using typename TypeGroupDynamic::TypeMemberStep;
	public:
		StepGroupDynamicHomogeneous(std::shared_ptr<TypeMemberStep> p_sMemberStep) : sMemberStep(p_sMemberStep) {
		}
	public:
		// Main
		TypeStateVectorDynamic operator()(const double* pState, const double& t) const override {
			return TypeGroupDynamic::operator()(pState, t);
		}

		void update(std::vector<double>& state, const double& t) override {
			TypeGroupDynamic::update(state, t);
		}

		unsigned int stateSize() const override {
			return TypeGroupDynamic::stateSize();
		}
	public:
		// member management
		void addMember(std::vector<double>& state, std::shared_ptr<TypeMemberStep> sMemberStep, const unsigned int& stateSize) {
			TypeGroupDynamic::addMember(state, sMemberStep, stateSize);
		}
		void removeMember(std::vector<double>& state, const unsigned int& memberIndex) override {
			TypeGroupDynamic::removeMember(state, memberIndex);
		}
		void removeMember(std::vector<double>& state) override {
			TypeGroupDynamic::removeMember(state);
		}
		void updateMemberSize(std::vector<double>& state, const unsigned int& memberIndex, const unsigned int& stateSize) override {
			TypeGroupDynamic::updateMemberSize(state, memberIndex, stateSize);
		}
	public:
		// getters
		std::size_t size() const override {
			return TypeGroupDynamic::size();
		}

		const double* cMemberState(const double* pState, const std::size_t& memberIndex) const {
			return TypeGroupDynamic::cMemberState(pState, memberIndex);
		}

		double* memberState(double* pState, const std::size_t& memberIndex) const {
			return TypeGroupDynamic::memberState(pState, memberIndex);
		}
	public:
		std::vector<TypeSpaceVector> positions(const double* pState) const override {
			return TypeGroupDynamic::positions(pState);
		}
	public:
		void addMember(std::vector<double>& state) override {
			TypeGroupDynamic::addMember(state, std::make_shared<TypeMemberStep>(*sMemberStep), sMemberStep->stateSize());
		}
	public:
		void registerState(const std::vector<double>& state) override {
			const std::size_t n = state.size() / sMemberStep->stateSize();
			// register and unregister if state size does not correspond to size
			int difference = n - TypeGroupDynamic::size();
			if(difference > 0) {
				for(unsigned int i = 0; i < difference; i++) {
					// register member
					if (not memberStateIndexs.empty()) {
						memberIndexs.push_back(memberIndexs.back() + 1);
						memberStateIndexs.push_back(memberStateIndexs.back() + memberStateSizes.back());
					} else {
						memberIndexs.push_back(0);
						memberStateIndexs.push_back(0);
					}
					sMemberSteps.push_back(std::make_shared<TypeMemberStep>(*sMemberStep));
					memberStateSizes.push_back(sMemberStep->stateSize());
				}
			} else if (difference < 0) {
				for(unsigned int i = 0; i < std::abs(difference); i++) {
					// unregister member
					memberIndexs.pop_back();
					sMemberSteps.pop_back();
					memberStateIndexs.pop_back();
					memberStateSizes.pop_back();
				}
			}
		}
	public:
		using TypeGroupDynamic::memberIndexs;
		using TypeGroupDynamic::sMemberSteps;
		using TypeGroupDynamic::memberStateSizes;
		using TypeGroupDynamic::memberStateIndexs;
	public:
		std::shared_ptr<TypeMemberStep> sMemberStep;
};

}

#endif
