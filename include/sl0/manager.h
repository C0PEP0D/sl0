#ifndef SL0_MANAGER_H
#define SL0_MANAGER_H
#pragma once

// std includes
#include "sl0/group/dynamic.h"

namespace sl0 {

template<template<int> typename TypeVector, unsigned int _DIM>
class StepManagerBase {
	public:
		static const unsigned int DIM = _DIM;
		using TypeSpaceVector = TypeVector<DIM>;
		using TypeStateVectorDynamic = TypeVector<-1>;
	public:
		StepManagerBase() {
		}
		// Returns dState
		virtual TypeStateVectorDynamic operator()(const double* pState, const double& t, const unsigned int index) const = 0;
		// Applies non linear changes to state
		virtual void update(std::vector<std::vector<double>>& states, const double& t) = 0;
		// Returns object's positions
		virtual std::vector<TypeSpaceVector> positions(const double* pState, const unsigned int index) const = 0;
		// Returns managed objects number
		virtual std::size_t number() const = 0;
	public:
		// managed management
		virtual void removeManaged(std::vector<std::vector<double>>& states, const unsigned int& managedIndex) = 0;
		virtual void removeManaged(std::vector<std::vector<double>>& states) = 0;
};

template<template<int> typename TypeVector, unsigned int _DIM, typename _TypeManagedStep>
class StepManager : public StepManagerBase<TypeVector, _DIM> {
	public:
		using Type = StepManagerBase<TypeVector, _DIM>;
		using Type::DIM;
		using typename Type::TypeSpaceVector;
		using typename Type::TypeStateVectorDynamic;
		using TypeManagedStep = _TypeManagedStep;
	public:
		StepManager() {
		}
		// Returns dState
		TypeStateVectorDynamic operator()(const double* pState, const double& t, const unsigned int index) const override {
			return (*sManagedSteps[index])(pState, t);
		}
		// Applies non linear changes to state
		void update(std::vector<std::vector<double>>& states, const double& t) override {
			for (unsigned int managedIndex = 0; managedIndex < states.size(); managedIndex++) {
				(*sManagedSteps[managedIndex]).update(states[managedIndex], t);
			}
		};
		// Returns object's positions
		std::vector<TypeSpaceVector> positions(const double* pState, const unsigned int index) const override {
			return sManagedSteps[index]->positions(pState);
		}
		// Returns managed objects number
		std::size_t number() const override {
			return sManagedSteps.size();
		}
	public:
		// member management
		void addManaged(std::vector<std::vector<double>>& states, std::shared_ptr<TypeManagedStep> sManagedStep) {
			states.emplace_back(sManagedStep->stateSize());
			// register managed
			if (not managedIndexs.empty()) {
				managedIndexs.push_back(managedIndexs.back() + 1);
			} else {
				managedIndexs.push_back(0);
			}
			sManagedSteps.push_back(sManagedStep);
		}

		void removeManaged(std::vector<std::vector<double>>& states, const unsigned int& managedIndex) override {
			states.erase(states.begin() + managedIndex);
			// unregister managed
			managedIndexs.pop_back();
			sManagedSteps.erase(sManagedSteps.begin() + managedIndex);
		}


		void removeManaged(std::vector<std::vector<double>>& states) override {
			states.pop_back();
			// unregister managed
			managedIndexs.pop_back();
			sManagedSteps.pop_back();
		}
	public:
		std::vector<unsigned int> managedIndexs;
		std::vector<std::shared_ptr<TypeManagedStep>> sManagedSteps;
};

template<typename TypeSolver, typename _TypeStep>
class Manager {
	public:
		using TypeStep = _TypeStep;
	public:
		Manager(const std::shared_ptr<TypeStep>& p_sStep) : sSolver(std::make_shared<TypeSolver>()), sStep(p_sStep), t(0.0) {
		}
		virtual void update(const double& dt) {
			for (unsigned int index = 0; index < states.size(); index++) {
				(*sSolver)([this, index](const double* pState, const double& t) { return (*sStep)(pState, t, index); }, states[index].data(), states[index].size(), t, dt);
			}
			(*sStep).update(states, t);
			t += dt;
		}
	public:
		// Parameters
		std::shared_ptr<TypeSolver> sSolver;
		std::shared_ptr<TypeStep> sStep;
		// State
		double t;
		std::vector<std::vector<double>> states;
};

template<template<int> typename TypeVector, unsigned int _DIM>
class StepManagerHomogeneousBase : public StepManagerBase<TypeVector, _DIM> {
	public:
		using Type = StepManagerBase<TypeVector, _DIM>;
		using Type::DIM;
		using typename Type::TypeSpaceVector;
		using typename Type::TypeStateVectorDynamic;
	public:
		StepManagerHomogeneousBase() {
		}
	public:
		virtual void addManaged(std::vector<std::vector<double>>& states) = 0;
	public:
		virtual void registerStates(std::vector<std::vector<double>>& state) = 0;
};

template<template<int> typename TypeVector, unsigned int DIM, typename _TypeManagedStep>
class StepManagerHomogeneous : public StepManager<TypeVector, DIM, _TypeManagedStep>, public StepManagerHomogeneousBase<TypeVector, DIM> {
	public:
		using TypeManager = StepManager<TypeVector, DIM, _TypeManagedStep>;
		using TypeManagerHomogeneousBase = StepManagerHomogeneousBase<TypeVector, DIM>;
		using typename TypeManager::TypeSpaceVector;
		using typename TypeManager::TypeStateVectorDynamic;
		using typename TypeManager::TypeManagedStep;
	public:
		StepManagerHomogeneous(std::shared_ptr<TypeManagedStep> p_sManagedStep) : sManagedStep(p_sManagedStep) {
		}
	public:
		// Returns dState
		TypeStateVectorDynamic operator()(const double* pState, const double& t, const unsigned int index) const override {
			return TypeManager::operator()(pState, t, index);
		}
		// Applies non linear changes to state
		void update(std::vector<std::vector<double>>& states, const double& t) override {
			return TypeManager::update(states, t);
		};
		// Returns object's positions
		std::vector<TypeSpaceVector> positions(const double* pState, const unsigned int index) const override {
			return TypeManager::positions(pState, index);
		}
		// Returns managed objects number
		std::size_t number() const override {
			return TypeManager::number();
		}
	public:
		// member management
		void addManaged(std::vector<std::vector<double>>& states, std::shared_ptr<TypeManagedStep> sManagedStep) {
			TypeManager::addManaged(states, sManagedStep);
		}
		void removeManaged(std::vector<std::vector<double>>& states, const unsigned int& managedIndex) override {
			TypeManager::removeManaged(states, managedIndex);
		}
		void removeManaged(std::vector<std::vector<double>>& states) override {
			TypeManager::removeManaged(states);
		}
	public:
		void addManaged(std::vector<std::vector<double>>& states) override {
			TypeManager::addManaged(states, std::make_shared<TypeManagedStep>(*sManagedStep));
		}
	public:
		void registerStates(std::vector<std::vector<double>>& states) override {
			// register and unregister if states sizes does not correspond to the current number of managed objects
			int difference = states.size() - TypeManager::number();
			if(difference > 0) {
				for(unsigned int i = 0; i < difference; i++) {
					// register managed
					if (not managedIndexs.empty()) {
						managedIndexs.push_back(managedIndexs.back() + 1);
					} else {
						managedIndexs.push_back(0);
					}
					sManagedSteps.push_back(std::make_shared<TypeManagedStep>(*sManagedStep));
				}
			} else if (difference < 0) {
				for(unsigned int i = 0; i < std::abs(difference); i++) {
					// unregister member
					sManagedSteps.pop_back();
					managedIndexs.pop_back();
				}
			}
			// register managed states
			for(unsigned int i = 0; i < TypeManager::number(); i++) {
				sManagedSteps[i]->registerState(states[i]);
			}
		}
	public:
		using TypeManager::managedIndexs;
		using TypeManager::sManagedSteps;
	public:
		std::shared_ptr<TypeManagedStep> sManagedStep;
};

}

#endif
