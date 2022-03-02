#ifndef SL0_CHAIN_MANAGER_H
#define SL0_CHAIN_MANAGER_H
#pragma once

// module includes
#include "sl0/manager.h"
#include "sl0/chain/dynamic.h"

namespace sl0 {

template<template<int> typename TypeVector, unsigned int DIM, template<typename...> class TypeView, template<typename...> class TypeRef, typename TypeMemberStep>
class StepChainManager : public StepManager<TypeVector, DIM> {
	public:
		using Type = StepManager<TypeVector, DIM>;
		using typename Type::TypeSpaceVector;
		using typename Type::TypeStateVectorDynamic;
	public:
		using TypeStepChainDynamic = StepChainDynamic<TypeVector, DIM, TypeView, TypeRef, TypeMemberStep>;
	public:
		StepChainManager() {
		}

		TypeStateVectorDynamic operator()(const double* pState, const double& t, const unsigned int index) const override {
			return (*sChainSteps[index])(pState, t);
		}

		void update(std::vector<std::vector<double>>& states, const double& t) override {
			// Solve intersections
			// TODO
			// Update all chains
			for (unsigned int index = 0; index < states.size(); index++) {
            	(*sChainSteps[index]).update(states[index], t);
            }
		}
	public:
		// member management
		void addChain(std::vector<std::vector<double>>& states, std::shared_ptr<TypeStepChainDynamic> sChainStep) {
			states.emplace_back(sChainStep->stateSize());
			sChainSteps.push_back(sChainStep);
		}

		void removeChain(std::vector<std::vector<double>>& states, const unsigned int& chainIndex) {
			states.erase(states.begin() + chainIndex);
			sChainSteps.erase(sChainSteps.begin() + chainIndex);
		}


		void removeChain(std::vector<std::vector<double>>& states) {
			states.pop_back();
			sChainSteps.pop_back();
		}
	public:
		std::vector<std::shared_ptr<TypeStepChainDynamic>> sChainSteps;
};

template<template<int> typename TypeVector, unsigned int DIM, template<typename...> class TypeView, template<typename...> class TypeRef, typename TypeMemberStep, typename TypeSolver>
class ChainManager : public Manager<TypeSolver, StepChainManager<TypeVector, DIM, TypeView, TypeRef, TypeMemberStep>> {
	public:
		using TypeStep = StepChainManager<TypeVector, DIM, TypeView, TypeRef, TypeMemberStep>;
	public:
		ChainManager() : Manager<TypeSolver, TypeStep>(std::make_shared<TypeStep>()) {
		}
	public:
		// Inherited
		using Manager<TypeSolver, TypeStep>::sSolver;
		using Manager<TypeSolver, TypeStep>::sStep;
		using Manager<TypeSolver, TypeStep>::states;
		using Manager<TypeSolver, TypeStep>::t;
};

}

#endif
