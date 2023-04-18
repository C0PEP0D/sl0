#ifndef SL0_CHAIN_MANAGER_H
#define SL0_CHAIN_MANAGER_H
#pragma once

// module includes
#include "sl0/manager.h"
#include "sl0/chain/dynamic.h"

namespace sl0 {

template<template<int> typename TypeVector, unsigned int DIM, template<typename...> class TypeView, template<typename...> class TypeRef, typename TypeMemberStep>
class StepChainManager : public StepManagerHomogeneous<TypeVector, DIM, StepChainDynamic<TypeVector, DIM, TypeView, TypeRef, TypeMemberStep>> {
	public:
		using Type = StepManagerHomogeneous<TypeVector, DIM, StepChainDynamic<TypeVector, DIM, TypeView, TypeRef, TypeMemberStep>>;
		using typename Type::TypeSpaceVector;
		using typename Type::TypeStateVectorDynamic;
		using typename Type::TypeManagedStep;
	public:
		StepChainManager(std::shared_ptr<TypeMemberStep> sMemberStep, const double& dl, const unsigned int& interpolationOrder, const bool closed) : Type(std::make_shared<TypeManagedStep>(sMemberStep, dl, interpolationOrder, closed)) {
		}

		void update(std::vector<std::vector<double>>& states, const double& t) override {
			Type::update(states, t);
			// solve intersections
			std::vector<bool> isChainInGroupB(states.size(), false); // Nice but not enough to segregate stuff
			// TODO: DEAL WITH j + 1 = 0
			for (unsigned int chainIndex = 0; chainIndex < states.size(); chainIndex++) {
				// std::cout << "intersection nb: " << (*sManagedSteps[chainIndex]).intersections(states[chainIndex].data()).size() << std::endl;
            	const auto firstIntersectionData = (*sManagedSteps[chainIndex]).firstIntersection(states[chainIndex].data());
            	if(firstIntersectionData.j != 0 && firstIntersectionData.i != 0) {
            		// debug
            		const TypeStateVectorDynamic intersectionState = 0.5 * (sManagedSteps[chainIndex]->cState(states[chainIndex].data(), firstIntersectionData.iS) + sManagedSteps[chainIndex]->cState(states[chainIndex].data(), firstIntersectionData.jS));
            		// add new chain
            		states.emplace_back(states[chainIndex].begin() + (firstIntersectionData.i + 1) * TypeManagedStep::TypeMemberStep::StateSize, states[chainIndex].begin() + (firstIntersectionData.j + 2) * TypeManagedStep::TypeMemberStep::StateSize);
            		isChainInGroupB.emplace_back(!isChainInGroupB[chainIndex]);
            		states[chainIndex].erase(states[chainIndex].begin() + (firstIntersectionData.i + 1) * TypeManagedStep::TypeMemberStep::StateSize, states[chainIndex].begin() + (firstIntersectionData.j + 1) * TypeManagedStep::TypeMemberStep::StateSize);
            		Type::registerStates(states);
            		const unsigned int newChainIndex = states.size() - 1;
            		// set intersections states
            		TypeView<TypeStateVectorDynamic>(sManagedSteps[newChainIndex]->memberState(states[newChainIndex].data(), sManagedSteps[newChainIndex]->size() - 1), intersectionState.size()) = intersectionState;
					TypeView<TypeStateVectorDynamic>(sManagedSteps[chainIndex]->memberState(states[chainIndex].data(), firstIntersectionData.i + 1), intersectionState.size()) = intersectionState;
            	}
            }
            // compute group length
            double lengthGroupA = 0.0;
            double lengthGroupB = 0.0;
            for (unsigned int chainIndex = 0; chainIndex < states.size(); chainIndex++) {
            	if (isChainInGroupB[chainIndex]) {
            		lengthGroupB += (*sManagedSteps[chainIndex]).length(states[chainIndex].data());
            	} else {
            		lengthGroupA += (*sManagedSteps[chainIndex]).length(states[chainIndex].data());
            	}
            }
            bool isGroupBVirtual = lengthGroupB < lengthGroupA;
            // remove all virtual chains and all states that are too small
            for (int chainIndex = states.size() - 1; chainIndex > -1; chainIndex--) {
				if(states[chainIndex].size() / TypeManagedStep::TypeMemberStep::StateSize < Type::sManagedStep->interpolationOrder + 1) {
					states.erase(states.begin() + chainIndex);
				} else {
					if (isGroupBVirtual) {
						if (isChainInGroupB[chainIndex]) {
							states.erase(states.begin() + chainIndex);
						}
					} else {
						if (not isChainInGroupB[chainIndex]) {
							states.erase(states.begin() + chainIndex);
						}
					}
				}
            }
            Type::registerStates(states);
		}
	public:
		using Type::sManagedSteps;
};

template<template<int> typename TypeVector, unsigned int DIM, template<typename...> class TypeView, template<typename...> class TypeRef, typename TypeMemberStep, typename TypeSolver>
class ChainManager : public Manager<TypeSolver, StepChainManager<TypeVector, DIM, TypeView, TypeRef, TypeMemberStep>> {
	public:
		using TypeStep = StepChainManager<TypeVector, DIM, TypeView, TypeRef, TypeMemberStep>;
	public:
		ChainManager(std::shared_ptr<TypeMemberStep> sMemberStep, const double& dl, const unsigned int& interpolationOrder) : Manager<TypeSolver, TypeStep>(std::make_shared<TypeStep>(sMemberStep, dl, interpolationOrder)) {
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
