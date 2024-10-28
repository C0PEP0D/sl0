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
			std::vector<bool> isChainVirtual(states.size(), false); // Nice but not enough to segregate stuff
			// TODO: DEAL WITH j + 1 = 0
			for (unsigned int chainIndex = 0; chainIndex < states.size(); chainIndex++) {
            	const auto firstIntersectionData = (*sManagedSteps[chainIndex]).firstIntersection(states[chainIndex].data());
            	if(firstIntersectionData.j != 0 && firstIntersectionData.i != 0) {
            		// debug
            		const TypeStateVectorDynamic intersectionState = 0.5 * (sManagedSteps[chainIndex]->cState(states[chainIndex].data(), firstIntersectionData.iS) + sManagedSteps[chainIndex]->cState(states[chainIndex].data(), firstIntersectionData.jS));
            		// add new chain
            		states.emplace_back(states[chainIndex].begin() + (firstIntersectionData.i + 1) * TypeManagedStep::TypeMemberStep::StateSize, states[chainIndex].begin() + (firstIntersectionData.j + 2) * TypeManagedStep::TypeMemberStep::StateSize);
            		states[chainIndex].erase(states[chainIndex].begin() + (firstIntersectionData.i + 1) * TypeManagedStep::TypeMemberStep::StateSize, states[chainIndex].begin() + (firstIntersectionData.j + 1) * TypeManagedStep::TypeMemberStep::StateSize);
            		Type::registerStates(states);
            		const unsigned int newChainIndex = states.size() - 1;
            		// set intersections states
            		TypeView<TypeStateVectorDynamic>(sManagedSteps[newChainIndex]->memberState(states[newChainIndex].data(), sManagedSteps[newChainIndex]->size() - 1), intersectionState.size()) = intersectionState;
					TypeView<TypeStateVectorDynamic>(sManagedSteps[chainIndex]->memberState(states[chainIndex].data(), firstIntersectionData.i + 1), intersectionState.size()) = intersectionState;
					// set virtual
					// if ((*sManagedSteps[chainIndex]).length(states[chainIndex].data()) < (*sManagedSteps[newChainIndex]).length(states[newChainIndex].data())) {
						// isChainVirtual[chainIndex] = true;
					// } else {
						// isChainVirtual[chainIndex] = false;
					// }
					// isChainVirtual.emplace_back(!isChainVirtual[chainIndex]);
					isChainVirtual.emplace_back(false);
            	}
            }
            // remove all states that are too small
            for (int chainIndex = states.size() - 1; chainIndex > -1; chainIndex--) {
				if(states[chainIndex].size() / TypeManagedStep::TypeMemberStep::StateSize < Type::sManagedStep->interpolationOrder + 1) {
					states.erase(states.begin() + chainIndex);
					isChainVirtual.pop_back();
				}
            }
            Type::registerStates(states);
            // TESTING: CHECKING WICH CHAIN IS VIRTUAL, TODO: ONLY "NEW ONES" SHOULD BE CHECKED
            for (unsigned int chainIndex = 0; chainIndex < states.size(); chainIndex++) {
            	int intersectionsNumber = 0;
            	// get points
				const double* pChainPointState0 = sManagedSteps[chainIndex]->cMemberState(states[chainIndex].data(), 2);
				const auto& chainPointX0 = sManagedSteps[chainIndex]->sMemberStep->cX(pChainPointState0);
				const double* pChainPointState1 = sManagedSteps[chainIndex]->cMemberState(states[chainIndex].data(), 3);
				const auto& chainPointX1 = sManagedSteps[chainIndex]->sMemberStep->cX(pChainPointState1);
				// compute segment
				const auto& chainSegment = chainPointX1 - chainPointX0;
				const auto& chainSegmentX = 0.5 * (chainPointX0 + chainPointX1);
				const auto& chainSegmentNormal = TypeSpaceVector(chainSegment(1), -chainSegment(0), chainSegment(2));
            	for (unsigned int otherChainIndex = 0; otherChainIndex < states.size(); otherChainIndex++) {
            		if (otherChainIndex != chainIndex) {
            			intersectionsNumber += sManagedSteps[otherChainIndex]->rayIntersections(states[otherChainIndex].data(), chainSegmentX, chainSegmentNormal).size();
            		} else {
            			// intersectionsNumber += sManagedSteps[otherChainIndex]->rayIntersections(states[otherChainIndex].data(), chainSegmentX, chainSegmentNormal, {2}).size();
            		}
            	}
            	isChainVirtual[chainIndex] = (intersectionsNumber % 2 == 1);
            }
            // remove all virtual chains
            for (int chainIndex = states.size() - 1; chainIndex > -1; chainIndex--) {
            	if (isChainVirtual[chainIndex]) {
            		states.erase(states.begin() + chainIndex);
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
