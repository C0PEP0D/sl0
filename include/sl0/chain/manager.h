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
		StepChainManager(std::shared_ptr<TypeMemberStep> sMemberStep, const double& dl, const unsigned int& interpolationOrder) : Type(std::make_shared<TypeManagedStep>(sMemberStep, dl, interpolationOrder)) {
		}

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
			}
			return intersectionsData;
		}

		void update(std::vector<std::vector<double>>& states, const double& t) override {
			Type::update(states, t);
			// solve intersections
			const unsigned int n = states.size();
			for (unsigned int chainIndex = 0; chainIndex < n; chainIndex++) {
				//std::cout << (*sManagedSteps[chainIndex]).intersections(states[chainIndex].data()).size() << std::endl;
            	for(auto intersectionData : (*sManagedSteps[chainIndex]).intersections(states[chainIndex].data())) {
            		//std::cout << intersectionData.point << std::endl;
            		//std::cout << intersectionData.i << std::endl;
            		//std::cout << intersectionData.j << std::endl;
            		//std::cout << intersectionData.iS << std::endl;
            		//std::cout << intersectionData.jS << std::endl;
            		const TypeStateVectorDynamic intersectionState = 0.5 * (sManagedSteps[chainIndex]->cState(states[chainIndex].data(), intersectionData.iS) + sManagedSteps[chainIndex]->cState(states[chainIndex].data(), intersectionData.jS));
            		// add new chain
            		states.emplace_back(states[chainIndex].begin() + intersectionData.i + 1, states[chainIndex].begin() + intersectionData.j + 1);
            		states[chainIndex].erase(states[chainIndex].begin() + intersectionData.i + 1, states[chainIndex].begin() + intersectionData.j + 1);
            		/*const unsigned int newChainIndex = states.size() - 1;
            		sManagedSteps.push_back(std::make_shared<TypeManagedStep>(std::make_shared<TypeMemberStep>(*sManagedSteps[chainIndex]->sMemberStep), sManagedSteps[chainIndex]->dl, sManagedSteps[chainIndex]->interpolationOrder));
            		//// add new chain members
            		sManagedSteps[newChainIndex]->addMember(states[newChainIndex]);
            		TypeView<TypeStateVectorDynamic>(sManagedSteps[newChainIndex]->memberState(states[newChainIndex].data(), 0), intersectionState.size()) = intersectionState;
            		for(unsigned int memberIndex = intersectionData.i + 1; memberIndex < intersectionData.j + 1; memberIndex++) {
            			sManagedSteps[newChainIndex]->addMember(states[newChainIndex]);
            				TypeView<TypeStateVectorDynamic>(sManagedSteps[newChainIndex]->memberState(states[newChainIndex].data(), sManagedSteps[newChainIndex]->size() - 1), intersectionState.size()) = TypeView<const TypeStateVectorDynamic>(sManagedSteps[chainIndex]->cMemberState(states[chainIndex].data(), memberIndex), intersectionState.size());
            		}
            		sManagedSteps[newChainIndex]->addMember(states[newChainIndex]);
            		TypeView<TypeStateVectorDynamic>(sManagedSteps[newChainIndex]->memberState(states[newChainIndex].data(), sManagedSteps[newChainIndex]->size() - 1), intersectionState.size()) = intersectionState; // TODO: should be a periodicly closed line
					// update original chain
					TypeView<TypeStateVectorDynamic>(sManagedSteps[chainIndex]->memberState(states[chainIndex].data(), intersectionData.i + 1), intersectionState.size()) = intersectionState;
					for(unsigned int memberIndex = intersectionData.j; memberIndex > intersectionData.i + 1; memberIndex--) {
            			sManagedSteps[chainIndex]->removeMember(states[newChainIndex], memberIndex);
            		}*/
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
