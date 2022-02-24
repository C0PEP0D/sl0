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
        StepChainManager(std::shared_ptr<TypeMemberStep> p_sMemberStep, const double& p_dl, const unsigned int& p_interpolationOrder) : sStepChainDynamic(std::make_shared<TypeStepChainDynamic>(p_sMemberStep, p_dl, p_interpolationOrder)) {
        }

        TypeStateVectorDynamic operator()(const double* pState, const double& t, const unsigned int index) const override {
        	return (*sStepChainDynamic)(pState, t);
        }

        void update(std::vector<std::vector<double>>& states, const double& t) override {
        }
    public:
        std::shared_ptr<TypeStepChainDynamic> sStepChainDynamic;
};

template<template<int> typename TypeVector, unsigned int DIM, template<typename...> class TypeView, template<typename...> class TypeRef, typename TypeMemberStep, typename TypeSolver>
class ChainManager : public Manager<TypeSolver, StepChainManager<TypeVector, DIM, TypeView, TypeRef, TypeMemberStep>> {
    public:
        using TypeStep = StepChainManager<TypeVector, DIM, TypeView, TypeRef, TypeMemberStep>;
    public:
        ChainDynamic(const std::shared_ptr<TypeMemberStep>& p_sMemberStep, const double& p_dl, const unsigned int& p_interpolationOrder) : Manager<TypeSolver, TypeStep>(std::make_shared<TypeStep>(p_sMemberStep, p_dl, p_interpolationOrder)) {
        }
    public:
        // Inherited
        using ObjectDynamic<TypeSolver, TypeStep>::sSolver;
        using ObjectDynamic<TypeSolver, TypeStep>::sStep;
        using ObjectDynamic<TypeSolver, TypeStep>::states;
        using ObjectDynamic<TypeSolver, TypeStep>::t;
};

}

#endif
