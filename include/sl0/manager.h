#ifndef SL0_MANAGER_H
#define SL0_MANAGER_H
#pragma once

// std includes
#include "sl0/object.h"

namespace sl0 {

template<template<int> typename TypeVector, unsigned int _DIM>
class StepManager {
    public:
        static const unsigned int DIM = _DIM;
        using TypeSpaceVector = TypeVector<DIM>;
        using TypeStateVectorDynamic = TypeVector<-1>;
    public:
        StepManager() {
        }
        // Returns dState
        virtual TypeStateVectorDynamic operator()(const double* pState, const double& t, const unsigned int index) const = 0;
        // Applies non linear changes to state
        virtual void update(std::vector<std::vector<double>>& states, const double& t) {
        };
        // Returns object's positions
        // virtual std::vector<TypeSpaceVector> positions(const double* pState) const = 0;
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

}

#endif
