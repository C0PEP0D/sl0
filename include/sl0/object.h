#ifndef SL0_OBJECT_H
#define SL0_OBJECT_H
#pragma once

// std includes
#include <memory> // shared_ptr
#include <vector>

namespace sl0 {

template<template<int> typename TypeVector, unsigned int _DIM>
class StepObject {
    public:
        static const unsigned int DIM = _DIM;
        using TypeSpaceVector = TypeVector<DIM>;
        using TypeStateVectorDynamic = TypeVector<-1>;
    public:
        StepObject() {
        }
        // Returns dState
        virtual TypeStateVectorDynamic operator()(const double* state, const double& t) const = 0;
        // Applies non linear changes to state
        virtual void update(std::vector<double>& state, const double& t) {
        };
        // Returns current state size
        virtual unsigned int stateSize() const = 0;
        // Returns object's positions
        virtual std::vector<TypeSpaceVector> positions(const double* pState) const = 0;
};

template<template<int> typename TypeVector, unsigned int DIM>
class StepObjectStaticBase : public StepObject<TypeVector, DIM> {
    public:
        using Type = StepObject<TypeVector, DIM>;
        using typename Type::TypeSpaceVector;
        using typename Type::TypeStateVectorDynamic;
    public:
        StepObjectStaticBase() {
        }
        // Applies non linear changes to state
        void update(std::vector<double>& state, const double& t) override {
            update(state.data(), t);
        }
        virtual void update(double* pState, const double& t) {
        };
};

template<template<int> typename TypeVector, unsigned int DIM, unsigned int _StateSize>
class StepObjectStatic : public StepObjectStaticBase<TypeVector, DIM> {
    public:
        using Type = StepObject<TypeVector, DIM>;
        using typename Type::TypeSpaceVector;
        using typename Type::TypeStateVectorDynamic;
        static const unsigned int StateSize = _StateSize;
        //using TypeStateVector = TypeVector<StateSize>;
    public:
        StepObjectStatic() {
        }
        // Applies non linear changes to state
        void update(std::vector<double>& state, const double& t) override {
            update(state.data(), t);
        }
        virtual void update(double* pState, const double& t) {
        };
        // Returns current state size
        unsigned int stateSize() const override {
            return StateSize;
        }
};

template<typename TypeSolver, typename _TypeStep>
class ObjectDynamic {
    public:
        using TypeStep = _TypeStep;
    public:
        ObjectDynamic(const std::shared_ptr<TypeStep>& p_sStep) : sSolver(std::make_shared<TypeSolver>()), sStep(p_sStep), t(0.0) {
        }
        virtual void update(const double& dt) {
            (*sSolver)(*sStep, state.data(), state.size(), t, dt);
            (*sStep).update(state, t);
            t += dt;
        }
    public:
        // Parameters
        std::shared_ptr<TypeSolver> sSolver;
        std::shared_ptr<TypeStep> sStep;
        // State
        double t;
        std::vector<double> state;
};


template<typename TypeSolver, typename _TypeStep>
class ObjectStatic {
    public:
        using TypeStep = _TypeStep;
    public:
        ObjectStatic(const std::shared_ptr<TypeStep>& p_sStep) : sSolver(std::make_shared<TypeSolver>()), sStep(p_sStep), t(0.0), state(sStep->stateSize(), 0.0) {
        }
        virtual void update(const double& dt) {
            (*sSolver)(*sStep, state.data(), state.size(), t, dt);
            (*sStep).update(state.data(), t);
            t += dt;
        }
    public:
        // Parameters
        std::shared_ptr<TypeSolver> sSolver;
        std::shared_ptr<TypeStep> sStep;
        // State
        double t;
        std::vector<double> state;
};

}

#endif
