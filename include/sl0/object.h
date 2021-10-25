#ifndef SL0_OBJECT_H
#define SL0_OBJECT_H
#pragma once

// std includes
#include <memory> // shared_ptr
#include <vector>

namespace sl0 {

template<template<int> typename TypeVector, unsigned int _DIM, template<typename...> class TypeRef>
class StepObject {
    public:
        using TypeStateDynamic = TypeVector<Eigen::Dynamic>;
        static const unsigned int DIM = _DIM;
        using TypeSpaceVector = TypeVector<DIM>;
    public:
        StepObject() {
        }
        // Returns dState
        virtual TypeStateDynamic operator()(const TypeRef<const TypeStateDynamic>& state, const double& t) const = 0;
        // Applies non linear changes to state
        virtual void update(TypeRef<TypeStateDynamic> state, const double& t) = 0;
        // Returns current state size
        virtual unsigned int stateSize() const = 0;
        // Returns object's positions
        virtual std::vector<TypeSpaceVector> positions(const TypeRef<const TypeStateDynamic>& state) const = 0;
};

template<template<int> typename TypeVector, unsigned int DIM, template<typename...> class TypeRef, unsigned int _StateSize>
class StepObjectStatic : public StepObject<TypeVector, DIM, TypeRef> {
    public:
        using typename StepObject<TypeVector, DIM, TypeRef>::TypeStateDynamic;
        using typename StepObject<TypeVector, DIM, TypeRef>::TypeSpaceVector;
        static const unsigned int StateSize = _StateSize;
        using TypeStateStatic = TypeVector<StateSize>;
    public:
        StepObjectStatic() {
        }
        // Returns dState
        TypeStateDynamic operator()(const TypeRef<const TypeStateDynamic>& state, const double& t) const override {
            return operator()(TypeRef<const TypeStateStatic>(state), t);
        }
        virtual TypeStateStatic operator()(const TypeRef<const TypeStateStatic>& state, const double& t) const = 0;
        // Applies non linear changes to state
        void update(TypeRef<TypeStateDynamic> state, const double& t) override {
            update(TypeRef<TypeStateStatic>(state), t);
        }
        virtual void update(TypeRef<TypeStateStatic> state, const double& t) {}
        // Returns current state size
        unsigned int stateSize() const override {
            return StateSize;
        }
        // Returns object's positions
        std::vector<TypeSpaceVector> positions(const TypeRef<const TypeStateDynamic>& state) const override {
            return positions(TypeRef<const TypeStateStatic>(state));
        }
        virtual std::vector<TypeSpaceVector> positions(const TypeRef<const TypeStateStatic>& state) const = 0;
};

template<typename TypeSolver, typename _TypeStep>
class ObjectDynamic {
    public:
        using TypeStep = _TypeStep;
    public:
        ObjectDynamic(const std::shared_ptr<TypeStep>& p_sStep) : sSolver(std::make_shared<TypeSolver>()), sStep(p_sStep), t(0.0) {
        }
        virtual void update(const double& dt) {
            state = (*sSolver)(*sStep, state, t, dt);
            (*sStep).update(state, t);
            t += dt;
        }
    public:
        // Parameters
        std::shared_ptr<TypeSolver> sSolver;
        std::shared_ptr<TypeStep> sStep;
        // State
        double t;
        typename TypeStep::TypeStateDynamic state;
};


template<template<typename...> class TypeRef, typename TypeSolver, typename _TypeStep>
class ObjectStatic {
    public:
        using TypeStep = _TypeStep;
    public:
        ObjectStatic(const std::shared_ptr<TypeStep>& p_sStep) : sSolver(std::make_shared<TypeSolver>()), sStep(p_sStep), t(0.0), state(TypeStep::TypeStateStatic::Zero()) {
        }
        virtual void update(const double& dt) {
            state = (*sSolver)(*sStep, state, t, dt);
            (*sStep).update(TypeRef<typename TypeStep::TypeStateStatic>(state), t);
            t += dt;
        }
    public:
        // Parameters
        std::shared_ptr<TypeSolver> sSolver;
        std::shared_ptr<TypeStep> sStep;
        // State
        double t;
        typename TypeStep::TypeStateStatic state;
};

}

#endif
