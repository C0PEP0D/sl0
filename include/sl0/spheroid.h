#ifndef SL0_SPHEROID_H
#define SL0_SPHEROID_H
#pragma once

#include <memory>
#include <cmath>
#include <vector>

#include "sl0/object.h"

namespace sl0 {

template<template<int> typename TypeVector, unsigned int DIM, template<typename...> class TypeRef, template<typename...> class TypeView>
class StepAxis : public StepObjectStatic<TypeVector, DIM, TypeRef, 2*DIM> {
    public:
        using StepObjectStatic<TypeVector, DIM, TypeRef, 2*DIM>::StateSize;
        using typename StepObjectStatic<TypeVector, DIM, TypeRef, StateSize>::TypeStateStatic;
        using typename StepObjectStatic<TypeVector, DIM, TypeRef, StateSize>::TypeSpaceVector;
    public:
        using StepObjectStatic<TypeVector, DIM, TypeRef, StateSize>::StepObjectStatic;
        void update(TypeRef<TypeStateStatic> state, const double& t) override {
            StepObjectStatic<TypeVector, DIM, TypeRef, StateSize>::update(state, t);
            // just orthonormalize axis
            normalizeAxis(state);
        }
        // sub update methods
        void normalizeAxis(TypeRef<TypeStateStatic> state) {
            axis(state).normalize();
        }
    public:
        virtual TypeView<const TypeSpaceVector> cAxis(const TypeRef<const TypeStateStatic>& state) const = 0;
        virtual TypeView<TypeSpaceVector> axis(TypeRef<TypeStateStatic> state) const = 0;
};

template<template<int> typename TypeVector, unsigned int DIM, template<typename...> class TypeRef, template<typename...> class TypeView, typename TypeFlow>
class StepSpheroid : public StepAxis<TypeVector, DIM, TypeRef, TypeView> {
    public:
        using TypeStepAxis = StepAxis<TypeVector, DIM, TypeRef, TypeView>;
        using TypeStepAxis::StateSize;
        using typename TypeStepAxis::TypeStateStatic;
        using typename TypeStepAxis::TypeSpaceVector;
    public:
        StepSpheroid(const std::shared_ptr<TypeFlow>& p_sFlow, const double& prop) : TypeStepAxis(), sFlow(p_sFlow) {
            setProportion(prop);
        }
        TypeStateStatic operator()(const TypeRef<const TypeStateStatic>& state, const double& t) const override {
            TypeStateStatic dState = TypeStateStatic::Zero();
            // Get data from state
            const TypeView<const TypeSpaceVector> sX = cX(state);
            const TypeView<const TypeSpaceVector> sAxis = cAxis(state);
            // Prepare to set data
            TypeView<TypeSpaceVector> dX = x(dState);
            TypeView<TypeSpaceVector> dAxis = axis(dState);
            // Compute linear velocity
            dX = sFlow->getVelocity(sX, t);
            // Compute rotaton velocity
            TypeSpaceVector omega = sFlow->getVorticity(sX, t) + factor * (sAxis.cross(sFlow->getStrain(sX, t) * sAxis));
            dAxis = omega.cross(sAxis);
            // Return 
            return dState;
        }
    public:
        TypeView<const TypeSpaceVector> cX(const TypeRef<const TypeStateStatic>& state) const {
            return TypeView<const TypeSpaceVector>(state.data());
        }
        TypeView<TypeSpaceVector> x(TypeRef<TypeStateStatic> state) const {
            return TypeView<TypeSpaceVector>(state.data());
        }
        TypeView<const TypeSpaceVector> cAxis(const TypeRef<const TypeStateStatic>& state) const override {
            return TypeView<const TypeSpaceVector>(state.data() + DIM);
        }
        TypeView<TypeSpaceVector> axis(TypeRef<TypeStateStatic> state) const override {
            return TypeView<TypeSpaceVector>(state.data() + DIM);
        }
    public:
        std::vector<TypeSpaceVector> positions(const TypeRef<const TypeStateStatic>& state) const override {
            return { cX(state) };
        }
    public:
        void setProportion(const double& prop) {
            factor = (std::pow(prop, 2) - std::pow(1.0, 2)) / (std::pow(prop, 2) + std::pow(1.0, 2));
        }
        double getFactor() const {
            return factor;
        }
    public:
        std::shared_ptr<TypeFlow> sFlow;
        double factor;
};

template<template<int> typename TypeVector, unsigned int DIM, template<typename...> class TypeRef, template<typename...> class TypeView, typename TypeFlow, typename TypeSolver>
class Spheroid : public ObjectStatic<TypeRef, TypeSolver, StepSpheroid<TypeVector, DIM, TypeRef, TypeView, TypeFlow>> {
    public:
        using TypeStep = StepSpheroid<TypeVector, DIM, TypeRef, TypeView, TypeFlow>;
    public:
        Spheroid(const std::shared_ptr<TypeFlow>& sFlow, const double& prop) : ObjectStatic<TypeRef, TypeSolver, TypeStep>::ObjectStatic(std::make_shared<TypeStep>(sFlow, prop)) {
        }
    public:
        using ObjectStatic<TypeRef, TypeSolver, TypeStep>::sSolver;
        using ObjectStatic<TypeRef, TypeSolver, TypeStep>::sStep;
        using ObjectStatic<TypeRef, TypeSolver, TypeStep>::state;
        using ObjectStatic<TypeRef, TypeSolver, TypeStep>::t;
};

}

#endif
