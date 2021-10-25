#ifndef SL0_POINT_H
#define SL0_POINT_H
#pragma once

// std includes
#include <memory> // shared_ptr
#include <vector>
// module includes
#include "sl0/object.h"

namespace sl0 {

template<template<int> typename TypeVector, unsigned int DIM, template<typename...> class TypeRef, template<typename...> class TypeView, typename TypeFlow>
class StepPoint : public StepObjectStatic<TypeVector, DIM, TypeRef, DIM> {
    public:
        using StepObjectStatic<TypeVector, DIM, TypeRef, DIM>::StateSize;
        using typename StepObjectStatic<TypeVector, DIM, TypeRef, StateSize>::TypeStateStatic;
        using typename StepObjectStatic<TypeVector, DIM, TypeRef, StateSize>::TypeSpaceVector;
    public:
        StepPoint(const std::shared_ptr<TypeFlow>& p_sFlow) : sFlow(p_sFlow) {

        }

        TypeStateStatic operator()(const TypeRef<const TypeStateStatic>& state, const double& t) const override {
            TypeStateStatic dState;
            TypeView<TypeSpaceVector> dX = x(dState);
            dX = sFlow->getVelocity(cX(state), t);
            return dState;
        }
    public:
        TypeView<const TypeSpaceVector> cX(const TypeRef<const TypeStateStatic>& state) const {
            return TypeView<const TypeSpaceVector>(state.data());
        }
        TypeView<TypeSpaceVector> x(TypeRef<TypeStateStatic> state) const {
            return TypeView<TypeSpaceVector>(state.data());
        }
    public:
        std::vector<TypeSpaceVector> positions(const TypeRef<const TypeStateStatic>& state) const override {
            return { cX(state) };
        }
    public:
        std::shared_ptr<TypeFlow> sFlow;
};

template<template<int> typename TypeVector, unsigned int DIM, template<typename...> class TypeRef, template<typename...> class TypeView, typename TypeFlow, typename TypeSolver>
class Point : public ObjectStatic<TypeRef, TypeSolver, StepPoint<TypeVector, DIM, TypeRef, TypeView, TypeFlow>> {
    public:
        using TypeStep = StepPoint<TypeVector, DIM, TypeRef, TypeView, TypeFlow>;
    public:
        Point(const std::shared_ptr<TypeFlow>& sFlow) : ObjectStatic<TypeRef, TypeSolver, TypeStep>::ObjectStatic(std::make_shared<TypeStep>(sFlow)) {
        }
    public:
        using ObjectStatic<TypeRef, TypeSolver, TypeStep>::sSolver;
        using ObjectStatic<TypeRef, TypeSolver, TypeStep>::sStep;
        using ObjectStatic<TypeRef, TypeSolver, TypeStep>::state;
        using ObjectStatic<TypeRef, TypeSolver, TypeStep>::t;
};

}

#endif
