#ifndef SL0_POINT_H
#define SL0_POINT_H
#pragma once

// std includes
#include <memory> // shared_ptr
#include <vector>
// module includes
#include "sl0/object.h"

namespace sl0 {

template<template<int> typename TypeVector, unsigned int DIM, template<typename...> class TypeView, typename TypeFlow>
class StepPoint : public StepObjectStatic<TypeVector, DIM, DIM> {
    public:
        using Type = StepObjectStatic<TypeVector, DIM, DIM>;
        using Type::StateSize;
        using typename Type::TypeSpaceVector;
        using typename Type::TypeStateVectorDynamic;
    public:
        StepPoint(const std::shared_ptr<TypeFlow>& p_sFlow) : sFlow(p_sFlow) {

        }

        TypeStateVectorDynamic operator()(const double* pState, const double& t) const override {
            TypeStateVectorDynamic dState(Type::stateSize());
            TypeView<TypeSpaceVector> dX = x(dState.data());
            dX = sFlow->getVelocity(cX(pState), t);
            return dState;
        }
    public:
        TypeView<const TypeSpaceVector> cX(const double* pState) const {
            return TypeView<const TypeSpaceVector>(pState);
        }
        TypeView<TypeSpaceVector> x(double* pState) const {
            return TypeView<TypeSpaceVector>(pState);
        }
    public:
        std::vector<TypeSpaceVector> positions(const double* pState) const override {
            return { cX(pState) };
        }
    public:
        std::shared_ptr<TypeFlow> sFlow;
};

template<template<int> typename TypeVector, unsigned int DIM, template<typename...> class TypeView, typename TypeFlow, typename TypeSolver>
class Point : public ObjectStatic<TypeSolver, StepPoint<TypeVector, DIM, TypeView, TypeFlow>> {
    public:
        using TypeStep = StepPoint<TypeVector, DIM, TypeView, TypeFlow>;
        using Type = ObjectStatic<TypeSolver, TypeStep>;
    public:
        Point(const std::shared_ptr<TypeFlow>& sFlow) : ObjectStatic<TypeSolver, TypeStep>::ObjectStatic(std::make_shared<TypeStep>(sFlow)) {
        }
    public:
        using Type::sSolver;
        using Type::sStep;
        using Type::state;
        using Type::t;
};

}

#endif
