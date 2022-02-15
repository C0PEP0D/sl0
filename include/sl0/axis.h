#ifndef SL0_AXIS_H
#define SL0_AXIS_H
#pragma once

#include <memory>
#include <cmath>
#include <vector>

#include "sl0/object.h"

namespace sl0 {

template<template<int> typename TypeVector, unsigned int DIM, template<typename...> class TypeView>
class StepAxis : public StepObjectStatic<TypeVector, DIM, 2*DIM> {
    public:
        using Type = StepObjectStatic<TypeVector, DIM, 2*DIM>;
        using Type::StateSize;
        using typename Type::TypeSpaceVector;
        using typename Type::TypeStateVectorDynamic;
    public:
        using StepObjectStatic<TypeVector, DIM, StateSize>::StepObjectStatic;
        void update(double* pState, const double& t) override {
            StepObjectStatic<TypeVector, DIM, StateSize>::update(pState, t);
            // just orthonormalize axis
            normalizeAxis(pState);
        }
        // sub update methods
        void normalizeAxis(double* pState) {
            axis(pState).normalize();
        }
    public:
        virtual TypeView<const TypeSpaceVector> cAxis(const double* pState) const = 0;
        virtual TypeView<TypeSpaceVector> axis(double* pState) const = 0;
};

}

#endif
