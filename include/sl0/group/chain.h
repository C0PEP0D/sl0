#ifndef SL0_GROUP_CHAIN_H
#define SL0_GROUP_CHAIN_H
#pragma once

// std includes
#include <memory> // shared_ptr
#include <vector>
#include <numeric> // iota
#include <execution>
// module includes
#include "sl0/group/static.h"
#include "p0l/interpolation.h"
#include "m0sh/regular.h"
#include "m0sh/structured_sub.h"

namespace sl0 {

template<template<int> typename TypeVector, unsigned int DIM, template<typename...> class TypeRef, unsigned int Size, template<typename...> class TypeView, typename TypeMemberStep>
class StepGroupChain : public StepGroupStatic<TypeVector, DIM, TypeRef, Size, TypeView, TypeMemberStep> {
    public:
        using TypeStepGroupStatic = StepGroupStatic<TypeVector, DIM, TypeRef, Size, TypeView, TypeMemberStep>;
        using typename TypeStepGroupStatic::TypeStateDynamic;
    public:
        template<typename ...Args>
        using TypeContainer = std::vector<Args...>;
        using TypeMesh = m0sh::Regular<TypeVector<1>, TypeRef, TypeContainer>;
        using TypeMeshSub = m0sh::StructuredSub<TypeVector<1>, TypeRef, TypeContainer>;
    public:
        StepGroupChain(const std::shared_ptr<TypeMemberStep>& p_sMemberStep, const double& p_length, const unsigned int& p_interpolationOrder) : TypeStepGroupStatic(p_sMemberStep), length(p_length), interpolationOrder(p_interpolationOrder), sMesh(std::make_shared<TypeMesh>(std::vector<long unsigned int>(1, Size), std::vector<double>(1, Size/(Size - 1.0)), TypeVector<1>(-0.5/(Size - 1.0)))) {
        }

        void update(TypeRef<TypeStateDynamic> state, const double& t) override {
            TypeStepGroupStatic::update(state, t);
            // input
            std::vector<double> x(Size);
            std::vector<double> y(Size);
            std::vector<double> z(Size);
            for(std::size_t i = 0; i < Size; i++) {
                const auto position = sMemberStep->cX(cMemberState(state, i));
                x[i] = position[0];
                y[i] = position[1];
                z[i] = position[2];
            }
            // interpolation
            const double aLength = actualLength(state);
            const double sStart = 0.5 * (aLength - length) / aLength;
            const double ds = length/aLength/(Size-1);
            for(std::size_t i = 0; i < Size; i++) {
                const TypeVector<1> s(sStart + i * ds);
                TypeView<TypeVector<DIM>> uX = sMemberStep->x(memberState(state, i));
                uX[0] = p0l::lagrangeMesh<TypeMesh, TypeContainer, double, TypeVector<1>, TypeRef, TypeMeshSub>(sMesh, x, s, interpolationOrder + 1, false);
                uX[1] = p0l::lagrangeMesh<TypeMesh, TypeContainer, double, TypeVector<1>, TypeRef, TypeMeshSub>(sMesh, y, s, interpolationOrder + 1, false);
                uX[2] = p0l::lagrangeMesh<TypeMesh, TypeContainer, double, TypeVector<1>, TypeRef, TypeMeshSub>(sMesh, z, s, interpolationOrder + 1, false);
            }
        }
    public:
        double actualLength(const TypeRef<const TypeStateDynamic>& state) const {
            double l = 0.0;
            for(std::size_t i = 1; i < Size; i++){ 
                const auto x1 = sMemberStep->cX(cMemberState(state, i));
                const auto x0 = sMemberStep->cX(cMemberState(state, i-1));
                l += (x1 - x0).norm();
            }
            return l;
        }
    public:
        double length;
        unsigned int interpolationOrder;
    public:
        std::shared_ptr<TypeMesh> sMesh;
    public:
        using TypeStepGroupStatic::cMemberState;
        using TypeStepGroupStatic::memberState;
    public:
        using TypeStepGroupStatic::sMemberStep;
};

template<template<int> typename TypeVector, unsigned int DIM, template<typename...> class TypeRef, unsigned int _Size, template<typename...> class TypeView, typename TypeMemberStep, typename TypeSolver>
class GroupChain : public ObjectDynamic<TypeSolver, StepGroupChain<TypeVector, DIM, TypeRef, _Size, TypeView, TypeMemberStep>> {
    public:
        using TypeStep = StepGroupChain<TypeVector, DIM, TypeRef, _Size, TypeView, TypeMemberStep>;
    public:
        GroupChain(const std::shared_ptr<TypeMemberStep>& p_sMemberStep, const double& p_length, const unsigned int& p_interpolationOrder) : ObjectDynamic<TypeSolver, TypeStep>(std::make_shared<TypeStep>(p_sMemberStep, p_length, p_interpolationOrder)) {
            state.resize(TypeStep::Size * TypeStep::MemberStateSize);
        }
    public:
        // Inherited
        using ObjectDynamic<TypeSolver, TypeStep>::sSolver;
        using ObjectDynamic<TypeSolver, TypeStep>::sStep;
        using ObjectDynamic<TypeSolver, TypeStep>::state;
        using ObjectDynamic<TypeSolver, TypeStep>::t;
};

}

#endif
