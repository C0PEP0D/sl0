#ifndef SL0_CHAIN_STATIC_H
#define SL0_CHAIN_STATIC_H
#pragma once

// std includes
#include <memory> // shared_ptr
#include <vector>
#include <numeric> // iota
#include <execution>
// module includes
#include "sl0/group/static.h"
#include "p0l/interpolation.h"
#include "m0sh/non_uniform.h"
#include "m0sh/structured_sub.h"

namespace sl0 {

template<template<int> typename TypeVector, unsigned int DIM, template<typename...> class TypeRef, unsigned int Size, template<typename...> class TypeView, typename TypeMemberStep>
class StepChainStatic : public StepGroupStatic<TypeVector, DIM, TypeRef, Size, TypeView, TypeMemberStep> {
    public:
        using TypeStepGroupStatic = StepGroupStatic<TypeVector, DIM, TypeRef, Size, TypeView, TypeMemberStep>;
        using typename TypeStepGroupStatic::TypeStateDynamic;
    public:
        template<typename ...Args>
        using TypeContainer = std::vector<Args...>;
        using TypeMesh = m0sh::NonUniform<TypeVector<1>, TypeRef, TypeContainer>;
        using TypeMeshSub = m0sh::StructuredSub<TypeVector<1>, TypeRef, TypeContainer>;
    public:
        StepChainStatic(const std::shared_ptr<TypeMemberStep>& p_sMemberStep, const double& p_length, const unsigned int& p_interpolationOrder) : TypeStepGroupStatic(p_sMemberStep), length(p_length), interpolationOrder(p_interpolationOrder), coordinates(DIM, TypeContainer<double>(Size)) {
        }

        void update(TypeRef<TypeStateDynamic> state, const double& t) override {
            TypeStepGroupStatic::update(state, t);
            // coordinates
            for(std::size_t i = 0; i < Size; i++) {
                const TypeVector<DIM> position = sMemberStep->cX(cMemberState(state, i));
                for(std::size_t j = 0; j < DIM; j++) {
                    coordinates[j][i] = position[j];
                }
            }
            // mesh
            double l = 0.0;
            TypeContainer<double> gridPoints(Size);
            gridPoints[0] = 0.0;
            for(std::size_t i = 1; i < Size; i++){ 
                const auto x1 = sMemberStep->cX(cMemberState(state, i));
                const auto x0 = sMemberStep->cX(cMemberState(state, i-1));
                l += (x1 - x0).norm();
                gridPoints[i] = l;
            }
            for(std::size_t i = 1; i < Size; i++){ 
                gridPoints[i] /= l;
            }
            sMesh = std::make_shared<TypeMesh>(gridPoints, false);
            // interpolation
            const double sStart = (0.5 * l - 0.5 * length) / l;
            const double ds = length/l/(Size-1);
            for(std::size_t i = 0; i < Size; i++) {
                const TypeVector<1> s(sStart + i * ds);
                TypeView<TypeVector<DIM>> x = sMemberStep->x(memberState(state, i));
                for(std::size_t j = 0; j < DIM; j++) {
                    x[j] = p0l::lagrangeMeshPoint<TypeMesh, TypeContainer, double, TypeVector<1>, TypeRef, TypeMeshSub>(sMesh, coordinates[j], s, interpolationOrder + 1, false);
                }
            }
        }
    public:
        TypeVector<DIM> cX(const TypeRef<const TypeStateDynamic> state, const double& p_s) const {
            const TypeVector<1> s(p_s);
            // coordinates
            TypeContainer<TypeContainer<double>> tmpCoordinates(DIM, TypeContainer<double>(Size)); // TODO: do this better
            for(std::size_t i = 0; i < Size; i++) {
                const TypeVector<DIM> position = sMemberStep->cX(cMemberState(state, i));
                for(std::size_t j = 0; j < DIM; j++) {
                    tmpCoordinates[j][i] = position[j];
                }
            }
            // mesh
            double l = 0.0;
            TypeContainer<double> gridPoints(Size);
            gridPoints[0] = 0.0;
            for(std::size_t i = 1; i < Size; i++){ 
                const auto x1 = sMemberStep->cX(cMemberState(state, i));
                const auto x0 = sMemberStep->cX(cMemberState(state, i-1));
                l += (x1 - x0).norm();
                gridPoints[i] = l;
            }
            for(std::size_t i = 1; i < Size; i++){ 
                gridPoints[i] /= l;
            }
            sMesh = std::make_shared<TypeMesh>(gridPoints, false);
            // interpolate for each coordinates
            TypeVector<DIM> x;
            for(std::size_t j = 0; j < DIM; j++) {
                x[j] = p0l::lagrangeMeshPoint<TypeMesh, TypeContainer, double, TypeVector<1>, TypeRef, TypeMeshSub>(sMesh, tmpCoordinates[j], s, interpolationOrder + 1, false);
            }
            return x;
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
        // shape description
        std::shared_ptr<TypeMesh> sMesh;
        std::vector<std::vector<double>> coordinates; // TODO: iterate avec all state, not just coordinates ?
    public:
        using TypeStepGroupStatic::cMemberState;
        using TypeStepGroupStatic::memberState;
    public:
        using TypeStepGroupStatic::sMemberStep;
};

template<template<int> typename TypeVector, unsigned int DIM, template<typename...> class TypeRef, unsigned int _Size, template<typename...> class TypeView, typename TypeMemberStep, typename TypeSolver>
class ChainStatic : public ObjectDynamic<TypeSolver, StepChainStatic<TypeVector, DIM, TypeRef, _Size, TypeView, TypeMemberStep>> {
    public:
        using TypeStep = StepChainStatic<TypeVector, DIM, TypeRef, _Size, TypeView, TypeMemberStep>;
    public:
        ChainStatic(const std::shared_ptr<TypeMemberStep>& p_sMemberStep, const double& p_length, const unsigned int& p_interpolationOrder) : ObjectDynamic<TypeSolver, TypeStep>(std::make_shared<TypeStep>(p_sMemberStep, p_length, p_interpolationOrder)) {
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
