#ifndef SL0_CHAIN_DYNAMIC_H
#define SL0_CHAIN_DYNAMIC_H
#pragma once

// std includes
#include <memory> // shared_ptr
#include <vector>
#include <numeric> // iota
#include <execution>
// module includes
#include "sl0/group/dynamic.h"
#include "p0l/interpolation.h"
#include "m0sh/non_uniform.h"
#include "m0sh/structured_sub.h"

namespace sl0 {

template<template<int> typename TypeVector, unsigned int DIM, template<typename...> class TypeView, template<typename...> class TypeRef, typename TypeMemberStep>
class StepChainDynamic : public StepGroupDynamicHomogeneous<TypeVector, DIM, TypeView, TypeMemberStep> {
    public:
        using Type = StepGroupDynamicHomogeneous<TypeVector, DIM, TypeView, TypeMemberStep>;
        using typename Type::TypeStateVectorDynamic;
    public:
        template<typename ...Args>
        using TypeContainer = std::vector<Args...>;
        using TypeMesh = m0sh::NonUniform<TypeVector<1>, TypeRef, TypeContainer>;
        using TypeMeshSub = m0sh::StructuredSub<TypeVector<1>, TypeRef, TypeContainer>;
    public:
        StepChainDynamic(std::shared_ptr<TypeMemberStep> p_sMemberStep, const double& p_dl, const unsigned int& p_interpolationOrder) : Type(p_sMemberStep), dl(p_dl), interpolationOrder(p_interpolationOrder), coordinates(DIM) {
        }

        void update(std::vector<double>& state, const double& t) override {
            Type::update(state, t);
            // coordinates
            for(std::size_t i = 0; i < DIM; i++) {
                coordinates[i].resize(Type::size());
            }
            for(std::size_t i = 0; i < Type::size(); i++) {
                const TypeVector<DIM> position = sMemberStep->cX(cMemberState(state.data(), i));
                for(std::size_t j = 0; j < DIM; j++) {
                    coordinates[j][i] = position[j];
                }
            }
            // mesh
            double l = 0.0;
            TypeContainer<double> gridPoints(Type::size());
            gridPoints[0] = 0.0;
            for(std::size_t i = 1; i < Type::size(); i++){ 
                const auto x1 = sMemberStep->cX(cMemberState(state.data(), i));
                const auto x0 = sMemberStep->cX(cMemberState(state.data(), i - 1));
                l += (x1 - x0).norm();
                gridPoints[i] = l;
            }
            for(std::size_t i = 1; i < Type::size(); i++){ 
                gridPoints[i] /= l;
            }
            sMesh = std::make_shared<TypeMesh>(gridPoints, false);
            // update size
            const double newLength = length(state.data());
            const std::size_t newSize = std::ceil(newLength/dl);
            const double newDs = 1.0 / (newSize - 1);
            // manage data
            const int difference = newSize - Type::size();
            if(difference < 0) {
                for(int i = -1; i >= difference; i--) {
                    Type::removeMember(state);
                }
            } else if(difference > 0) {
                for(int i = 0; i < difference; i++) {
                    // create new member
                    std::shared_ptr<TypeMemberStep> sNewMemberStep = std::make_shared<TypeMemberStep>(*sMemberStep);
                    // add member
                    Type::addMember(state);
                }
            }
            // interpolation
            std::for_each(std::execution::par_unseq, memberIndexs.cbegin(), memberIndexs.cend(), [this, pState, t, &dState](const unsigned int& memberIndex){ 
                const TypeVector<1> s(memberIndex * newDs);
                TypeView<TypeVector<DIM>> x = sMemberStep->x(memberState(state.data(), i));
                for(std::size_t i = 0; i < DIM; i++) {
                    x[i] = p0l::lagrangeMeshPoint<TypeMesh, TypeContainer, double, TypeVector<1>, TypeRef, TypeMeshSub>(sMesh, coordinates[i], s, interpolationOrder + 1, false);
                }
            });
        }
    public:
        TypeVector<DIM> cX(const double* pState, const double& p_s) const {
            const TypeVector<1> s(p_s);
            // coordinates
            TypeContainer<TypeContainer<double>> tmpCoordinates(DIM, TypeContainer<double>(Type::size())); // TODO: do this better
            for(std::size_t i = 0; i < Type::size(); i++) {
                const TypeVector<DIM> position = sMemberStep->cX(cMemberState(pState, i));
                for(std::size_t j = 0; j < DIM; j++) {
                    tmpCoordinates[j][i] = position[j];
                }
            }
            // mesh
            double l = 0.0;
            TypeContainer<double> gridPoints(Type::size());
            gridPoints[0] = 0.0;
            for(std::size_t i = 1; i < Type::size(); i++){ 
                const auto x1 = sMemberStep->cX(cMemberState(pState, i));
                const auto x0 = sMemberStep->cX(cMemberState(pState, i - 1));
                l += (x1 - x0).norm();
                gridPoints[i] = l;
            }
            for(std::size_t i = 1; i < Type::size(); i++){ 
                gridPoints[i] /= l;
            }
            std::shared_ptr<TypeMesh> sTmpMesh = std::make_shared<TypeMesh>(gridPoints, false);
            // interpolate for each coordinates
            TypeVector<DIM> x;
            for(std::size_t j = 0; j < DIM; j++) {
                x[j] = p0l::lagrangeMeshPoint<TypeMesh, TypeContainer, double, TypeVector<1>, TypeRef, TypeMeshSub>(sTmpMesh, tmpCoordinates[j], s, interpolationOrder + 1, false);
            }
            return x;
        }
    public:
        double length(const double* pState) const {
            double l = 0.0;
            for(std::size_t i = 1; i < Type::size(); i++){ 
                const auto x1 = sMemberStep->cX(cMemberState(pState, i));
                const auto x0 = sMemberStep->cX(cMemberState(pState, i-1));
                l += (x1 - x0).norm();
            }
            return l;
        }
    public:
        // parameters
        double dl;
        unsigned int interpolationOrder;
        // shape description
        std::shared_ptr<TypeMesh> sMesh;
        std::vector<std::vector<double>> coordinates; // TODO: iterate avec all state, not just coordinates ?
    public:
        using Type::sMemberStep;
    public:
        using Type::cMemberState;
        using Type::memberState;
    public:
        using Type::sMemberSteps;
};

template<template<int> typename TypeVector, unsigned int DIM, template<typename...> class TypeView, template<typename...> class TypeRef, typename TypeMemberStep, typename TypeSolver>
class ChainDynamic : public ObjectDynamic<TypeSolver, StepChainDynamic<TypeVector, DIM, TypeView, TypeRef, TypeMemberStep>> {
    public:
        using TypeStep = StepChainDynamic<TypeVector, DIM, TypeView, TypeRef, TypeMemberStep>;
    public:
        template<typename ...Args>
        using TypeContainer = std::vector<Args...>;
    public:
        ChainDynamic(const std::shared_ptr<TypeMemberStep>& p_sMemberStep, const double& p_dl, const unsigned int& p_interpolationOrder) : ObjectDynamic<TypeSolver, TypeStep>(std::make_shared<TypeStep>(p_sMemberStep, p_dl, p_interpolationOrder)) {
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
