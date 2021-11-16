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

template<template<int> typename TypeVector, unsigned int DIM, template<typename...> class TypeRef, template<typename...> class TypeView, typename TypeMemberStep>
class StepChainDynamic : public StepGroupDynamic<TypeVector, DIM, TypeRef, TypeView, TypeMemberStep> {
    public:
        using TypeStepGroupDynamic = StepGroupDynamic<TypeVector, DIM, TypeRef, TypeView, TypeMemberStep>;
        using typename TypeStepGroupDynamic::TypeStateDynamic;
    public:
        template<typename ...Args>
        using TypeContainer = std::vector<Args...>;
        using TypeMesh = m0sh::NonUniform<TypeVector<1>, TypeRef, TypeContainer>;
        using TypeMeshSub = m0sh::StructuredSub<TypeVector<1>, TypeRef, TypeContainer>;
    public:
        StepChainDynamic(const double& p_dl, const unsigned int& p_interpolationOrder) : TypeStepGroupDynamic(), dl(p_dl), interpolationOrder(p_interpolationOrder), coordinates(DIM) {
        }

        void update(TypeRef<TypeStateDynamic> state, const double& t) override {
            TypeStepGroupDynamic::update(state, t);
            // coordinates
            for(std::size_t i = 0; i < DIM; i++) {
                coordinates[i].resize(TypeStepGroupDynamic::size());
            }
            for(std::size_t i = 0; i < TypeStepGroupDynamic::size(); i++) {
                const TypeVector<DIM> position = sMemberSteps[i]->cX(cMemberState(state, i));
                for(std::size_t j = 0; j < DIM; j++) {
                    coordinates[j][i] = position[j];
                }
            }
            // mesh
            double l = 0.0;
            TypeContainer<double> gridPoints(TypeStepGroupDynamic::size());
            gridPoints[0] = 0.0;
            for(std::size_t i = 1; i < TypeStepGroupDynamic::size(); i++){ 
                const auto x1 = sMemberSteps[i]->cX(cMemberState(state, i));
                const auto x0 = sMemberSteps[i]->cX(cMemberState(state, i - 1));
                l += (x1 - x0).norm();
                gridPoints[i] = l;
            }
            for(std::size_t i = 1; i < TypeStepGroupDynamic::size(); i++){ 
                gridPoints[i] /= l;
            }
            sMesh = std::make_shared<TypeMesh>(gridPoints, false);
            // update size
            newLength = length(state);
            newSize = std::ceil(newLength/dl);
            newDs = 1.0 / (newSize - 1);
            // interpolation
            for(std::size_t i = 0; i < TypeStepGroupDynamic::size(); i++) {
                const TypeVector<1> s(i * newDs);
                TypeView<TypeVector<DIM>> x = sMemberSteps[i]->x(memberState(state, i));
                for(std::size_t j = 0; j < DIM; j++) {
                    x[j] = p0l::lagrangeMeshPoint<TypeMesh, TypeContainer, double, TypeVector<1>, TypeRef, TypeMeshSub>(sMesh, coordinates[j], s, interpolationOrder + 1, false);
                }
            }
        }
    public:
        TypeVector<DIM> cX(const TypeRef<const TypeStateDynamic> state, const double& p_s) const {
            const TypeVector<1> s(p_s);
            // coordinates
            TypeContainer<TypeContainer<double>> tmpCoordinates(DIM, TypeContainer<double>(TypeStepGroupDynamic::size())); // TODO: do this better
            for(std::size_t i = 0; i < TypeStepGroupDynamic::size(); i++) {
                const TypeVector<DIM> position = sMemberSteps[i]->cX(cMemberState(state, i));
                for(std::size_t j = 0; j < DIM; j++) {
                    tmpCoordinates[j][i] = position[j];
                }
            }
            // mesh
            double l = 0.0;
            TypeContainer<double> gridPoints(TypeStepGroupDynamic::size());
            gridPoints[0] = 0.0;
            for(std::size_t i = 1; i < TypeStepGroupDynamic::size(); i++){ 
                const auto x1 = sMemberSteps[i]->cX(cMemberState(state, i));
                const auto x0 = sMemberSteps[i]->cX(cMemberState(state, i - 1));
                l += (x1 - x0).norm();
                gridPoints[i] = l;
            }
            for(std::size_t i = 1; i < TypeStepGroupDynamic::size(); i++){ 
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
        double length(const TypeRef<const TypeStateDynamic>& state) const {
            double l = 0.0;
            for(std::size_t i = 1; i < TypeStepGroupDynamic::size(); i++){ 
                const auto x1 = sMemberSteps[i]->cX(cMemberState(state, i));
                const auto x0 = sMemberSteps[i]->cX(cMemberState(state, i-1));
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
        double newLength;
        std::size_t newSize;
        double newDs;
    public:
        using TypeStepGroupDynamic::cMemberState;
        using TypeStepGroupDynamic::memberState;
    public:
        using TypeStepGroupDynamic::sMemberSteps;
};

template<template<int> typename TypeVector, unsigned int DIM, template<typename...> class TypeRef, template<typename...> class TypeView, typename TypeMemberStep, typename TypeSolver>
class ChainDynamic : public ObjectDynamic<TypeSolver, StepChainDynamic<TypeVector, DIM, TypeRef, TypeView, TypeMemberStep>> {
    public:
        using TypeStep = StepChainDynamic<TypeVector, DIM, TypeRef, TypeView, TypeMemberStep>;
    public:
        template<typename ...Args>
        using TypeContainer = std::vector<Args...>;
    public:
        ChainDynamic(const std::shared_ptr<TypeMemberStep>& p_sMemberStep, const double& p_dl, const unsigned int& p_interpolationOrder) : ObjectDynamic<TypeSolver, TypeStep>(std::make_shared<TypeStep>(p_dl, p_interpolationOrder)), sMemberStep(p_sMemberStep) {
        }
    public:
        // member management
        void addMember(std::shared_ptr<TypeMemberStep> p_sMemberStep, const unsigned int& stateSize) {
            state.conservativeResize(state.size() + stateSize);
            sStep->registerMember(p_sMemberStep, stateSize);
        }
        void removeMember(const unsigned int& memberIndex) {
            std::copy(state.begin() + memberIndex + sStep->memberStateSizes[memberIndex], state.end(), state.begin() + memberIndex);
            state.conservativeResize(state.size() - sStep->memberStateSizes[memberIndex]);
            sStep->unregisterMember(memberIndex);
        }
        void updateMemberSize(const unsigned int& memberIndex, const unsigned int& stateSize) {
            const int difference = stateSize - sStep->memberStateSizes[memberIndex];
            if(difference > 0) {
                state.conservativeResize(state.size() + difference);
                std::copy(state.begin() + memberIndex + sStep->memberStateSizes[memberIndex], state.end() - difference, state.begin() + stateSize);
            } else {
                std::copy(state.begin() + memberIndex + sStep->memberStateSizes[memberIndex], state.end(), state.begin() + stateSize);
                state.conservativeResize(state.size() + difference);
            }
            sStep->updateMemberSize(memberIndex, difference);
        }
    public:
        void update(const double& dt) override {
            ObjectDynamic<TypeSolver, TypeStep>::update(dt);
            // manage data
            const int difference = sStep->newSize - sStep->size();
            if(difference < 0) {
                for(int i = -1; i >= difference; i--) {
                    removeMember(sStep->size() + i);
                }
            } else if(difference > 0) {
                for(int i = 0; i < difference; i++) {
                    // create new member
                    std::shared_ptr<TypeMemberStep> sNewMemberStep = std::make_shared<TypeMemberStep>(*sMemberStep);
                    // add member
                    addMember(sNewMemberStep, sNewMemberStep->stateSize());
                    // interpolate position of the member
                    const TypeVector<1> s((sStep->size() - 1) * sStep->newDs);
                    TypeView<TypeVector<DIM>> x = sNewMemberStep->x(sStep->memberState(state, sStep->size() - 1));
                    for(std::size_t j = 0; j < DIM; j++) {
                        x[j] = p0l::lagrangeMeshPoint<typename TypeStep::TypeMesh, TypeContainer, double, TypeVector<1>, TypeRef, typename TypeStep::TypeMeshSub>(sStep->sMesh, sStep->coordinates[j], s, sStep->interpolationOrder + 1, false);
                    }
                }
            }
        }
    public:
        std::shared_ptr<TypeMemberStep> sMemberStep;
    public:
        // Inherited
        using ObjectDynamic<TypeSolver, TypeStep>::sSolver;
        using ObjectDynamic<TypeSolver, TypeStep>::sStep;
        using ObjectDynamic<TypeSolver, TypeStep>::state;
        using ObjectDynamic<TypeSolver, TypeStep>::t;
};

}

#endif
