#ifndef SL0_GROUP_STATIC_H
#define SL0_GROUP_STATIC_H
#pragma once

// std includes
#include <memory> // shared_ptr
#include <vector>
#include <numeric> // iota
#include <execution>
// module includes
#include "sl0/object.h"

// TODO: add possibility to change execution policy, seq, par, par_unseq

namespace sl0 {

template<template<int> typename TypeVector, unsigned int DIM, template<typename...> class TypeView, typename _TypeMemberStep, unsigned int Size>
class StepGroupStatic : public StepObjectStatic<TypeVector, DIM, Size * _TypeMemberStep::StateSize> {
    public:
        using Type = StepObjectStatic<TypeVector, DIM, Size * _TypeMemberStep::StateSize>;
        using typename Type::StateSize;
        using typename Type::TypeSpaceVector;
        using typename Type::TypeStateVectorDynamic;
        using TypeMemberStep = _TypeMemberStep;
    public:
        StepGroupStatic(std::shared_ptr<TypeMemberStep> p_sMemberStep) : sMemberStep(p_sMemberStep), sMemberSteps(Size, std::make_shared<TypeMemberStep>(*p_sMemberStep)), memberIndexs(Size) {
            std::iota(memberIndexs.begin(), memberIndexs.end(), 0);
        }
    public:
        // Main
        TypeStateVectorDynamic operator()(const double* pState, const double& t) const override {
            // Init dState
            TypeStateVectorDynamic dState;
            dState.resize(stateSize());
            // Set dState
            std::for_each(std::execution::par_unseq, memberIndexs.cbegin(), memberIndexs.cend(), [this, pState, t, &dState](const unsigned int& memberIndex){ 
                TypeView<TypeStateVectorDynamic> dMemberState(memberState(dState.data(), memberIndex), TypeMemberStep::StateSize); 
                dMemberState = (*sMemberSteps[memberIndex])(cMemberState(pState, memberIndex), t); 
            });
            // Return result
            return dState;
        }

        void update(double* pState, const double& t) override {
            std::for_each(std::execution::par_unseq, memberIndexs.cbegin(), memberIndexs.cend(), [this, pState, t](const unsigned int& memberIndex){
                sMemberSteps[memberIndex]->update(memberState(pState, memberIndex), t);
            });
        }

        unsigned int stateSize() const override {
            return Size * sMemberStep->stateSize();
        }
    public:
        // getters
        std::size_t size() const {
            return Size;
        }

        const double* cMemberState(const double* pState, const std::size_t& memberIndex) const {
            return pState + memberIndex * sMemberStep->stateSize();
        }

        double* memberState(double* pState, const std::size_t& memberIndex) const {
            return pState + memberIndex * sMemberStep->stateSize();
        }
    public:
        std::vector<TypeSpaceVector> positions(const double* pState) const override {
            std::vector<TypeSpaceVector> result;
            //result.reserve(size() * sMemberSteps[0]->positions(cMemberState(state, 0)).size());
            for(unsigned int i = 0; i < size(); i++) {
                const std::vector<TypeSpaceVector> memberPositions = sMemberSteps[i]->positions(cMemberState(pState, i));
                result.insert(result.end(), memberPositions.begin(), memberPositions.end());
            }
            return result;
        }
    public:
        std::vector<unsigned int> memberIndexs;
        std::vector<std::shared_ptr<TypeMemberStep>> sMemberSteps;
    public:
        std::shared_ptr<TypeMemberStep> sMemberStep;
};

template<template<int> typename TypeVector, unsigned int DIM, template<typename...> class TypeView, typename TypeMemberStep, unsigned int Size, typename TypeSolver>
class GroupStatic : public ObjectStatic<TypeSolver, StepGroupStatic<TypeVector, DIM, TypeView, TypeMemberStep, Size>> {
    public:
        using TypeStep = StepGroupStatic<TypeVector, DIM, TypeView, TypeMemberStep, Size>;
        using Type = ObjectStatic<TypeSolver, TypeStep>;
    public:
        GroupStatic(std::shared_ptr<TypeMemberStep> sMemberStep) : Type(std::make_shared<TypeStep>(sMemberStep)) {
        }
    public:
        // Inherited properties
        using Type::sSolver;
        using Type::sStep;
        using Type::state;
        using Type::t;
};

}

#endif
