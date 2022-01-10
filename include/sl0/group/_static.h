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

// TODO: add possibility to change execution policy

namespace sl0 {

template<template<int> typename TypeVector, unsigned int DIM, template<typename...> class TypeRef, unsigned int _Size, template<typename...> class TypeView, typename TypeMemberStep>
class StepGroupStatic : public StepObject<TypeVector, DIM, TypeRef> {
    public:
        using typename StepObject<TypeVector, DIM, TypeRef>::TypeStateDynamic;
        using typename StepObject<TypeVector, DIM, TypeRef>::TypeSpaceVector;
        // Member
        static const unsigned int MemberStateSize = TypeMemberStep::StateSize;
        using TypeMemberState = typename TypeMemberStep::TypeStateStatic;
        // More
        static const unsigned int Size = _Size;
    public:
        StepGroupStatic(const std::shared_ptr<TypeMemberStep>& p_sMemberStep) : sMemberStep(p_sMemberStep), memberIndexs(Size) {
            std::iota(memberIndexs.begin(), memberIndexs.end(), 0);
        }

        TypeStateDynamic operator()(const TypeRef<const TypeStateDynamic>& state, const double& t) const override {
            // Init dState
            TypeStateDynamic dState(Size * MemberStateSize);
            // Set dState
            std::for_each(std::execution::par_unseq, memberIndexs.cbegin(), memberIndexs.cend(), [this, state, t, &dState](const unsigned int& memberIndex){ 
                memberState(dState, memberIndex) = (*sMemberStep)(cMemberState(state, memberIndex), t); 
            });
            // Return result
            return dState;
        }
        
        void update(TypeRef<TypeStateDynamic> state, const double& t) override {
            std::for_each(std::execution::par_unseq, memberIndexs.cbegin(), memberIndexs.cend(), [this, &state, t](const unsigned int& memberIndex){ 
                sMemberStep->update(TypeRef<TypeMemberState>(memberState(state, memberIndex)), t);
            });
        }

        unsigned int stateSize() const override {
            return Size * MemberStateSize;
        }
    public:
        std::size_t size() const {
            return Size;
        }
        TypeView<const TypeMemberState> cMemberState(const TypeRef<const TypeStateDynamic>& state, const std::size_t& n) const {
            return TypeView<const TypeMemberState>(&(*(state.begin() + MemberStateSize * n)));
        }
        TypeView<TypeMemberState> memberState(TypeRef<TypeStateDynamic> state, const std::size_t& n) const {
            return TypeView<TypeMemberState>(&(*(state.begin() + MemberStateSize * n)));
        }
    public:
        std::vector<TypeSpaceVector> positions(const TypeRef<const TypeStateDynamic>& state) const {
            std::vector<TypeSpaceVector> result;
            result.reserve(Size * sMemberStep->positions(cMemberState(state, 0)).size());
            for(unsigned int i = 0; i < Size; i++) {
                const std::vector<TypeSpaceVector> memberPositions = sMemberStep->positions(cMemberState(state, i));
                result.insert(result.end(), memberPositions.begin(), memberPositions.end());
            }
            return result;
        }
    public:
        std::shared_ptr<TypeMemberStep> sMemberStep;
        std::vector<unsigned int> memberIndexs;
};

template<template<int> typename TypeVector, unsigned int DIM, template<typename...> class TypeRef, unsigned int _Size, template<typename...> class TypeView, typename TypeMemberStep, typename TypeSolver>
class GroupStatic : public ObjectDynamic<TypeSolver, StepGroupStatic<TypeVector, DIM, TypeRef, _Size, TypeView, TypeMemberStep>> {
    public:
        using TypeStep = StepGroupStatic<TypeVector, DIM, TypeRef, _Size, TypeView, TypeMemberStep>;
    public:
        GroupStatic(const std::shared_ptr<TypeMemberStep>& p_sMemberStep) : ObjectDynamic<TypeSolver, TypeStep>(std::make_shared<TypeStep>(p_sMemberStep)) {
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
