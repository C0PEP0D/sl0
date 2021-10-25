#ifndef SL0_GROUP_DYNAMIC_H
#define SL0_GROUP_DYNAMIC_H
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

template<template<int> typename TypeVector, unsigned int DIM, template<typename...> class TypeRef, template<typename...> class TypeView, typename TypeMemberStep>
class StepGroupDynamic : public StepObject<TypeVector, DIM, TypeRef> {
    public:
        using typename StepObject<TypeVector, DIM, TypeRef>::TypeStateDynamic;
        using typename StepObject<TypeVector, DIM, TypeRef>::TypeSpaceVector;
    public:
        StepGroupDynamic() {
        }
    public:
        // Main
        TypeStateDynamic operator()(const TypeRef<const TypeStateDynamic>& state, const double& t) const override {
            // Init dState
            TypeStateDynamic dState;
            dState.resize(stateSize());
            // Set dState
            std::vector<unsigned int> memberIndexs(sMemberSteps.size());
            std::iota(memberIndexs.begin(), memberIndexs.end(), 0);
            std::for_each(std::execution::par_unseq, memberIndexs.cbegin(), memberIndexs.cend(), [this, state, t, &dState](const unsigned int& memberIndex){ 
                memberState(dState, memberIndex) = (*sMemberSteps[memberIndex])(cMemberState(state, memberIndex), t); 
            });
            // Return result
            return dState;
        }

        void update(TypeRef<TypeStateDynamic> state, const double& t) override {
            std::vector<unsigned int> memberIndexs(sMemberSteps.size());
            std::iota(memberIndexs.begin(), memberIndexs.end(), 0);
            std::for_each(std::execution::par_unseq, memberIndexs.cbegin(), memberIndexs.cend(), [this, &state, t](const unsigned int& memberIndex){
                sMemberSteps[memberIndex]->update(TypeRef<TypeStateDynamic>(memberState(state, memberIndex)), t);
            });
        }

        unsigned int stateSize() const override {
            if (not memberStateIndexs.empty()) {
                return memberStateIndexs.back() + memberStateSizes.back();
            } else {
                return 0;
            }
        }
    public:
        // registering
        virtual void registerMember(std::shared_ptr<TypeMemberStep> sStep, const unsigned int& stateSize) {
            if (not memberStateIndexs.empty()) {
                memberIndexs.push_back(memberIndexs.back() + 1);
                memberStateIndexs.push_back(memberStateIndexs.back() + memberStateSizes.back());
            } else {
                memberIndexs.push_back(0);
                memberStateIndexs.push_back(0);
            }
            sMemberSteps.push_back(sStep);
            memberStateSizes.push_back(stateSize);
        }

        virtual void unregisterMember(const unsigned int& memberIndex) {
            sMemberSteps.erase(sMemberSteps.begin() + memberIndex);
            // TODO STL
            for(unsigned int i = memberIndex + 1; i < sMemberSteps.size(); i++) {
                memberStateIndexs[i] -= memberStateSizes[memberIndex];
            }
            memberStateIndexs.erase(memberStateIndexs.begin() + memberIndex);
            memberStateSizes.erase(memberStateSizes.begin() + memberIndex);
        }

        virtual void updateMemberSize(const unsigned int& memberIndex, const int& difference) {
            memberStateSizes[memberIndex] += difference;
            // TODO STL
            for(unsigned int i = memberIndex + 1; i < sMemberSteps.size(); i++) {
                memberStateIndexs[i] += difference;
            }
        }
    public:
        // getters
        std::size_t size() const {
            return sMemberSteps.size();
        }

        TypeView<const TypeStateDynamic> cMemberState(const TypeRef<const TypeStateDynamic>& state, const std::size_t& memberIndex) const {
            return TypeView<const TypeStateDynamic>(state.data() + memberStateIndexs[memberIndex], memberStateSizes[memberIndex]);
        }

        TypeView<TypeStateDynamic> memberState(TypeRef<TypeStateDynamic> state, const std::size_t& memberIndex) const {
            return TypeView<TypeStateDynamic>(state.data() + memberStateIndexs[memberIndex], memberStateSizes[memberIndex]);
        }
    public:
        std::vector<TypeSpaceVector> positions(const TypeRef<const TypeStateDynamic>& state) const {
            std::vector<TypeSpaceVector> result;
            result.reserve(size() * sMemberSteps[0]->positions(cMemberState(state, 0)).size());
            for(unsigned int i = 0; i < size(); i++) {
                const std::vector<TypeSpaceVector> memberPositions = sMemberSteps[i]->positions(cMemberState(state, i));
                result.insert(result.end(), memberPositions.begin(), memberPositions.end());
            }
            return result;
        }
    public:
        std::vector<unsigned int> memberIndexs;
        std::vector<std::shared_ptr<TypeMemberStep>> sMemberSteps;
        std::vector<unsigned int> memberStateSizes;
        std::vector<unsigned int> memberStateIndexs;
};

template<template<int> typename TypeVector, unsigned int DIM, template<typename...> class TypeRef, template<typename...> class TypeView, typename TypeMemberStep, typename TypeSolver>
class GroupDynamic : public ObjectDynamic<TypeSolver, StepGroupDynamic<TypeVector, DIM, TypeRef, TypeView, TypeMemberStep>> {
    public:
        using TypeStep = StepGroupDynamic<TypeVector, DIM, TypeRef, TypeView, TypeMemberStep>;
    public:
        GroupDynamic() : ObjectDynamic<TypeSolver, TypeStep>(std::make_shared<TypeStep>()) {
        }
    public:
        // member management
        void addMember(std::shared_ptr<TypeMemberStep> sMemberStep, const unsigned int& stateSize) {
            state.conservativeResize(state.size() + stateSize);
            sStep->registerMember(sMemberStep, stateSize);
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
        // Inherited properties
        using ObjectDynamic<TypeSolver, TypeStep>::sSolver;
        using ObjectDynamic<TypeSolver, TypeStep>::sStep;
        using ObjectDynamic<TypeSolver, TypeStep>::state;
        using ObjectDynamic<TypeSolver, TypeStep>::t;
};

}

#endif
