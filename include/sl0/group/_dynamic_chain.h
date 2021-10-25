#ifndef SL0_GROUP_CHAIN_H
#define SL0_GROUP_CHAIN_H
#pragma once

// std includes
#include <memory> // shared_ptr
#include <vector>
#include <numeric> // iota
#include <execution>
// module includes
#include "sl0/group/dynamic.h"

// TODO: add possibility to change execution policy, seq, par, par_unseq

namespace sl0 {

template<typename TypeStateDynamic, template<typename...> class TypeRef, template<typename...> class TypeView, typename TypeMemberStep>
class StepGroupChain : public StepGroupDynamic<TypeStateDynamic, TypeRef, TypeView, TypeMemberStep> {
    public:
        using TypeStepGroupDynamic = StepGroupDynamic<TypeStateDynamic, TypeRef, TypeView, TypeMemberStep>;
    public:
        StepGroupChain() : TypeStepGroupDynamic() {
        }

        void registerMember(std::shared_ptr<TypeMemberStep> sStep, const unsigned int& stateSize) override {
            registerMember(sStep, stateSize, TypeStepGroupDynamic::size());
        }

        void registerMember(std::shared_ptr<TypeMemberStep> sStep, const unsigned int& stateSize, const unsigned int& position) {
            TypeStepGroupDynamic::registerMember(sStep, stateSize);
            memberIndexsSorted.insert(memberIndexsSorted.begin() + position, TypeStepGroupDynamic::memberIndexs.back());
        }

        void unregisterMember(const unsigned int& memberIndex) override {
            TypeStepGroupDynamic::unregisterMember(memberIndex);
        }

        void updateMemberSize(const unsigned int& memberIndex, const int& difference) override {
            TypeStepGroupDynamic::updateMemberSize(memberIndex, difference);
        }
    public:
        // using
        //using TypeStepGroupDynamic::operator();
        //using TypeStepGroupDynamic::update();
        //using TypeStepGroupDynamic::stateSize;
        //using TypeStepGroupDynamic::registerMember;
        //using TypeStepGroupDynamic::unregisterMember;
        //using TypeStepGroupDynamic::updateMemberSize;
        //using TypeStepGroupDynamic::size;
        //using TypeStepGroupDynamic::cMemberState;
        //using TypeStepGroupDynamic::memberState;
    public:
        double dlMax;
        double lMax;
        std::vector<unsigned int> memberIndexsSorted;
        // using
        //using TypeStepGroupDynamic::memberIndexs;
        //using TypeStepGroupDynamic::sMemberSteps;
        //using TypeStepGroupDynamic::memberStateSizes;
        //using TypeStepGroupDynamic::memberStateIndexs;
};

template<typename TypeStateDynamic, template<typename...> class TypeRef, template<typename...> class TypeView, typename TypeMemberStep, typename TypeSolver>
class GroupChain : public ObjectDynamic<TypeStateDynamic, TypeRef, TypeSolver, StepGroupChain<TypeStateDynamic, TypeRef, TypeView, TypeMemberStep>> {
    public:
        using TypeStep = StepGroupChain<TypeStateDynamic, TypeRef, TypeView, TypeMemberStep>;
    public:
        GroupChain() : ObjectDynamic<TypeStateDynamic, TypeRef, TypeSolver, TypeStep>(std::make_shared<TypeStep>()) {
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
        void update(const double& dt) override {
            // refine if necessary
            for(std::size_t position = 1; position < sStep->memberIndexsSorted.size(); position++){ 
                const std::size_t index1 = sStep->memberIndexsSorted[position];
                const std::size_t index0 = sStep->memberIndexsSorted[position-1];
                const auto x1 = sStep->sMemberSteps[index1]->cX(sStep->cMemberState(state, index1));
                const auto x0 = sStep->sMemberSteps[index0]->cX(sStep->cMemberState(state, index0));
                const double dl = (x1 - x0).norm();
                if(dl > sStep->dlMax) {
                    // add point in the middle
                    state.conservativeResize(state.size() + sStep->sMemberSteps[index0]->stateSize());
                    sStep->registerMember(std::make_shared<TypeMemberStep>(*(sStep->sMemberSteps[index0])), sStep->sMemberSteps[index0]->stateSize(), position);
                    sStep->sMemberSteps[sStep->memberIndexs.back()]->x(sStep->memberState(state, sStep->memberIndexs.back())) = 0.5 * (x0 + x1);
                    // repass over new segments
                    position--;
                }
            }
            // classic object update
            ObjectDynamic<TypeStateDynamic, TypeRef, TypeSolver, TypeStep>::update(dt);
        }
    public:
        // Inherited properties
        using ObjectDynamic<TypeStateDynamic, TypeRef, TypeSolver, TypeStep>::sSolver;
        using ObjectDynamic<TypeStateDynamic, TypeRef, TypeSolver, TypeStep>::sStep;
        using ObjectDynamic<TypeStateDynamic, TypeRef, TypeSolver, TypeStep>::state;
        using ObjectDynamic<TypeStateDynamic, TypeRef, TypeSolver, TypeStep>::t;
};

}

#endif
