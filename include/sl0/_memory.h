#ifndef SL0_MEMORY_H
#define SL0_MEMORY_H
#pragma once

// std includes
#include <memory> // shared_ptr
// module includes
#include "sl0/object.h"

namespace sl0 {

template<typename TypeState, template<typename...> class TypeRef, template<typename...> class TypeView, typename TypeStepObject, typename TypeStateObject>
class StepMemory : public StepObject<TypeState, TypeRef> {
    public:
        StepMemory(const std::shared_ptr<TypeStepObject>& sStepObject, const double& Te);
        TypeState operator()(const TypeRef<const TypeState>& state, const double& t) const override;
    public:
        TypeView<const TypeStateObject> cMemoryState(const TypeRef<const TypeState>& state, const std::size_t& n) const;
        TypeView<TypeStateObject> memoryState(TypeRef<TypeState> state, const std::size_t& n) const;
        double cMemoryTime(const TypeRef<const TypeState>& state, const std::size_t& n) const;
        double& memoryTime(TypeRef<TypeState> state, const std::size_t& n) const;
    public:
        std::shared_ptr<TypeStepObject> sStepObject;
        double Te;
};

template<typename TypeState, template<typename...> class TypeRef, template<typename...> class TypeView, template<typename> class TypeSolver, typename TypeStepObject, typename TypeStateObject>
class ObjectMemory : public Object<TypeState, TypeRef, TypeSolver<StepMemory<TypeState, TypeRef, TypeView, TypeStepObject, TypeStateObject>>, StepMemory<TypeState, TypeRef, TypeView, TypeStepObject, TypeStateObject>> {
    public:
        using TypeStep = StepMemory<TypeState, TypeRef, TypeView, TypeStepObject, TypeStateObject>;
    public:
        ObjectMemory(const std::shared_ptr<TypeStepObject>& sStepObject, const double& Te);
    public:
        // Inherited
        using Object<TypeState, TypeRef, TypeSolver<TypeStep>, TypeStep>::sSolver;
        using Object<TypeState, TypeRef, TypeSolver<TypeStep>, TypeStep>::sStep;
        using Object<TypeState, TypeRef, TypeSolver<TypeStep>, TypeStep>::state;
        using Object<TypeState, TypeRef, TypeSolver<TypeStep>, TypeStep>::t;
};

// StepMemory class

template<typename TypeState, template<typename...> class TypeRef, template<typename...> class TypeView, typename TypeStepObject, typename TypeStateObject>
StepMemory<TypeState, TypeRef, TypeView, TypeStepObject, TypeStateObject>::StepMemory(const std::shared_ptr<TypeStepObject>& p_sStepObject, const double& p_Te) : StepObject<TypeState, TypeRef>(), sStepObject(p_sStepObject), Te(p_Te) {

}

// MemoryStates

template<typename TypeState, template<typename...> class TypeRef, template<typename...> class TypeView, typename TypeStepObject, typename TypeStateObject>
TypeView<const TypeStateObject> StepMemory<TypeState, TypeRef, TypeView, TypeStepObject, TypeStateObject>::cMemoryState(const TypeRef<const TypeState>& state, const std::size_t& n) const {
    return TypeView<const TypeStateObject>(&(*(state.begin() + (TypeStateObject::SizeAtCompileTime + 1) * n)));
}

template<typename TypeState, template<typename...> class TypeRef, template<typename...> class TypeView, typename TypeStepObject, typename TypeStateObject>
TypeView<TypeStateObject> StepMemory<TypeState, TypeRef, TypeView, TypeStepObject, TypeStateObject>::memoryState(TypeRef<TypeState> state, const std::size_t& n) const {
    return TypeView<TypeStateObject>(&(*(state.begin() + (TypeStateObject::SizeAtCompileTime + 1) * n)));
}

// Time

template<typename TypeState, template<typename...> class TypeRef, template<typename...> class TypeView, typename TypeStepObject, typename TypeStateObject>
double StepMemory<TypeState, TypeRef, TypeView, TypeStepObject, TypeStateObject>::cMemoryTime(const TypeRef<const TypeState>& state, const std::size_t& n) const {
    return *(state.begin() + (TypeStateObject::SizeAtCompileTime + 1) * n + TypeStateObject::SizeAtCompileTime);
}

template<typename TypeState, template<typename...> class TypeRef, template<typename...> class TypeView, typename TypeStepObject, typename TypeStateObject>
double& StepMemory<TypeState, TypeRef, TypeView, TypeStepObject, TypeStateObject>::memoryTime(TypeRef<TypeState> state, const std::size_t& n) const {
    return *(state.begin() + (TypeStateObject::SizeAtCompileTime + 1) * n + TypeStateObject::SizeAtCompileTime);
}

// operator()

template<typename TypeState, template<typename...> class TypeRef, template<typename...> class TypeView, typename TypeStepObject, typename TypeStateObject>
TypeState StepMemory<TypeState, TypeRef, TypeView, TypeStepObject, TypeStateObject>::operator()(const TypeRef<const TypeState>& state, const double& p_t) const {
    // Init dState
    TypeState dState;
    dState.fill(0.0);
    // Set dState
    for(std::size_t n = state.size()/(TypeStateObject::SizeAtCompileTime + 1) - 1; n > 0; n--) {
        const double dt = cMemoryTime(state, n - 1) - cMemoryTime(state, n);
        if(dt > Te) {
            memoryTime(dState, n) = 1.0;
            memoryState(dState, n) = (cMemoryState(state, n-1) - cMemoryState(state, n)) / dt;
        }
        memoryTime(dState, 0) = 1.0;
        memoryState(dState, 0) = (*sStepObject)(cMemoryState(state, 0), p_t);
    }
    // Return result
    return dState;
}

// ObjectMemory class

template<typename TypeState, template<typename...> class TypeRef, template<typename...> class TypeView, template<typename> class TypeSolver, typename TypeStepObject, typename TypeStateObject>
ObjectMemory<TypeState, TypeRef, TypeView, TypeSolver, TypeStepObject, TypeStateObject>::ObjectMemory(const std::shared_ptr<TypeStepObject>& sStepObject, const double& Te) : Object<TypeState, TypeRef, TypeSolver<TypeStep>, TypeStep>(std::make_shared<TypeStep>(sStepObject, Te)) {

}

}

#endif
