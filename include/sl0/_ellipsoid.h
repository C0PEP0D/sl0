#ifndef SL0_ELLIPSOID_H
#define SL0_ELLIPSOID_H
#pragma once

#include <memory>
#include <cmath>
#include <iostream>

#include "sl0/object.h"

namespace sl0 {

template<typename TypeState, template<typename...> class TypeRef, template<typename...> class TypeView, typename TypeMatrix>
class StepBasis : public StepObject<TypeState, TypeRef> {
    public:
        using StepObject<TypeState, TypeRef>::StepObject;
        void update(TypeRef<TypeState> state, const double& t) override;
        // sub update methods
        void orthoNormalizeBasis(TypeRef<TypeState> state);
    public:
        virtual TypeView<const TypeMatrix> cBasis(const TypeRef<const TypeState>& p_state) const = 0;
        virtual TypeView<TypeMatrix> basis(TypeRef<TypeState> p_state) const = 0;
};

template<typename TypeState, template<typename...> class TypeRef, template<typename...> class TypeView, typename TypeVector, typename TypeMatrix, typename TypeFlow>
class StepEllipsoid : public StepBasis<TypeState, TypeRef, TypeView, TypeMatrix> {
    public:
        StepEllipsoid(const std::shared_ptr<TypeFlow>& sFlow, const TypeVector& props);
        TypeState operator()(const TypeRef<const TypeState>& state, const double& t) const override;
    public:
        TypeView<const TypeVector> cX(const TypeRef<const TypeState>& state) const;
        TypeView<TypeVector> x(TypeRef<TypeState> state) const;
        TypeView<const TypeMatrix> cBasis(const TypeRef<const TypeState>& p_state) const override;
        TypeView<TypeMatrix> basis(TypeRef<TypeState> p_state) const override;
    public:
        void setProportions(const TypeVector& props);
        TypeView<const TypeVector> getProportions() const;
    public:
        std::shared_ptr<TypeFlow> sFlow;
        TypeVector factors;
};

template<typename TypeState, template<typename...> class TypeRef, template<typename...> class TypeView, template<typename> class TypeSolver, typename TypeVector, typename TypeMatrix, typename TypeFlow>
class Ellipsoid : public Object<TypeState, TypeRef, TypeSolver<StepEllipsoid<TypeState, TypeRef, TypeView, TypeVector, TypeMatrix, TypeFlow>>, StepEllipsoid<TypeState, TypeRef, TypeView, TypeVector, TypeMatrix, TypeFlow>> {
    public:
        using TypeStep = StepEllipsoid<TypeState, TypeRef, TypeView, TypeVector, TypeMatrix, TypeFlow>;
    public:
        Ellipsoid(const std::shared_ptr<TypeFlow>& sFlow, const TypeVector& props);
    public:
        using Object<TypeState, TypeRef, TypeSolver<TypeStep>, TypeStep>::sSolver;
        using Object<TypeState, TypeRef, TypeSolver<TypeStep>, TypeStep>::sStep;
        using Object<TypeState, TypeRef, TypeSolver<TypeStep>, TypeStep>::state;
        using Object<TypeState, TypeRef, TypeSolver<TypeStep>, TypeStep>::t;
};

// StepBasis

// Gramm-Schmit
template<typename TypeState, template<typename...> class TypeRef, template<typename...> class TypeView, typename TypeMatrix>
void StepBasis<TypeState, TypeRef, TypeView, TypeMatrix>::orthoNormalizeBasis(TypeRef<TypeState> p_state) {
    auto b = basis(p_state);
    // First vector just gets normalized
    b.col(0).normalize();
    b.col(1) -= b.col(0) * b.col(1).dot(b.col(0));
    b.col(1).normalize();
    b.col(2) = b.col(0).cross(b.col(1));
}

template<typename TypeState, template<typename...> class TypeRef, template<typename...> class TypeView, typename TypeMatrix>
void StepBasis<TypeState, TypeRef, TypeView, TypeMatrix>::update(TypeRef<TypeState> p_state, const double& t) {
    StepObject<TypeState, TypeRef>::update(p_state, t);
    // Just orthonormalize basis
    orthoNormalizeBasis(p_state);
}


// StepEllispoid

template<typename TypeState, template<typename...> class TypeRef, template<typename...> class TypeView, typename TypeVector, typename TypeMatrix, typename TypeFlow>
StepEllipsoid<TypeState, TypeRef, TypeView, TypeVector, TypeMatrix, TypeFlow>::StepEllipsoid(const std::shared_ptr<TypeFlow>& p_sFlow, const TypeVector& props) : StepBasis<TypeState, TypeRef, TypeView, TypeMatrix>(), sFlow(p_sFlow) {
    setProportions(props);
}

// proportions

template<typename TypeState, template<typename...> class TypeRef, template<typename...> class TypeView, typename TypeVector, typename TypeMatrix, typename TypeFlow>
void StepEllipsoid<TypeState, TypeRef, TypeView, TypeVector, TypeMatrix, TypeFlow>::setProportions(const TypeVector& props) {
    factors[0] = (std::pow(props[1], 2) - std::pow(props[2], 2)) / (std::pow(props[1], 2) + std::pow(props[2], 2));
    factors[1] = (std::pow(props[2], 2) - std::pow(props[0], 2)) / (std::pow(props[2], 2) + std::pow(props[0], 2));
    factors[2] = (std::pow(props[0], 2) - std::pow(props[1], 2)) / (std::pow(props[0], 2) + std::pow(props[1], 2));
}

template<typename TypeState, template<typename...> class TypeRef, template<typename...> class TypeView, typename TypeVector, typename TypeMatrix, typename TypeFlow>
TypeView<const TypeVector> StepEllipsoid<TypeState, TypeRef, TypeView, TypeVector, TypeMatrix, TypeFlow>::getProportions() const {
    return TypeView<const TypeVector>(factors.data());
}

// x

template<typename TypeState, template<typename...> class TypeRef, template<typename...> class TypeView, typename TypeVector, typename TypeMatrix, typename TypeFlow>
TypeView<const TypeVector> StepEllipsoid<TypeState, TypeRef, TypeView, TypeVector, TypeMatrix, TypeFlow>::cX(const TypeRef<const TypeState>& p_state) const {
    return TypeView<const TypeVector>(p_state.data());
}

template<typename TypeState, template<typename...> class TypeRef, template<typename...> class TypeView, typename TypeVector, typename TypeMatrix, typename TypeFlow>
TypeView<TypeVector> StepEllipsoid<TypeState, TypeRef, TypeView, TypeVector, TypeMatrix, TypeFlow>::x(TypeRef<TypeState> p_state) const {
    return TypeView<TypeVector>(p_state.data());
}

// basis

template<typename TypeState, template<typename...> class TypeRef, template<typename...> class TypeView, typename TypeVector, typename TypeMatrix, typename TypeFlow>
TypeView<const TypeMatrix> StepEllipsoid<TypeState, TypeRef, TypeView, TypeVector, TypeMatrix, TypeFlow>::cBasis(const TypeRef<const TypeState>& p_state) const {
    return TypeView<const TypeMatrix>(p_state.data() + 3);
}

template<typename TypeState, template<typename...> class TypeRef, template<typename...> class TypeView, typename TypeVector, typename TypeMatrix, typename TypeFlow>
TypeView<TypeMatrix> StepEllipsoid<TypeState, TypeRef, TypeView, TypeVector, TypeMatrix, TypeFlow>::basis(TypeRef<TypeState> p_state) const {
    return TypeView<TypeMatrix>(p_state.data() + 3);
}

// operator()

template<typename TypeState, template<typename...> class TypeRef, template<typename...> class TypeView, typename TypeVector, typename TypeMatrix, typename TypeFlow>
TypeState StepEllipsoid<TypeState, TypeRef, TypeView, TypeVector, TypeMatrix, TypeFlow>::operator()(const TypeRef<const TypeState>& p_state, const double& p_t) const {
    TypeState dState; dState.fill(0.0);
    // Get data from state
    const TypeView<const TypeVector> sX = cX(p_state);
    const TypeView<const TypeMatrix> sBasis = cBasis(p_state);
    // Prepare to set data
    TypeView<TypeVector> dX = x(dState);
    TypeView<TypeMatrix> dBasis = basis(dState);
    // Compute linear velocity
    dX = sFlow->getVelocity(sX, p_t);
    // Compute rotaton velocity
    // Axysymetric solution : TypeVector omega = sFlow->getVorticity(x, p_t) + factors(2) * (basis.col(0).cross(sFlow->getStrain(x, p_t) * basis.col(0)));
    TypeVector omega = sFlow->getVorticity(sX, p_t);
    for(std::size_t i = 0; i < 3; i++) {
        omega += factors(i) * (sBasis.col((i + 1) % 3).dot(sFlow->getStrain(sX, p_t) * sBasis.col((i + 2) % 3))) * sBasis.col(i);
    }
    //// Express the derivatives of the basis vectors
    for(size_t i = 0; i < 3; i++){
        dBasis.col(i) = omega.cross(sBasis.col(i));
    }
    // Return 
    return dState;
}

// Ellipsoid

template<typename TypeState, template<typename...> class TypeRef, template<typename...> class TypeView, template<typename> class TypeSolver, typename TypeVector, typename TypeMatrix, typename TypeFlow>
Ellipsoid<TypeState, TypeRef, TypeView, TypeSolver, TypeVector, TypeMatrix, TypeFlow>::Ellipsoid(const std::shared_ptr<TypeFlow>& p_sFlow, const TypeVector& props) : Object<TypeState, TypeRef, TypeSolver<TypeStep>, TypeStep>::Object(std::make_shared<TypeStep>(p_sFlow, props)) {

}

}

#endif
