#ifndef SL0_ZERMELO_H
#define SL0_ZERMELO_H
#pragma once

#include "sl0/axis.h"

namespace sl0 {

template<template<int> typename TypeVector, unsigned int DIM, template<typename...> class TypeView, typename TypeFlow>
class StepZermelo : public StepAxis<TypeVector, DIM, TypeView> {
	public:
		using TypeStepAxis = StepAxis<TypeVector, DIM, TypeView>;
		using TypeStepAxis::StateSize;
		using typename TypeStepAxis::TypeSpaceVector;
		using typename TypeStepAxis::TypeStateVectorDynamic;
	public:
		StepZermelo(const std::shared_ptr<TypeFlow>& p_sFlow, const double& p_factor) : TypeStepAxis(), sFlow(p_sFlow), factor(p_factor) {
		}
		TypeStateVectorDynamic operator()(const double* pState, const double& t) const override {
			TypeStateVectorDynamic dState(TypeStepAxis::stateSize());
			// Get data from state
			const TypeView<const TypeSpaceVector> sX = cX(pState);
			const TypeView<const TypeSpaceVector> sAxis = cAxis(pState);
			// Prepare to set data
			TypeView<TypeSpaceVector> dX = x(dState.data());
			TypeView<TypeSpaceVector> dAxis = axis(dState.data());
			// Compute linear velocity
			dX = sFlow->getVelocity(sX, t);
			// Compute axis velocity (rotaton velocity)
			auto velocityGradients = sFlow->getVelocityGradients(sX, t);
			dAxis = -factor * velocityGradients.transpose() * sAxis;
			// Return
			return dState;
		}
	public:
		TypeView<const TypeSpaceVector> cX(const double* pState) const {
			return TypeView<const TypeSpaceVector>(pState);
		}
		TypeView<TypeSpaceVector> x(double* pState) const {
			return TypeView<TypeSpaceVector>(pState);
		}
		TypeView<const TypeSpaceVector> cAxis(const double* pState) const override {
			return TypeView<const TypeSpaceVector>(pState + DIM);
		}
		TypeView<TypeSpaceVector> axis(double* pState) const override {
			return TypeView<TypeSpaceVector>(pState + DIM);
		}
	public:
		std::vector<TypeSpaceVector> positions(const double* pState) const override {
			return { cX(pState) };
		}
	public:
		std::shared_ptr<TypeFlow> sFlow;
		double factor;
};

template<template<int> typename TypeVector, unsigned int DIM, template<typename...> class TypeView, typename TypeFlow, typename TypeSolver>
class Zermelo : public ObjectStatic<TypeSolver, StepZermelo<TypeVector, DIM, TypeView, TypeFlow>> {
	public:
		using TypeStep = StepZermelo<TypeVector, DIM, TypeView, TypeFlow>;
		using Type = ObjectStatic<TypeSolver, TypeStep>;
	public:
		Zermelo(const std::shared_ptr<TypeFlow>& sFlow, const double& factor) : Type::ObjectStatic(std::make_shared<TypeStep>(sFlow, factor)) {
		}
	public:
		using Type::sSolver;
		using Type::sStep;
		using Type::state;
		using Type::t;
};

}

#endif
