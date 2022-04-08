#ifndef SL0_SPHEROID_H
#define SL0_SPHEROID_H
#pragma once

#include "sl0/axis.h"

namespace sl0 {

template<template<int> typename TypeVector, unsigned int DIM, template<typename...> class TypeView, typename TypeFlow>
class StepSpheroid : public StepAxis<TypeVector, DIM, TypeView> {
	public:
		using TypeStepAxis = StepAxis<TypeVector, DIM, TypeView>;
		using TypeStepAxis::StateSize;
		using typename TypeStepAxis::TypeSpaceVector;
		using typename TypeStepAxis::TypeStateVectorDynamic;
	public:
		StepSpheroid(const std::shared_ptr<TypeFlow>& p_sFlow, const double& prop) : TypeStepAxis(), sFlow(p_sFlow) {
			setProportion(prop);
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
			auto skewVelocityGradients = 0.5 * (velocityGradients - velocityGradients.transpose());
			auto symVelocityGradients = 0.5 * (velocityGradients + velocityGradients.transpose());
			dAxis = skewVelocityGradients * sAxis + factor * (symVelocityGradients * sAxis - (sAxis.dot(symVelocityGradients * sAxis)) * sAxis);
			//TypeSpaceVector omega = sFlow->getVorticity(sX, t) + factor * sAxis.cross(sFlow->getStrain(sX, t) * sAxis);
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
		void setProportion(const double& prop) {
			factor = (std::pow(prop, 2) - std::pow(1.0, 2)) / (std::pow(prop, 2) + std::pow(1.0, 2));
		}
		double getFactor() const {
			return factor;
		}
	public:
		std::shared_ptr<TypeFlow> sFlow;
		double factor;
};

template<template<int> typename TypeVector, unsigned int DIM, template<typename...> class TypeView, typename TypeFlow, typename TypeSolver>
class Spheroid : public ObjectStatic<TypeSolver, StepSpheroid<TypeVector, DIM, TypeView, TypeFlow>> {
	public:
		using TypeStep = StepSpheroid<TypeVector, DIM, TypeView, TypeFlow>;
		using Type = ObjectStatic<TypeSolver, TypeStep>;
	public:
		Spheroid(const std::shared_ptr<TypeFlow>& sFlow, const double& prop) : Type::ObjectStatic(std::make_shared<TypeStep>(sFlow, prop)) {
		}
	public:
		using Type::sSolver;
		using Type::sStep;
		using Type::state;
		using Type::t;
};

}

#endif
