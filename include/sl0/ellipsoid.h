#ifndef SL0_ELLIPSOID_H
#define SL0_ELLIPSOID_H
#pragma once

#include "sl0/basis.h"

namespace sl0 {

template<template<int> typename TypeVector, template<int, int> typename TypeMatrix, unsigned int DIM, template<typename...> class TypeView, typename TypeFlow>
class StepEllipsoid : public StepBasis<TypeVector, TypeMatrix, DIM, TypeView> {
	public:
		using TypeStepBasis = StepBasis<TypeVector, TypeMatrix, DIM, TypeView>;
		using TypeStepBasis::StateSize;
		using typename TypeStepBasis::TypeSpaceVector;
		using typename TypeStepBasis::TypeStateVectorDynamic;
	public:
		StepEllipsoid(const std::shared_ptr<TypeFlow>& p_sFlow, const std::vector<double>& dimensions) : TypeStepBasis(), sFlow(p_sFlow) {
			setDimensions(dimensions);
		}
		TypeStateVectorDynamic operator()(const double* pState, const double& t) const override {
			TypeStateVectorDynamic dState(TypeStepBasis::stateSize());
			// Get data from state
			const TypeView<const TypeSpaceVector> sX = cX(pState);
			const TypeView<const TypeSpaceVector> sAxis0 = cAxis(pState, 0);
			const TypeView<const TypeSpaceVector> sAxis1 = cAxis(pState, 1);
			const TypeSpaceVector sAxis2 = TypeStepBasis::cAxis2(pState);
			// Prepare to set data
			TypeView<TypeSpaceVector> dX = x(dState.data());
			TypeView<TypeSpaceVector> dAxis0 = axis(dState.data(), 0);
			TypeView<TypeSpaceVector> dAxis1 = axis(dState.data(), 1);
			// Compute linear velocity
			dX = sFlow->getVelocity(sX, t);
			// Compute axis velocity (rotaton velocity)
			auto velocityGradients = sFlow->getVelocityGradients(sX, t);
			auto skewVelocityGradients = 0.5 * (velocityGradients - velocityGradients.transpose());
			auto symVelocityGradients = 0.5 * (velocityGradients + velocityGradients.transpose());
			// Axis 0
			dAxis0 = skewVelocityGradients * sAxis0;
			dAxis0 += factors[2] * (sAxis1 * sAxis1.transpose()) * symVelocityGradients * sAxis0;
			dAxis0 -= factors[1] * (sAxis2 * sAxis2.transpose()) * symVelocityGradients * sAxis0;
			// Axis 1
			dAxis1 = skewVelocityGradients * sAxis1;
			dAxis1 -= factors[2] * (sAxis0 * sAxis0.transpose()) * symVelocityGradients * sAxis1;
			dAxis1 += factors[0] * (sAxis2 * sAxis2.transpose()) * symVelocityGradients * sAxis1;
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
		TypeView<const TypeSpaceVector> cAxis(const double* pState, const unsigned int index) const override {
			return TypeView<const TypeSpaceVector>(pState + (index + 1) * DIM);
		}
		TypeView<TypeSpaceVector> axis(double* pState, const unsigned int index) const override {
			return TypeView<TypeSpaceVector>(pState + (index + 1) * DIM);
		}
	public:
		std::vector<TypeSpaceVector> positions(const double* pState) const override {
			return { cX(pState) };
		}
	public:
		void setDimensions(const std::vector<double>& dimensions) {
			factors.resize(dimensions.size());
			factors[0] = (std::pow(dimensions[1], 2) - std::pow(dimensions[2], 2)) / (std::pow(dimensions[2], 2) + std::pow(dimensions[1], 2));
			factors[1] = (std::pow(dimensions[2], 2) - std::pow(dimensions[0], 2)) / (std::pow(dimensions[0], 2) + std::pow(dimensions[2], 2));
			factors[2] = (std::pow(dimensions[0], 2) - std::pow(dimensions[1], 2)) / (std::pow(dimensions[1], 2) + std::pow(dimensions[0], 2));
		}
	public:
		std::shared_ptr<TypeFlow> sFlow;
		std::vector<double> factors;
};

template<template<int> typename TypeVector, template<int, int> typename TypeMatrix, unsigned int DIM, template<typename...> class TypeView, typename TypeFlow, typename TypeSolver>
class Ellipsoid : public ObjectStatic<TypeSolver, StepEllipsoid<TypeVector, TypeMatrix, DIM, TypeView, TypeFlow>> {
	public:
		using TypeStep = StepEllipsoid<TypeVector, TypeMatrix, DIM, TypeView, TypeFlow>;
		using Type = ObjectStatic<TypeSolver, TypeStep>;
	public:
		Ellipsoid(const std::shared_ptr<TypeFlow>& sFlow, const std::vector<double>& dimensions) : Type::ObjectStatic(std::make_shared<TypeStep>(sFlow, dimensions)) {
		}
	public:
		using Type::sSolver;
		using Type::sStep;
		using Type::state;
		using Type::t;
};

}

#endif
