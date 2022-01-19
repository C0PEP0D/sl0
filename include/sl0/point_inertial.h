#ifndef SL0_POINT_INERTIAL_H
#define SL0_POINT_INERTIAL_H
#pragma once

// std includes
#include <memory> // shared_ptr
#include <vector>
// module includes
#include "sl0/object.h"

namespace sl0 {

template<template<int> typename TypeVector, unsigned int DIM, template<typename...> class TypeView, typename TypeFlow>
class StepPointInertial : public StepObjectStatic<TypeVector, DIM, 2*DIM> {
	public:
		using Type = StepObjectStatic<TypeVector, DIM, 2*DIM>;
		using Type::StateSize;
		using typename Type::TypeSpaceVector;
		using typename Type::TypeStateVectorDynamic;
	public:
		StepPointInertial(const std::shared_ptr<TypeFlow>& p_sFlow, const double p_delay) : sFlow(p_sFlow), delay(p_delay) {
		}

		TypeStateVectorDynamic operator()(const double* pState, const double& t) const override {
			TypeStateVectorDynamic dState(Type::stateSize());
			TypeView<TypeSpaceVector> dX = x(dState.data());
			TypeView<TypeSpaceVector> dU = u(dState.data());
			dX = cU(pState);
			dU = (sFlow->getVelocity(cX(pState), t) - cU(pState)) / delay;
			return dState;
		}
	public:
		TypeView<const TypeSpaceVector> cX(const double* pState) const {
			return TypeView<const TypeSpaceVector>(pState);
		}
		TypeView<TypeSpaceVector> x(double* pState) const {
			return TypeView<TypeSpaceVector>(pState);
		}

		TypeView<const TypeSpaceVector> cU(const double* pState) const {
			return TypeView<const TypeSpaceVector>(pState + DIM);
		}
		TypeView<TypeSpaceVector> u(double* pState) const {
			return TypeView<TypeSpaceVector>(pState + DIM);
		}
	public:
		std::vector<TypeSpaceVector> positions(const double* pState) const override {
			return { cX(pState) };
		}
	public:
		std::shared_ptr<TypeFlow> sFlow;
		double delay;
};

template<template<int> typename TypeVector, unsigned int DIM, template<typename...> class TypeView, typename TypeFlow, typename TypeSolver>
class PointInertial : public ObjectStatic<TypeSolver, StepPointInertial<TypeVector, DIM, TypeView, TypeFlow>> {
	public:
		using TypeStep = StepPointInertial<TypeVector, DIM, TypeView, TypeFlow>;
		using Type = ObjectStatic<TypeSolver, TypeStep>;
	public:
		PointInertial(const std::shared_ptr<TypeFlow>& sFlow, const double delay) : ObjectStatic<TypeSolver, TypeStep>::ObjectStatic(std::make_shared<TypeStep>(sFlow, delay)) {
		}
	public:
		using Type::sSolver;
		using Type::sStep;
		using Type::state;
		using Type::t;
};

}

#endif
