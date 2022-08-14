#ifndef SL0_POINT_INERTIAL_KINEMATIC_H
#define SL0_POINT_INERTIAL_KINEMATIC_H
#pragma once

// std includes
#include <memory> // shared_ptr
#include <vector>
#include <unsupported/Eigen/MatrixFunctions>
// module includes
#include "sl0/object.h"

namespace sl0 {

template<template<int> typename TypeVector, unsigned int DIM, template<typename...> class TypeView, typename TypeFlow>
class StepPointInertialKinematic : public StepObjectStatic<TypeVector, DIM, DIM> {
	public:
		using Type = StepObjectStatic<TypeVector, DIM, DIM>;
		using Type::StateSize;
		using typename Type::TypeSpaceVector;
		using typename Type::TypeStateVectorDynamic;
	public:
		StepPointInertialKinematic(const std::shared_ptr<TypeFlow>& p_sFlow, const double p_delay, const TypeSpaceVector p_velocity) : sFlow(p_sFlow), delay(p_delay), velocity(p_velocity) {
		}

		TypeStateVectorDynamic operator()(const double* pState, const double& t) const override {
			TypeStateVectorDynamic dState(Type::stateSize());
			TypeView<TypeSpaceVector> dX = x(dState.data());
			dX = cU(pState, t);
			return dState;
		}
	public:
		TypeView<const TypeSpaceVector> cX(const double* pState) const {
			return TypeView<const TypeSpaceVector>(pState);
		}
		TypeView<TypeSpaceVector> x(double* pState) const {
			return TypeView<TypeSpaceVector>(pState);
		}

		TypeSpaceVector cU(const double* pState, const double& t) const {
			//return (TypeSpaceMatrix::Identity() + delay * sFlow->getVelocityGradients(cX(pState), t)).inverse() * (sFlow->getVelocity(cX(pState), t) + velocity - delay * sFlow->getAcceleration(cX(pState), t));
			const auto velocityGradients = sFlow->getVelocityGradients(cX(pState), t);
			const auto identity = velocityGradients.Identity();
			return (identity + delay * velocityGradients).inverse() * (sFlow->getVelocity(cX(pState), t) + velocity);
		}
	public:
		std::vector<TypeSpaceVector> positions(const double* pState) const override {
			return { cX(pState) };
		}
	public:
		std::shared_ptr<TypeFlow> sFlow;
		double delay;
		TypeSpaceVector velocity;
};

template<template<int> typename TypeVector, unsigned int DIM, template<typename...> class TypeView, typename TypeFlow, typename TypeSolver>
class PointInertialKinematic : public ObjectStatic<TypeSolver, StepPointInertialKinematic<TypeVector, DIM, TypeView, TypeFlow>> {
	public:
		using TypeStep = StepPointInertialKinematic<TypeVector, DIM, TypeView, TypeFlow>;
		using Type = ObjectStatic<TypeSolver, TypeStep>;
		using typename Type::TypeSpaceVector;
	public:
		PointInertialKinematic(const std::shared_ptr<TypeFlow>& sFlow, const double delay, const TypeSpaceVector velocity) : ObjectStatic<TypeSolver, TypeStep>::ObjectStatic(std::make_shared<TypeStep>(sFlow, delay, velocity)) {
		}
	public:
		using Type::sSolver;
		using Type::sStep;
		using Type::state;
		using Type::t;
};

}

#endif
