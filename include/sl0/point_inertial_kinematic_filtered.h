#ifndef SL0_POINT_INERTIAL_KINEMATIC_FILTERED_H
#define SL0_POINT_INERTIAL_KINEMATIC_FILTERED_H
#pragma once

// std includes
#include <memory> // shared_ptr
#include <vector>
#include <unsupported/Eigen/MatrixFunctions>
// module includes
#include "sl0/object.h"

namespace sl0 {

template<template<int> typename TypeVector, template<int, int> typename TypeMatrix, unsigned int DIM, template<typename...> class TypeView, typename TypeFlow>
class StepPointInertialKinematicFiltered : public StepObjectStatic<TypeVector, DIM, DIM + DIM> {
	public:
		using Type = StepObjectStatic<TypeVector, DIM, DIM + DIM>;
		using Type::StateSize;
		using typename Type::TypeSpaceVector;
		using typename Type::TypeStateVectorDynamic;
	public:
		StepPointInertialKinematicFiltered(const std::shared_ptr<TypeFlow>& p_sFlow, const double p_delay, const TypeSpaceVector p_velocity, const double p_dt) : sFlow(p_sFlow), delay(p_delay), velocity(p_velocity), dt(p_dt) {
		}

		TypeStateVectorDynamic operator()(const double* pState, const double& t) const override {
			TypeStateVectorDynamic dState(Type::stateSize());
			TypeView<TypeSpaceVector> dX = x(dState.data());
			TypeView<TypeSpaceVector> dFilteredVelocity = filteredVelocity(dState.data());
			// set dstate
			dX = cU(pState, t);
			dFilteredVelocity = (sFlow->getVelocity(cX(pState), t) - cFilteredVelocity(pState)) / delay;
			// return dstate
			return dState;
		}
	public:
		TypeView<const TypeSpaceVector> cX(const double* pState) const {
			return TypeView<const TypeSpaceVector>(pState);
		}
		TypeView<TypeSpaceVector> x(double* pState) const {
			return TypeView<TypeSpaceVector>(pState);
		}

		TypeView<const TypeSpaceVector> cFilteredVelocity(const double* pState) const {
			return TypeView<const TypeSpaceVector>(pState + DIM);
		}
		TypeView<TypeSpaceVector> filteredVelocity(double* pState) const {
			return TypeView<TypeSpaceVector>(pState + DIM);
		}

		TypeSpaceVector cU(const double* pState, const double& t) const {
			return velocity + cFilteredVelocity(pState);
		}
	public:
		std::vector<TypeSpaceVector> positions(const double* pState) const override {
			return { cX(pState) };
		}
	public:
		std::shared_ptr<TypeFlow> sFlow;
		double delay;
		TypeSpaceVector velocity;
		double dt; // used to compute eulerian acceleartion

		// TODO:
		// If this is kept remove unused TypeMatrix and dt
};

// template<template<int> typename TypeVector, template<int, int> typename TypeMatrix, unsigned int DIM, template<typename...> class TypeView, typename TypeFlow>
// class StepPointInertialKinematicFiltered : public StepObjectStatic<TypeVector, DIM, DIM + DIM * DIM + DIM> {
	// public:
		// using Type = StepObjectStatic<TypeVector, DIM, DIM + DIM * DIM + DIM>;
		// using Type::StateSize;
		// using typename Type::TypeSpaceVector;
		// using typename Type::TypeStateVectorDynamic;
		// using TypeSpaceMatrix = TypeMatrix<DIM, DIM>;
	// public:
		// StepPointInertialKinematicFiltered(const std::shared_ptr<TypeFlow>& p_sFlow, const double p_delay, const TypeSpaceVector p_velocity, const double p_dt) : sFlow(p_sFlow), delay(p_delay), velocity(p_velocity), dt(p_dt) {
		// }
// 
		// TypeStateVectorDynamic operator()(const double* pState, const double& t) const override {
			// TypeStateVectorDynamic dState(Type::stateSize());
			// TypeView<TypeSpaceVector> dX = x(dState.data());
			// TypeView<TypeSpaceMatrix> dFilteredVelocityGradients = filteredVelocityGradients(dState.data());
			// TypeView<TypeSpaceVector> dFilteredEulerianAcceleration = filteredEulerianAcceleration(dState.data());
			// // set dstate
			// dX = cU(pState, t);
			// dFilteredVelocityGradients = (sFlow->getVelocityGradients(cX(pState), t) - cFilteredVelocityGradients(pState)) / delay;
			// dFilteredEulerianAcceleration = ((sFlow->getVelocity(cX(pState), t + 0.5 * dt) - sFlow->getVelocity(cX(pState), t - 0.5 * dt)) / dt - cFilteredEulerianAcceleration(pState)) / delay;
			// // return dstate
			// return dState;
		// }
	// public:
		// TypeView<const TypeSpaceVector> cX(const double* pState) const {
			// return TypeView<const TypeSpaceVector>(pState);
		// }
		// TypeView<TypeSpaceVector> x(double* pState) const {
			// return TypeView<TypeSpaceVector>(pState);
		// }
// 
		// TypeView<const TypeSpaceMatrix> cFilteredVelocityGradients(const double* pState) const {
			// return TypeView<const TypeSpaceMatrix>(pState + DIM);
		// }
		// TypeView<TypeSpaceMatrix> filteredVelocityGradients(double* pState) const {
			// return TypeView<TypeSpaceMatrix>(pState + DIM);
		// }
// 
		// TypeView<const TypeSpaceVector> cFilteredEulerianAcceleration(const double* pState) const {
			// return TypeView<const TypeSpaceVector>(pState + DIM + DIM * DIM);
		// }
		// TypeView<TypeSpaceVector> filteredEulerianAcceleration(double* pState) const {
			// return TypeView<TypeSpaceVector>(pState + DIM + DIM * DIM);
		// }
// 
		// TypeSpaceVector cU(const double* pState, const double& t) const {
			// const auto identity = TypeSpaceMatrix::Identity();
			// return (identity + delay * cFilteredVelocityGradients(pState)).colPivHouseholderQr().solve(sFlow->getVelocity(cX(pState), t) + velocity + delay * cFilteredEulerianAcceleration(pState));
		// }
	// public:
		// std::vector<TypeSpaceVector> positions(const double* pState) const override {
			// return { cX(pState) };
		// }
	// public:
		// std::shared_ptr<TypeFlow> sFlow;
		// double delay;
		// TypeSpaceVector velocity;
		// double dt; // used to compute eulerian acceleartion
// };

template<template<int> typename TypeVector, template<int, int> typename TypeMatrix, unsigned int DIM, template<typename...> class TypeView, typename TypeFlow, typename TypeSolver>
class PointInertialKinematicFiltered : public ObjectStatic<TypeSolver, StepPointInertialKinematicFiltered<TypeVector, TypeMatrix, DIM, TypeView, TypeFlow>> {
	public:
		using TypeStep = StepPointInertialKinematicFiltered<TypeVector, TypeMatrix, DIM, TypeView, TypeFlow>;
		using Type = ObjectStatic<TypeSolver, TypeStep>;
		using typename Type::TypeSpaceVector;
	public:
		PointInertialKinematicFiltered(const std::shared_ptr<TypeFlow>& sFlow, const double delay, const TypeSpaceVector velocity, const double dt) : ObjectStatic<TypeSolver, TypeStep>::ObjectStatic(std::make_shared<TypeStep>(sFlow, delay, velocity, dt)) {
		}
	public:
		using Type::sSolver;
		using Type::sStep;
		using Type::state;
		using Type::t;
};

}

#endif
