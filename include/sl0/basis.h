#ifndef SL0_BASIS_H
#define SL0_BASIS_H
#pragma once

#include <memory>
#include <cmath>
#include <vector>

#include "sl0/object.h"

namespace sl0 {

template<template<int> typename TypeVector, template<int, int> typename TypeMatrix, unsigned int DIM, template<typename...> class TypeView>
class StepBasis : public StepObjectStatic<TypeVector, DIM, 3*DIM> {
	public:
		using Type = StepObjectStatic<TypeVector, DIM, 3*DIM>;
		using Type::StateSize;
		using typename Type::TypeSpaceVector;
		using typename Type::TypeStateVectorDynamic;
		using TypeSpaceMatrix = TypeMatrix<DIM, DIM>;
	public:
		using StepObjectStatic<TypeVector, DIM, StateSize>::StepObjectStatic;
		void update(double* pState, const double& t) override {
			StepObjectStatic<TypeVector, DIM, StateSize>::update(pState, t);
			orthoNormalizeBasis(pState);
		}
		// sub update methods
		void orthoNormalizeBasis(double* pState) {
			axis(pState, 0).normalize();
			axis(pState, 1) = (cAxis(pState, 1) - cAxis(pState, 1).dot(cAxis(pState, 0)) * cAxis(pState, 0)).normalized();
		}
		// additional getters
		TypeSpaceVector cAxis2(const double* pState) const {
   			return cAxis(pState, 0).cross(cAxis(pState, 1));
   		}
		TypeSpaceMatrix cBasis(const double* pState) const {
			TypeSpaceMatrix basis;
			basis.col(0) = cAxis(pState, 0);
			basis.col(1) = cAxis(pState, 1);
			basis.col(2) = cAxis2(pState);
   			return basis;
   		}
	public:
		virtual TypeView<const TypeSpaceVector> cAxis(const double* pState, const unsigned int index) const = 0;
		virtual TypeView<TypeSpaceVector> axis(double* pState, const unsigned int index) const = 0;
};

}

#endif
