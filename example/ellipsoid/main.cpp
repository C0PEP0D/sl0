// Std includes
#include <cmath>
#include <iostream>
// Thirdparties includes
#include <Eigen/Dense>
// Lib includes
#include "Eigen/src/Core/Matrix.h"
#include "s0s/runge_kutta_fehlberg.h"
#include "sl0/ellipsoid.h"
// Simple includes
#include "flow.h"

// Types
using TypeScalar = double;
template<int Size>
using TypeVector = Eigen::Matrix<TypeScalar, Size, 1>;
template<int Nx, int Ny>
using TypeMatrix = Eigen::Matrix<TypeScalar, Nx, Ny>;
// Space
constexpr unsigned int DIM = 3;
using TypeSpaceVector = Eigen::Matrix<TypeScalar, DIM, 1>;
using TypeSpaceMatrix = Eigen::Matrix<TypeScalar, DIM, DIM>;
// Ref and View
template<typename ...Args>
using TypeRef = Eigen::Ref<Args...>;
template<typename ...Args>
using TypeView = Eigen::Map<Args...>;
// Solver
using TypeSolver = s0s::SolverRungeKuttaFehlberg<TypeVector<Eigen::Dynamic>, TypeView>;
// Flow
using TypeFlow = Flow<TypeSpaceVector, TypeSpaceMatrix, TypeRef>;

int main () { 
	TypeSpaceVector x0;
	x0 << 0,
		  0,
		  0;
	TypeSpaceVector p0;
	p0 << 1,
		  0,
		  0;
	TypeSpaceVector p1;
	p1 << 0,
		  1,
		  0;
	double t0 = 0.0;
	double dt = 1e-2;
	double tEnd = 1.0;
	// Create ellipsoid
	sl0::Ellipsoid<TypeVector, TypeMatrix, DIM, TypeView, TypeFlow, TypeSolver> ellipsoid(std::make_shared<TypeFlow>(), std::vector<double>({1.0, 1.0, 1.0}));
	// Set initial state
	ellipsoid.sStep->x(ellipsoid.state.data()) = x0;
	ellipsoid.sStep->axis(ellipsoid.state.data(), 0) = p0;
	ellipsoid.sStep->axis(ellipsoid.state.data(), 1) = p1;
	ellipsoid.t = t0;
	// Compute
	std::cout << "Ellipsoid in a simple shear flow : \n";
	std::cout << "\n";
	std::cout << "Orientation after each step : \n";
	for(std::size_t i = 0; i < (tEnd - t0)/dt; i++) {
		ellipsoid.update(dt);
		std::cout << ellipsoid.sStep->cBasis(ellipsoid.state.data()) << "\n";
	}
	// out
	std::cout << "Ellipsoid in a simple shear flow : \n";
	std::cout << "\n";
	std::cout << "Final orientation : \n";
	std::cout << ellipsoid.sStep->cBasis(ellipsoid.state.data()) << "\n";
	std::cout << std::endl;
}
