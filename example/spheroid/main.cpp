// Std includes
#include <cmath>
#include <iostream>
// Thirdparties includes
#include <Eigen/Dense>
// Lib includes
#include "Eigen/src/Core/Matrix.h"
#include "s0s/runge_kutta_fehlberg.h"
#include "sl0/spheroid.h"
// Simple includes
#include "flow.h"

// Types
using TypeScalar = double;
template<int Size>
using TypeVector = Eigen::Matrix<TypeScalar, Size, 1>;
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
using TypeSolver = s0s::SolverRungeKuttaFehlberg;
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
    double t0 = 0.0;
    double dt = 1e-2;
    double tEnd = 1.0;
    // Create ellipsoid
    sl0::Spheroid<TypeVector, DIM, TypeRef, TypeView, TypeFlow, TypeSolver> spheroid(std::make_shared<TypeFlow>(), 1.0);
    // Set initial state
    spheroid.sStep->x(spheroid.state) = x0;
    spheroid.sStep->axis(spheroid.state) = p0;
    spheroid.t = t0;
    // Compute
    std::cout << "Spheroid in a simple shear flow : \n";
    std::cout << "\n";
    std::cout << "Orientation after each step : \n";
    for(std::size_t i = 0; i < (tEnd - t0)/dt; i++) {
        spheroid.update(dt);
        std::cout << spheroid.sStep->axis(spheroid.state) << "\n";
    }
    // out
    std::cout << "Spheroid in a simple shear flow : \n";
    std::cout << "\n";
    std::cout << "Final orientation : \n";
    std::cout << spheroid.sStep->axis(spheroid.state) << "\n";
    std::cout << std::endl;
}
