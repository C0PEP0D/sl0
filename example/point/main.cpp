// Std includes
#include <cmath>
#include <iostream>
// Thirdparties includes
#include <Eigen/Dense>
// Lib includes
#include "s0s/runge_kutta_fehlberg.h"
#include "sl0/point.h"
// Simple includes
#include "flow.h"

using TypeScalar = double;
// Linear Algebra
template<int Size>
using TypeVector = Eigen::Matrix<TypeScalar, Size, 1>;
// Space
constexpr unsigned int DIM = 3;
using TypeSpaceVector = Eigen::Matrix<TypeScalar, DIM, 1>;
// Ref and View
template<typename ...Args>
using TypeRef = Eigen::Ref<Args...>;
template<typename ...Args>
using TypeView = Eigen::Map<Args...>;
// Solver
using TypeSolver = s0s::SolverRungeKuttaFehlberg<TypeVector<Eigen::Dynamic>, TypeView>;
// Flow
using TypeFlow = Flow<TypeSpaceVector, TypeRef>;

int main () { 
    TypeSpaceVector x0 = TypeSpaceVector::Constant(1.0);
    double t0 = 0.0;
    double dt = 1e-3;
    double tEnd = 1.0;
    unsigned int nt = std::round((tEnd - t0) / dt);
    // Create point
    sl0::Point<TypeVector, DIM, TypeView, TypeFlow, TypeSolver> point(std::make_shared<TypeFlow>());
    // Set initial state
    point.sStep->x(point.state.data()) = x0;
    point.t = t0;
    // Computation
    for(std::size_t i = 0; i < nt; i++) {
        point.update(dt);
    }
    // out
    std::cout << "\n";
    std::cout << "Point advected following a an exponential flow, exp(" << point.t << ") = " << "\n";
    std::cout << "\n";
    std::cout << "Point Final Position : " << "\n" << point.sStep->x(point.state.data()) << "\n";
    std::cout << std::endl;
}
