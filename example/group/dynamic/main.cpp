// Std includes
#include <cmath>
#include <iostream>
#include <memory>
// Thirdparties includes
#include <Eigen/Dense>
// Lib includes
#include "s0s/runge_kutta_fehlberg.h"
#include "sl0/point.h"
#include "sl0/group/dynamic.h"
// Simple includes
#include "flow.h"

using TypeScalar = double;
// State
template<int StateSize>
using TypeState = Eigen::Matrix<TypeScalar, StateSize, 1>;
using TypeStateDynamic = Eigen::Matrix<TypeScalar, Eigen::Dynamic, 1>;
// Space
constexpr unsigned int DIM = 3;
using TypeVector = Eigen::Matrix<TypeScalar, DIM, 1>;
// Ref and View
template<typename ...Args>
using TypeRef = Eigen::Ref<Args...>;
template<typename ...Args>
using TypeView = Eigen::Map<Args...>;
// Group Parameters
using TypeStepPoint = sl0::StepPoint<TypeState, DIM, TypeRef, TypeView, Flow>;
// Solver
using TypeSolver = s0s::SolverRungeKuttaFehlberg;

int main () { 
    // Parameters
    TypeVector x0 = TypeVector::Constant(1.0);
    TypeScalar t0 = 0.0;
    TypeScalar dt = 1e-3;
    unsigned int nt = std::round(1.0 / dt);
    unsigned int np = 1000;
    // Create group
    sl0::GroupDynamic<TypeState, DIM, TypeRef, TypeView, sl0::StepObject<TypeState, DIM, TypeRef>, TypeSolver> group;
    // Ceate points
    std::cout << "Point creation" << "\n";
    std::shared_ptr<TypeStepPoint> sStepPoint = std::make_shared<TypeStepPoint>(std::make_shared<Flow>());
    for(std::size_t i = 0; i < np; i++) {
        group.addMember(std::make_shared<TypeStepPoint>(*sStepPoint), DIM);
        sStepPoint->x(group.sStep->memberState(group.state, i)) = x0;
    }
    // Set initial state
    std::cout << "Point init" << "\n";
    for(std::size_t n = 0; n < group.sStep->size(); n++) {
        sStepPoint->x(group.sStep->memberState(group.state, n)) = x0;
    }
    group.t = t0;
    // Computation
    std::cout << "Computing" << "\n";
    for(std::size_t i = 0; i < nt; i++) {
        group.update(dt);
    }
    // out
    std::cout << "\n";
    std::cout << "Group advected following a an exponential flow, exp(" << group.t << ") = " << "\n";
    std::cout << "\n";
    for(std::size_t n = 0; n < group.sStep->size(); n++) {
        std::cout << "Group Final Position : " << "\n" << sStepPoint->x(group.sStep->memberState(group.state, n)) << "\n";
    }
    std::cout << std::endl;
}
