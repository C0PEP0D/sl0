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
// Group Parameters
using TypeStepPoint = sl0::StepPoint<TypeVector, DIM, TypeView, Flow>;
// Solver
using TypeSolver = s0s::SolverRungeKuttaFehlberg<TypeVector<Eigen::Dynamic>, TypeView>;

int main () { 
    // Parameters
    TypeSpaceVector x0 = TypeSpaceVector::Constant(1.0);
    TypeScalar t0 = 0.0;
    TypeScalar dt = 1e-3;
    unsigned int nt = std::round(1.0 / dt);
    unsigned int np = 1000;
    // Create group
    sl0::GroupDynamic<TypeVector, DIM, TypeView, TypeSolver> group;
    // Ceate points
    std::cout << "Point creation" << "\n";
    std::shared_ptr<TypeStepPoint> sStepPoint = std::make_shared<TypeStepPoint>(std::make_shared<Flow>());
    for(std::size_t i = 0; i < np; i++) {
        group.sStep->addMember(group.state, std::make_shared<TypeStepPoint>(*sStepPoint), DIM);
        sStepPoint->x(group.sStep->memberState(group.state.data(), i)) = x0;
    }
    // Set initial state
    std::cout << "Point init" << "\n";
    for(std::size_t n = 0; n < group.sStep->size(); n++) {
        sStepPoint->x(group.sStep->memberState(group.state.data(), n)) = x0;
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
        std::cout << "Group Final Position : " << "\n" << sStepPoint->x(group.sStep->memberState(group.state.data(), n)) << "\n";
    }
    std::cout << std::endl;
}
