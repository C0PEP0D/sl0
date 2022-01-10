// Std includes
#include <cmath>
#include <iostream>
#include <memory>
// Thirdparties includes
#include <Eigen/Dense>
// Lib includes
#include "s0s/runge_kutta_fehlberg.h"
#include "sl0/point.h"
#include "sl0/group/static.h"
// Simple includes
#include "flow.h"


using TypeScalar = double;
// State
template<int StateSize>
using TypeVector = Eigen::Matrix<TypeScalar, StateSize, 1>;
using TypeVectorDynamic = TypeVector<Eigen::Dynamic>;
// Space
constexpr unsigned int DIM = 3;
using TypeSpaceVector = Eigen::Matrix<TypeScalar, DIM, 1>;
// Ref and View
template<typename ...Args>
using TypeRef = Eigen::Ref<Args...>;
template<typename ...Args>
using TypeView = Eigen::Map<Args...>;
// Group Parameters
using TypeMemberStep = sl0::StepPoint<TypeVector, DIM, TypeView, Flow>;
const unsigned int MemberNb = 128; // number of objects in group
// Solver
using TypeSolver = s0s::SolverRungeKuttaFehlberg<TypeVector<-1>, TypeView>;

int main () { 
    // Parameters
    TypeSpaceVector x0 = TypeSpaceVector::Constant(1.0);
    double t0 = 0.0;
    double dt = 1e-2;
    double tEnd = 1.0;
    unsigned int nt = std::round(tEnd / dt);
    // Create group
    sl0::GroupStatic<TypeVector, DIM, TypeView, TypeMemberStep, MemberNb, TypeSolver> group(std::make_shared<TypeMemberStep>(std::make_shared<Flow>()));
    // Set initial state
    for(std::size_t n = 0; n < group.sStep->size(); n++) {
        group.sStep->sMemberStep->x(group.sStep->memberState(group.state.data(), n)) = x0;
    }
    group.t = t0;
    // Computation
    for(std::size_t i = 0; i < nt; i++) {
        group.update(dt);
    }
    // out
    std::cout << "\n";
    std::cout << "Group advected following a an exponential flow, exp(" << group.t << ") = " << "\n";
    std::cout << "\n";
    for(std::size_t n = 0; n < group.sStep->size(); n++) {
        std::cout << "Group Final Position : " << "\n" << group.sStep->sMemberStep->x(group.sStep->memberState(group.state.data(), n)) << "\n";
    }
    std::cout << std::endl;
}
