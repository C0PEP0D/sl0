// Std includes
#include <memory>
#include <vector>
#include <cmath>
#include <iostream>
// Thirdparties includes
#include <Eigen/Dense>
// Lib includes
#include "s0s/solver.h"
#include "sl0/ellipsoid.h"
#include "sa0/actuators.h"
#include "sa0/behaviours.h"
#include "sa0/sensors/local.h"
#include "sa0/sensors/computer.h"
#include "sa0/agents.h"
#include "fl0w/simple_shear.h"

using TypeScalar = double;
using TypeVector = Eigen::Matrix<TypeScalar, 3, 1>;
using TypeMatrix = Eigen::Matrix<TypeScalar, 3, 3>;
using TypeState = Eigen::Matrix<TypeScalar, 12, 1>;
template<typename ...Args>
using TypeView = Eigen::Map<Args...>;
template<typename A, template<typename ...Args> class B, typename C>
using TypeSolver = s0s::SolverRungeKuttaFehlberg<A, B, C>;
using TypeFlow = fl0w::SimpleShear<TypeVector, TypeMatrix>;

// Agent
using TypeStepObject = sl0::StepEllipsoid<TypeState, TypeView, TypeVector, TypeMatrix, TypeFlow>;
using TypeObject = sl0::ObjectBasis<TypeState, TypeView, TypeSolver, sl0::sa0::StepActuators<TypeState, TypeStepObject, std::vector>>;

using TypeObserver = sl0::sa0::Observer<std::vector<double>, TypeState, TypeStepObject>;
using TypeBehaviour = sl0::sa0::BehaviourLinearOnOff<std::vector<double>, std::vector<double>>;
using TypeAgent = sl0::sa0::ObjectActive<TypeObject, TypeBehaviour, TypeObserver, std::vector>;
using TypeStepSwim = sl0::sa0::StepSwim<TypeState, TypeStepObject, TypeVector>;

using TypeSensorVector = sl0::sa0::Sensor<TypeVector, TypeState, TypeStepObject, TypeVector>;
using TypeSensorMatrix = sl0::sa0::Sensor<TypeMatrix, TypeState, TypeStepObject, TypeVector>;
// Local sensors
using TypeSensorDirection = sl0::sa0::SensorLocalVector<sl0::sa0::SensorDirection, TypeVector, std::vector<double>, TypeState, TypeStepObject>;
using TypeSensorJacobian = sl0::sa0::SensorLocalMatrix<sl0::sa0::SensorJacobian, TypeVector, TypeMatrix, std::vector<double>, TypeState, TypeStepObject, TypeVector, TypeFlow>;
using TypeSensorStrain = sl0::sa0::SensorLocalMatrix<sl0::sa0::SensorStrain, TypeVector, TypeMatrix, std::vector<double>, TypeState, TypeStepObject, TypeVector, TypeFlow>;
// Computers
using TypeSensorTranspose = sl0::sa0::SensorTranspose<TypeMatrix, std::vector<double>, TypeState, TypeStepObject>;
using TypeSensorMult = sl0::sa0::SensorMult<TypeVector, std::vector<double>, TypeState, TypeStepObject, TypeMatrix>;
using TypeSensorAxis = sl0::sa0::SensorAxis<std::vector<double>, TypeState, TypeStepObject, TypeVector>;

class Agent : public TypeAgent {
    public:
        using TypeAgent::TypeAgent;
    public:
        void init() override {
            TypeAgent::init();
            // Sensors
            // // Direction
            sSensorDir = std::make_shared<TypeSensorDirection>();
            sSensorDir->pointSensor.direction << 1.0, 0.0, 0.0;
            // // Jacobian
            sSensorJac = std::make_shared<TypeSensorJacobian>();
            sSensorJac->position << 0.0, 0.0, 0.0;
            sSensorJac->pointSensor.sFlow = TypeObject::step.sFlow;
            // // Transpose
            sSensorJacTranspose = std::make_shared<TypeSensorTranspose>();
            sSensorJacTranspose->sSensor = sSensorJac;
            // // Gradient
            sSensorGrad = std::make_shared<TypeSensorMult>();
            sSensorGrad->sA = sSensorJacTranspose;
            sSensorGrad->sB = sSensorDir;
            // // Strain
            sSensorStrain = std::make_shared<TypeSensorStrain>();
            sSensorStrain->position << 0.0, 0.0, 0.0;
            sSensorStrain->pointSensor.sFlow = TypeObject::step.sFlow;
            sSensorStrainTranspose = std::make_shared<TypeSensorTranspose>();
            sSensorStrainTranspose->sSensor = sSensorStrain;
            sSensorStrainGrad = std::make_shared<TypeSensorMult>();
            sSensorStrainGrad->sA = sSensorStrainTranspose;
            sSensorStrainGrad->sB = sSensorDir;
            // Observers
            sObserveDir = std::make_shared<TypeSensorAxis>();
            sObserveDir->axis = 0;
            sObserveDir->sSensor = sSensorDir;
            sObserveGrad = std::make_shared<TypeSensorAxis>();
            sObserveGrad->axis = 0;
            sObserveGrad->sSensor = sSensorGrad;
            // Add to behaviour observers
            observers.push_back(sObserveDir);
            observers.push_back(sObserveGrad);
        }
    public:
        // Sensors
        std::shared_ptr<TypeSensorDirection> sSensorDir;
        std::shared_ptr<TypeSensorJacobian> sSensorJac;
        std::shared_ptr<TypeSensorTranspose> sSensorJacTranspose;
        std::shared_ptr<TypeSensorMult> sSensorGrad;
        std::shared_ptr<TypeSensorStrain> sSensorStrain;
        std::shared_ptr<TypeSensorTranspose> sSensorStrainTranspose;
        std::shared_ptr<TypeSensorMult> sSensorStrainGrad;
        // Observers
        std::shared_ptr<TypeSensorAxis> sObserveDir;
        std::shared_ptr<TypeSensorAxis> sObserveGrad;
};

int main () { 
    TypeVector x0;
    x0 << 0,
          0,
          0;
    TypeMatrix b0;
    b0 << 1, 0, 0,
          0, 1, 0,
          0, 0, 1;
    double t0 = 0.0;
    double dt = 1e-2;
    // Building the agent
    Agent agent;
    agent.step.setProportions(TypeVector(3.0, 2.0, 1.0));
    agent.step.x(agent.state) = x0;
    agent.step.basis(agent.state) = b0;
    agent.solver.t = t0;
    agent.solver.dt = dt;
    //// Building flow
    agent.step.sFlow = std::make_shared<TypeFlow>();
    agent.step.sFlow->gamma = 1.0;
    //// Building Swim
    std::shared_ptr<TypeStepSwim> sStepSwim = std::make_shared<TypeStepSwim>();
    sStepSwim->velocity << 0.0,
                           0.0,
                           0.0;
    sStepSwim->intensity = 0.0;
    agent.step.stepActuators.push_back(sStepSwim);
    //// Building behaviour
    agent.behaviour.threshold = 0.0;
    agent.behaviour.brain = {0.0, 0.0};
    // Init the agent
    agent.init();
    // Compute
    std::cout << std::endl;
    std::cout << "Start : " << std::endl;
    for(std::size_t i = 0; i < 1000; i++) {
        agent.update();
    }
    std::cout << "Dir agent local " << std::endl;
    std::cout << agent.sSensorDir->operator()(agent.step, agent.state, agent.solver.t) << std::endl;
    std::cout << "Dir agent global " << std::endl;
    std::cout << agent.step.basis(agent.state).transpose() * agent.sSensorDir->operator()(agent.step, agent.state, agent.solver.t) << std::endl;
    std::cout << "Dir agent point " << std::endl;
    std::cout << agent.sSensorDir->pointSensor(agent.step, agent.state, agent.solver.t) << std::endl;
    std::cout << "Dir true  global " << std::endl;
    std::cout << agent.sSensorDir->pointSensor.direction << std::endl;
    std::cout << "Jac agent local " << std::endl;
    std::cout << agent.sSensorJac->operator()(agent.step, agent.state, agent.solver.t) << std::endl;
    std::cout << "Jac agent global " << std::endl;
    std::cout << agent.step.basis(agent.state).transpose() * agent.sSensorJac->operator()(agent.step, agent.state, agent.solver.t) * agent.step.basis(agent.state) << std::endl;
    std::cout << "Jac agent point " << std::endl;
    std::cout << agent.sSensorJac->pointSensor(agent.step, agent.state, agent.solver.t) << std::endl;
    std::cout << "Jac true  global " << std::endl;
    std::cout << agent.step.sFlow->getJacobian(agent.step.x(agent.state), agent.solver.t) << std::endl;
    std::cout << "Grad agent local " << std::endl;
    std::cout << agent.sSensorGrad->operator()(agent.step, agent.state, agent.solver.t) << std::endl;
    std::cout << "Grad agent global " << std::endl;
    std::cout << agent.step.basis(agent.state).transpose() * agent.sSensorGrad->operator()(agent.step, agent.state, agent.solver.t) << std::endl;
    std::cout << "Grad true  global " << std::endl;
    std::cout << agent.step.sFlow->getJacobian(agent.step.x(agent.state), agent.solver.t).transpose().col(0) << std::endl;
    std::cout << "Strain agent local " << std::endl;
    std::cout << agent.sSensorStrain->operator()(agent.step, agent.state, agent.solver.t) << std::endl;
    std::cout << "Strain agent global " << std::endl;
    std::cout << agent.step.basis(agent.state).transpose() * agent.sSensorStrain->operator()(agent.step, agent.state, agent.solver.t) * agent.step.basis(agent.state) << std::endl;
    std::cout << "Strain agent point " << std::endl;
    std::cout << agent.sSensorStrain->pointSensor(agent.step, agent.state, agent.solver.t) << std::endl;
    std::cout << "Strain true  global " << std::endl;
    std::cout << agent.step.sFlow->getStrain(agent.step.x(agent.state), agent.solver.t) << std::endl;
    std::cout << "Grad Strain agent local " << std::endl;
    std::cout << agent.sSensorStrainGrad->operator()(agent.step, agent.state, agent.solver.t) << std::endl;
    std::cout << "Grad Strain agent global " << std::endl;
    std::cout << agent.step.basis(agent.state).transpose() * agent.sSensorStrainGrad->operator()(agent.step, agent.state, agent.solver.t) << std::endl;
    // out
    std::cout << std::endl;
    std::cout << "Final Position : " << std::endl;
    std::cout << agent.step.x(agent.state) << std::endl;
    std::cout << std::endl;
}
