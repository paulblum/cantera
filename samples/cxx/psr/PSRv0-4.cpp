/// @file PSRv0-4.cpp
/// Perfectly Stirred Reactor Solver (v0.4)

#include "cantera/zerodim.h"

using namespace Cantera;

const char *MECHANISM = "gri30.yaml";
const double PRESSURE = Cantera::OneAtm;

int main()
{
    try
    {
        auto sol = newSolution("gri30.yaml", "gri30", "None");
        auto gas = sol->thermo();

        Reservoir inlet;
        gas->setState_TPX(300, PRESSURE, "O2:1 CH4:0.5 N2:3.76");
        inlet.insert(sol);

        Reactor combustor;
        gas->setState_TPX(300, PRESSURE, "O2:1 CH4:0.5 N2:3.76");
        gas->equilibrate("HP");
        combustor.insert(sol);
        combustor.setInitialVolume(1.0);

        Reservoir exhaust;
        exhaust.insert(sol);

        MassFlowController m1;
        m1.install(inlet, combustor);
        m1.setMassFlowRate(1000);

        Valve v;
        v.install(combustor, exhaust);
        v.setValveCoeff(1.0);

        ReactorNet sim;
        sim.addReactor(combustor);

        sim.initialize();
        sim.solveSteady();
        std::cout << combustor.temperature();
        return 0;
    }
    catch (Cantera::CanteraError &err)
    {
        std::cerr << err.what() << std::endl;
        return -1;
    }
}
