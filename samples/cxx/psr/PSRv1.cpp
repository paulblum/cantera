/// @file PSRv1.cpp
/// Perfectly Stirred Reactor (v0.1)

#include "BoundaryValueProblem.h"

const char *MECHANISM = "gri30.yaml";
const double PRESSURE = Cantera::OneAtm;

class PSR : public BVP::BoundaryValueProblem
{
public:
    std::shared_ptr<Cantera::ThermoPhase> inletThermo;
    std::shared_ptr<Cantera::ThermoPhase> reactorThermo;
    std::shared_ptr<Cantera::Kinetics> reactorKinetics;

    const double *initialVector;
    double mdot;
    double volume;

    // constructor
    PSR(std::shared_ptr<Cantera::Solution> inletSol, std::shared_ptr<Cantera::Solution> reactorSol, double mdot_input, double volume_input)
        : BVP::BoundaryValueProblem(inletSol->thermo()->nSpecies(), 1, 0.0, 1)
    {
        inletThermo = inletSol->thermo();
        reactorThermo = reactorSol->thermo();
        reactorKinetics = reactorSol->kinetics();
        initialVector = reactorSol->thermo()->massFractions();
        mdot = mdot_input;
        volume = volume_input;

        BVP::Component initializer;
        initializer.lower = 0;
        initializer.upper = 1;
        initializer.rtol = 1.0e-12;
        initializer.atol = 1.0e-15;

        for (int i = 0; i < inletThermo->nSpecies(); i++)
        {
            initializer.name = inletThermo->speciesName(i);
            setComponent(i, initializer);
        }
    }

    // destructor
    virtual ~PSR() {}

    // specify guesses for the initial values. These can be anything
    // that leads to a converged solution.
    virtual doublereal initialValue(size_t n, size_t j)
    {
        return initialVector[n];
    }

    // Specify the residual function. This is where the ODE system and boundary
    // conditions are specified. The solver will attempt to find a solution
    // x so that rsd is zero.
    void eval(size_t jg, double *x, double *rsd, int *diag, double rdt)
    {
        size_t j = 0;

        reactorThermo->setMassFractions_NoNorm(x);
        reactorThermo->setState_HP(inletThermo->enthalpy_mass(), PRESSURE);

        Cantera::vector_fp wdot(inletThermo->nSpecies());
        reactorKinetics->getNetProductionRates(wdot.data());

        for (int i = 0; i < inletThermo->nSpecies(); i++)
        {
            rsd[index(i, j)] = wdot.at(i) * inletThermo->molecularWeight(i) * volume + mdot * (inletThermo->massFraction(i) - value(x, i, j));
        }
    }
};

int main()
{
    try
    {
        double mdot;
        double volume;

        std::cout << "Enter mdot: ";
        std::cin >> mdot;
        std::cout << "\nEnter volume: ";
        std::cin >> volume;

        std::shared_ptr<Cantera::Solution> mixInlet(Cantera::newSolution(MECHANISM));
        std::shared_ptr<Cantera::Solution> mixOutlet(Cantera::newSolution(MECHANISM));

        mixInlet->thermo()->setState_TPX(298, PRESSURE, "H2:2 O2:1");
        mixOutlet->thermo()->setState_TPX(298, PRESSURE, "H2:2 O2:1");
        mixOutlet->thermo()->equilibrate("HP");

        PSR reactor(mixInlet, mixOutlet, mdot, volume);

        reactor.solve();

        std::cout << mixOutlet->thermo()->report() << "\n";

        return 0;
    }
    catch (Cantera::CanteraError &err)
    {
        std::cerr << err.what() << std::endl;
        return -1;
    }
}
