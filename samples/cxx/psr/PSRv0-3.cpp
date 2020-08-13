/// @file PSRv0-3.cpp
/// Perfectly Stirred Reactor Solver (v0.3)
/// Changes since v0.2.1:
///     - Created a separate module for nonlinear system solution (Cantera1D_NonLinSol)
///     - Removed nonlinear solver capability from PSR

#include "cantera/numerics/Cantera_NonLinSol.h" //a nonlinear algebraic system solver

const char *MECHANISM = "gri30.yaml";
const double PRESSURE = Cantera::OneAtm;

class PSR : public Cantera_NonLinSol //inherit nonlinear solver functionality
{
public:
    // Constructor
    PSR(std::shared_ptr<Cantera::Solution> inletSol, std::shared_ptr<Cantera::Solution> reactorSol, double mdot, double volume)
    // ------------------------ INITIALIZE GLOBAL VARS ------------------------
    : inletThermo(inletSol->thermo()),
      reactorThermo(reactorSol->thermo()),
      reactorKinetics(reactorSol->kinetics()),
      nSpecies(inletSol->thermo()->nSpecies()),
      Y_in(inletSol->thermo()->massFractions()),               // mass fractions at the inlet
      molecularWeight(inletSol->thermo()->molecularWeights()), // molecular weight of each species (these never change)
      volume(volume),
      mdot(mdot)
    {
    }

    /// Providing implementations for the solver's virtual functions:

    // Specify guesses for the initial values.
    // Note: called during Sim1D initialization
    doublereal ctNLS_initialValue(size_t i, size_t j)
    {
        return reactorThermo->massFraction(i); // mass fractions of the provided intial reactorSol (aka the initial guess)
    }

    // Number of equations (state variables) for this reactor
    size_t ctNLS_nEqs()
    {
        return nSpecies;
    }

    // Advance the PSR to steady state
    void solveSteady()
    {
        Cantera_NonLinSol::solve();
    }

    // Specify the residual function for the system
    //  sol - iteration solution vector (input)
    //  rsd - residual vector (output)
    void ctNLS_residFunction(double *sol, double *rsd)
    {
        // ----------------------- UPDATE REACTOR STATE -----------------------
        reactorThermo->setMassFractions_NoNorm(sol);
        reactorThermo->setState_HP(inletThermo->enthalpy_mass(), PRESSURE); // keep total enthalpy constant, allow Cantera to control reactor temperature

        // ----------------------- GET REQ'D PROPERTIES -----------------------
        doublereal wdot[nSpecies];
        reactorKinetics->getNetProductionRates(wdot);

        // ----------------------- SPECIES CONSERVATION -----------------------
        for (int k = 0; k < nSpecies; k++)
        {
            rsd[k] = wdot[k] * molecularWeight[k] * volume + mdot * (Y_in[k] - sol[k]); //PSR residual equation
        }
    }

private:
    std::shared_ptr<Cantera::ThermoPhase> inletThermo;
    std::shared_ptr<Cantera::ThermoPhase> reactorThermo;
    std::shared_ptr<Cantera::Kinetics> reactorKinetics;
    int nSpecies;
    const double *Y_in;
    Cantera::vector_fp molecularWeight;
    double volume;
    double mdot;
};

int main()
{
    try
    {
        std::shared_ptr<Cantera::Solution> mixInlet(Cantera::newSolution(MECHANISM));
        std::shared_ptr<Cantera::Solution> mixOutlet(Cantera::newSolution(MECHANISM));
        mixInlet->thermo()->setState_TPX(300, PRESSURE, "O2:1 CH4:0.5 N2:3.76");
        mixOutlet->thermo()->setState_TPX(300, PRESSURE, "O2:1 CH4:0.5 N2:3.76");
        mixOutlet->thermo()->equilibrate("HP");

        double mdot;
        double volume;
        std::cout << "Enter mdot: ";
        std::cin >> mdot;
        std::cout << "Enter volume: ";
        std::cin >> volume;

        PSR reactor(mixInlet, mixOutlet, mdot, volume);
        reactor.solveSteady();

        std::cout << mixOutlet->thermo()->report() << "\n";

        return 0;
    }
    catch (Cantera::CanteraError &err)
    {
        std::cerr << err.what() << std::endl;
        return -1;
    }
}
