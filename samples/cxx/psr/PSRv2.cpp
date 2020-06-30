/// @file PSRv2.cpp
/// Perfectly Stirred Reactor Solver (v0.2)

#include "BoundaryValueProblem.h" // PSR uses the simplified BoundaryValueProblem interface to Cantera's solver (for now)

const char *MECHANISM = "gri30.yaml";
const double PRESSURE = Cantera::OneAtm;

class PSR : public BVP::BoundaryValueProblem
{
public:
    // constructor
    PSR(std::shared_ptr<Cantera::Solution> inletSol, std::shared_ptr<Cantera::Solution> reactorSol, double mdot, double volume, bool doEnergy = false, double fixedTemp = 0)
    // ------------------------ INITIALIZE GLOBAL VARS ------------------------
    : BVP::BoundaryValueProblem(inletSol->thermo()->nSpecies() + doEnergy, 1, 0.0, 1), // initialize BVP interface. 1st arg is number of solution components, 2nd arg is number of grid points (1 point for 0D...), 3rd and 4th are arbitrary (not used in this problem)
      inletThermo(inletSol->thermo()),
      reactorThermo(reactorSol->thermo()),
      reactorKinetics(reactorSol->kinetics()),
      nSpecies(inletSol->thermo()->nSpecies()),
      index_T(inletSol->thermo()->nSpecies()),                 // location of the "Temperature" component in the solution vector
      Y_in(inletSol->thermo()->massFractions()),               // mass fractions at the inlet
      molecularWeight(inletSol->thermo()->molecularWeights()), // molecular weight of each species (these never change)
      volume(volume),
      mdot(mdot),
      fixedTemp(fixedTemp * !doEnergy), // bool: true if fixed reactor temperature provided, *always* false if doEnergy = true
      doEnergy(doEnergy)

    {
        // ------------------ INITIALIZE SOLUTION COMPONENTS ------------------
        // solution vector contains nSpecies reactor mass fractions, ordered by species index.
        // if energy is enabled, reactor temp is added as the final solution component.
        // i.e.:
        //     sol = | Y_0 | Y_1 | Y_2 | ... | Y_(nSpecies-1) | Temp |
        BVP::Component initializer;
        initializer.lower = -1.0e-3;
        initializer.upper = 1.01;
        initializer.rtol = 1.0e-4;
        initializer.atol = 1.0e-9;
        initializer.refine = false;

        for (int k = 0; k < nSpecies; k++)
        {
            initializer.name = inletThermo->speciesName(k);
            setComponent(k, initializer);
        }

        if (doEnergy)
        {
            initializer.lower = 200;
            initializer.upper = 6000;
            initializer.rtol = 1;
            initializer.atol = 1;
            initializer.name = "Temperature";
            setComponent(index_T, initializer);
        }
    }

    // destructor
    virtual ~PSR() {}

    // Specify guesses for the initial values.
    virtual doublereal initialValue(size_t i, size_t j)
    {
        if (i < nSpecies)
            return reactorThermo->massFraction(i);
        return reactorThermo->temperature();
    }

    // Specify the residual function. The solver will attempt to find a solution
    // sol so that rsd is zero.
    void eval(size_t jg, double *sol, double *rsd, int *diag, double rdt)
    {
        size_t j = 0; // current grid point (only one point in a 0D problem...)

        // ----------------------- UPDATE REACTOR STATE -----------------------
        reactorThermo->setMassFractions_NoNorm(sol);
        // reactorThermo->setMassFractions(sol); -> alternative method, will compare performance between the two

        if (doEnergy)
            reactorThermo->setState_TP(sol[index_T], PRESSURE); // use the reactor temperature based on the energy equation
        else if (fixedTemp)
            reactorThermo->setState_TP(fixedTemp, PRESSURE); // use a fixed reactor temperature (if provided)
        else
            reactorThermo->setState_HP(inletThermo->enthalpy_mass(), PRESSURE); // keep total enthalpy constant, allow Cantera to control reactor temperature

        // ----------------------- GET REQ'D PROPERTIES -----------------------
        double reactorDensity = reactorThermo->density();
        double resTime = volume * reactorDensity / mdot;

        doublereal wdot[nSpecies];
        reactorKinetics->getNetProductionRates(wdot);

        // ----------------------- SPECIES CONSERVATION -----------------------
        for (int k = 0; k < nSpecies; k++)
        {
            //rsd[index(i, j)] = wdot[k] * molecularWeight[k] * volume + mdot * (Y_in[k] - sol[k])
            rsd[k] = wdot[k] * (molecularWeight[k] / reactorDensity) - (sol[k] - Y_in[k]) / resTime;

            // -------------------- TIME INTEGRATION --------------------------
            if (rdt != 0) // rdt is the reciprocal of the time step "dt".
            {
                rsd[k] -= rdt * (sol[k] - prevSoln(k, j)); // backward euler method (will be extracted from the residual)
                diag[k] = 1;                               // enable time stepping for this solution component (automatically resets each iteration)
            }
        }

        // -------------------------- ENERGY EQUATION -------------------------
        if (doEnergy)
        {
            // ------------------ GET REQ'D PROPERTIES ------------------------
            doublereal reactorCp = reactorThermo->cp_mass();

            // calculate the specific enthalpy of each species:
            doublereal h_in[nSpecies]; // inlet specific enthalpies
            inletThermo->getEnthalpy_RT(h_in);
            doublereal h[nSpecies]; // reactor specific enthalpies
            reactorThermo->getEnthalpy_RT(h);
            for (int k = 0; k < nSpecies; k++)
            {
                h_in[k] *= inletThermo->RT() * inletThermo->moleFraction(k) / inletThermo->meanMolecularWeight();
                h[k] *= reactorThermo->RT() * reactorThermo->moleFraction(k) / reactorThermo->meanMolecularWeight();
            }

            // ----------------------- EVAL ENERGY ----------------------------
            double sum1 = 0;
            double sum2 = 0;
            for (int k = 0; k < nSpecies; k++)
            {
                sum1 += Y_in[k] * (h_in[k] - h[k]);
                sum2 += h[k] * molecularWeight[k] * wdot[k];
            }
            rsd[index_T] = sum1 / (reactorCp * resTime) - sum2 / (reactorDensity * reactorCp);

            // -------------------- TIME INTEGRATION --------------------------
            if (rdt != 0)
            {
                rsd[index_T] -= rdt * (sol[index_T] - prevSoln(index_T, j));
                diag[index_T] = 1;
            }
        }
    }

private:
    std::shared_ptr<Cantera::ThermoPhase> inletThermo;
    std::shared_ptr<Cantera::ThermoPhase> reactorThermo;
    std::shared_ptr<Cantera::Kinetics> reactorKinetics;
    const int nSpecies;
    const int index_T;
    const double *Y_in;
    Cantera::vector_fp molecularWeight;
    bool doEnergy;
    double fixedTemp;
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
