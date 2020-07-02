/// @file PSRv2-1.cpp
/// Perfectly Stirred Reactor Solver (v0.2.1)
/// Changes since v0.2:
///     - Remove dependence on BVP interface
///     - Remove unused energy and fixed-temp solvers for better readability

#include "cantera/onedim.h" // general header for all 1D reacting flow problems (1D solver is used for this problem)

const char *MECHANISM = "gri30.yaml";
const double PRESSURE = Cantera::OneAtm;

class PSR : public Cantera::Domain1D
{
public:
    // constructor
    PSR(std::shared_ptr<Cantera::Solution> inletSol, std::shared_ptr<Cantera::Solution> reactorSol, double mdot, double volume)
    // ------------------------ INITIALIZE GLOBAL VARS ------------------------
    : Cantera::Domain1D(inletSol->thermo()->nSpecies()),
      inletThermo(inletSol->thermo()),
      reactorThermo(reactorSol->thermo()),
      reactorKinetics(reactorSol->kinetics()),
      nSpecies(inletSol->thermo()->nSpecies()),
      Y_in(inletSol->thermo()->massFractions()),               // mass fractions at the inlet
      molecularWeight(inletSol->thermo()->molecularWeights()), // molecular weight of each species (these never change)
      volume(volume),
      mdot(mdot)
    {
        // ------------------ INITIALIZE DOMAINS + SIMULATOR ------------------
        m_left = new Cantera::Empty1D;  // dummy terminator
        m_right = new Cantera::Empty1D; // dummy terminator

        std::vector<Cantera::Domain1D *> domains{m_left, this, m_right};
        m_sim = new Cantera::Sim1D(domains); // Sim1D instance that will control the solution process

        // ------------------ INITIALIZE SOLUTION COMPONENTS ------------------
        // solution vector contains nSpecies reactor mass fractions, ordered by species index.
        // i.e.:
        //     sol = | Y_0 | Y_1 | Y_2 | ... | Y_(nSpecies-1) |

        // bounds and error tolerances for solution mass fractions
        double lower = -1.0e-3;
        double upper = 1.01;
        double rtol = 1.0e-4;
        double atol = 1.0e-9;

        for (int k = 0; k < nSpecies; k++)
        {
            setBounds(k, lower, upper);
            setSteadyTolerances(rtol, atol, k);
            setTransientTolerances(rtol, atol, k);
            setComponentName(k, inletThermo->speciesName(k));
        }
    }

    /// Destructor. Deletes the dummy terminator domains, and the solver.
    virtual ~PSR()
    {
        delete m_left;
        delete m_right;
        delete m_sim;
    }

    // Specify guesses for the initial values.
    virtual doublereal initialValue(size_t i, size_t j)
    {
        return reactorThermo->massFraction(i); // mass fractions of the provided intial reactorSol (aka the initial guess)
    }

    /**
     * Solve the PSR problem.
     * @param loglevel controls amount of diagnostic output.
     */
    void solve(int loglevel = 0)
    {
        m_sim->solve(loglevel); // call sim1D's solver
    }

    // Specify the residual function. The solver will attempt to find a solution
    // sol so that rsd is zero.
    void eval(size_t jg, double *sol, double *rsd, int *diag, double rdt)
    {
        size_t j = 0; // current grid point (only one point in a 0D problem...)

        // ----------------------- UPDATE REACTOR STATE -----------------------
        reactorThermo->setMassFractions_NoNorm(sol);
        // reactorThermo->setMassFractions(sol); -> alternative method, will compare performance between the two

        reactorThermo->setState_HP(inletThermo->enthalpy_mass(), PRESSURE); // keep total enthalpy constant, allow Cantera to control reactor temperature

        // ----------------------- GET REQ'D PROPERTIES -----------------------
        double reactorDensity = reactorThermo->density();
        double resTime = volume * reactorDensity / mdot;

        doublereal wdot[nSpecies];
        reactorKinetics->getNetProductionRates(wdot);

        // ----------------------- SPECIES CONSERVATION -----------------------
        for (int k = 0; k < nSpecies; k++)
        {
            rsd[k] = wdot[k] * (molecularWeight[k] / reactorDensity) - (sol[k] - Y_in[k]) / resTime;
            //rsd[k] = wdot[k] * molecularWeight[k] * volume + mdot * (Y_in[k] - sol[k]) -> alternative method (should be equivalent)

            // -------------------- TIME INTEGRATION --------------------------
            if (rdt != 0) // rdt is the reciprocal of the time step "dt".
            {
                rsd[k] -= rdt * (sol[k] - prevSoln(k, j)); // backward euler method (will be extracted from the residual)
                diag[k] = 1;                               // enable time stepping for this solution component (automatically resets each iteration)
            }
        }
    }

private:
    Cantera::Domain1D *m_left;
    Cantera::Domain1D *m_right;
    Cantera::Sim1D *m_sim;
    std::shared_ptr<Cantera::ThermoPhase> inletThermo;
    std::shared_ptr<Cantera::ThermoPhase> reactorThermo;
    std::shared_ptr<Cantera::Kinetics> reactorKinetics;
    const int nSpecies;
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
