/* @file Cantera_NonLinSol.h
 * A nonlinear algebraic system solver, built upon Cantera's 1D multi-domain damped newton solver.
 *
 * Cantera_NonLinSol is a simplified interface to the Cantera solver, designed for solving standard
 * systems of nonlinear algebraic equations. Systems are solved by Cantera as one-domain, one-point
 * problems. */

#include "cantera/onedim.h"
#include <fstream>

class Cantera_NonLinSol : public Cantera::Domain1D
{
public:
    Cantera_NonLinSol() : outputFile("YOLO.csv") {
        // outputFile << "time, pressure, enthalpy, temperature, mass, volume, internal energy";
        // for (int k = 0; k < gas->nSpecies(); k++)
        //     outputFile << "Y_" << gas->speciesName(k) << ", ";
        // outputFile << std::endl;
    }
/* TO USE THIS SOLVER:
 *     1. Include this header file in your code.
 *          #include "path/to/Cantera_NonLinSol.h"
 *     2. Create a Cantera_NonLinSol child class, where you'll set up your problem.
 *          class YourClass : public Cantera_NonLinSol
 *     3. In the child class, provide implementations for problem-specific functions.
 *          void residFunction(args) { ... }
 *          doublereal initialValue(size_t i, size_t j) { ... }
 *          size_t neq() { ... }
 *     4. Initialize your class and call its inherited solve() function.
 *          YourClassObject.solve();
 *     5. (Optional) Reconfigure solver settings. 
 *          YourClassObject.reconfigure(neq, lowerBound, upperBound, rtol, atol);
 */

/// IMPLEMENT THESE FUNCTIONS:

    // Specify the residual function for the system
    //  sol - iteration solution vector (input)
    //  rsd - residual vector (output)
    virtual void ctNLS_residFunction(double *sol, double *rsd) = 0;

    // Specify guesses for the initial values.
    //  Note: called during Sim1D initialization
    virtual doublereal ctNLS_initialValue(size_t i) = 0;

    // Number of equations (state variables) for this reactor
    virtual size_t ctNLS_nEqs() = 0;

/// CALLABLE FUNCTIONS:

    void reconfigure(int neq, double lowerBound = -1.0e-3, double upperBound = 1.01,
                     double rtol = 1.0e-4, double atol = 1.0e-9)
    {
        Domain1D::resize(neq, 1);
        for (int i = 0; i < neq; i++)
        {
            Domain1D::setBounds(i, lowerBound, upperBound);
            Domain1D::setSteadyTolerances(rtol, atol, i);
            Domain1D::setTransientTolerances(rtol, atol, i);
            Domain1D::setComponentName(i, std::to_string(i));
        }
        Domain1D::setBounds(0, 0, 100);
        Domain1D::setBounds(1, 0, 5);
        Domain1D::setBounds(2, -10000000, 10000000);
    }

    /**
     * Solve the nonlinear algebraic system.
     * @param loglevel controls amount of diagnostic output.
     */
    void solve(int loglevel = 0)
    {
        if (Domain1D::nComponents() != ctNLS_nEqs())
            reconfigure(ctNLS_nEqs());

        std::vector<Cantera::Domain1D *> domains{this};
        Cantera::Sim1D(domains).solve(loglevel);
    }

    std::ofstream outputFile;
    int counter = 1;
private:
/// INTERNAL FUNCTIONS:
    

    // Implementing the residual function for the Cantera solver to use. Handles time integration, calls subclass residFunction for residuals.
    void eval(size_t jg, double *sol, double *rsd, int *timeintMask, double rdt)
    {
        // std::cout << "\n\n----- sol components -----\n";
        // for (int i = 0; i < m_nv; i++)
        //     std::cout << sol[i] << " ";

        ctNLS_residFunction(sol, rsd); // call subclass residFunction() implementation to update rsd
        counter++;
        // std::cout << "\n\n----- iter" << counter << " sol components -----\n";
        // for (int i = 0; i < m_nv; i++)
        //     std::cout << sol[i] << " ";
        for (int i = 0; i < m_nv; i++)
            outputFile << sol[i] << ", ";
        // std::cout << "\n\n----- iter" << counter << " rsd components -----\n";
        // for (int i = 0; i < m_nv; i++)
        //     std::cout << rsd[i] << " ";
        // for (int i = 0; i < m_nv; i++)
        //     outputFile << rsd[i] << ", ";
        outputFile << std::endl;
        

        if (rdt == 0) {
            //outputFile << "no" << std::endl;
            //std::cout << "no" << std::endl;
            return; // rdt is the reciprocal of the time step "dt"; rdt != 0 for time integration...
        }
        
        //outputFile << "yes" << std::endl;
        //std::cout << "yes" << std::endl;
        // -------------------- TIME INTEGRATION --------------------------
        for (int i = 0; i < ctNLS_nEqs(); i++)
        {
            std::cout << "time step " << rdt * (sol[i] - prevSoln(i, 0)) << "\n";
            if(i == 0)
                std::cout << "\ntimestep:\nsol[0]: " << sol[0] << " prevsoln[0]: " << prevSoln(0,0) << " rsd " << rsd[i] << " rsd -= " << rdt * (sol[i] - prevSoln(i, 0));
            rsd[i] -= rdt * (sol[i] - prevSoln(i, 0)); // backward euler method (result will be appropriately extracted from the residual)
            if(i==0)
                std::cout << " rsd = " << rsd[i] << "\n";
            timeintMask[i] = 1;                        // enable time stepping for this solution component (automatically resets each iteration)
        }
    }

    // Implementing the initial value function for the Cantera solver. Grid point j is always 0 (the only point in the single-point simulation), and thus unneeded
    doublereal initialValue(size_t n, size_t j)
    {
        return ctNLS_initialValue(n);
    }
};
