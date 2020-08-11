/* @file Cantera_NonLinSol.h
 * A nonlinear algebraic system solver, built upon Cantera's 1D multi-domain damped newton solver.
 *
 * Cantera_NonLinSol is a simplified interface to the Cantera solver, designed for solving standard
 * systems of nonlinear algebraic equations. Systems are solved by Cantera as one-domain, one-point
 * problems. */

#include "cantera/onedim.h"

class Cantera_NonLinSol : public Cantera::Domain1D
{
public:
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
    virtual void residFunction(double *sol, double *rsd) = 0;

    // Specify guesses for the initial values.
    //  Note: called during Sim1D initialization
    virtual doublereal initialValue(size_t i, size_t j) = 0;

    // Number of equations (state variables) for this reactor
    virtual size_t neq() = 0;

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
    }

    /**
     * Solve the nonlinear algebraic system.
     * @param loglevel controls amount of diagnostic output.
     */
    void solve(int loglevel = 0)
    {
        if (Domain1D::nComponents() != neq())
            reconfigure(neq());

        std::vector<Cantera::Domain1D *> domains{this};
        Cantera::Sim1D(domains).solve(loglevel);
    }

private:
/// INTERNAL FUNCTION:

    // Implementing the residual function for the Cantera solver to use. Handles time integration, calls subclass residFunction for residuals.
    void eval(size_t jg, double *sol, double *rsd, int *timeintMask, double rdt)
    {
        residFunction(sol, rsd); // call subclass residFunction() implementation to update rsd

        if (rdt == 0)
            return; // rdt is the reciprocal of the time step "dt"; rdt != 0 for time integration...

        // -------------------- TIME INTEGRATION --------------------------
        for (int i = 0; i < neq(); i++)
        {
            rsd[i] -= rdt * (sol[i] - prevSoln(i, 0)); // backward euler method (result will be appropriately extracted from the residual)
            timeintMask[i] = 1;                        // enable time stepping for this solution component (automatically resets each iteration)
        }
    }
};
