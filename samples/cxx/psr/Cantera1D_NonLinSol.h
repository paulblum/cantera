/// @file Cantera1D_NonLinSol.h
/// A nonlinear algebraic system solver, built on top of Cantera's 1D solver.

#include "cantera/onedim.h"
#include "NonLinSol_user.h"

class Cantera1D_NonLinSol : public Cantera::Domain1D
{
public:
    /// Constructor
    Cantera1D_NonLinSol(NonLinSol_user *user)
        : Cantera::Domain1D(user->neq()), //initialize base class
          m_user(user),                   //pointer to user class, to call user defined functions
          m_neq(user->neq())              //number of equations
    {
        // ------------------ INITIALIZE DOMAINS + SIMULATOR ------------------
        m_left = new Cantera::Empty1D;                                   // dummy terminator
        m_right = new Cantera::Empty1D;                                  // dummy terminator
        std::vector<Cantera::Domain1D *> domains{m_left, this, m_right}; // vector of pointers to domains to be linked (reqd to init sim1D)
        m_sim = new Cantera::Sim1D(domains);                             // Sim1D instance that will control the solution process

        // ------------------ INITIALIZE SOLUTION COMPONENTS ------------------
        // bounds and error tolerances for all solution components
        double lower = -1.0e-3;
        double upper = 1.01;
        double rtol = 1.0e-4;
        double atol = 1.0e-9;

        for (int i = 0; i < m_neq; i++)
        {
            setBounds(i, lower, upper); //methods inherited from Domain1D
            setSteadyTolerances(rtol, atol, i);
            setTransientTolerances(rtol, atol, i);
            setComponentName(i, std::to_string(i));
        }
    }

    /// Destructor. Deletes the dummy terminator domains, and the solver.
    ~Cantera1D_NonLinSol()
    {
        delete m_left;
        delete m_right;
        delete m_sim;
    }

    doublereal initialValue(size_t i, size_t j)
    {
        return m_user->initialGuess(i, j); //use pointer to user class to call user-supplied initial function
    }

    /**
     * Solve the nonlinear algebraic system.
     * @param loglevel controls amount of diagnostic output.
     */
    void solve(int loglevel = 0)
    {
        m_sim->solve(loglevel); // call sim1D's solver
    }

    // Specify the residual function. The solver will attempt to find a solution
    // sol so that rsd is zero.
    void eval(size_t jg, double *sol, double *rsd, int *timeintMask, double rdt)
    {
        m_user->residFunction(sol, rsd); // call subclass residFunction() implementation to update rsd

        if (rdt == 0)
            return; // rdt is the reciprocal of the time step "dt"; rdt != 0 for time integration...

        // -------------------- TIME INTEGRATION --------------------------
        for (int i = 0; i < m_neq; i++)
        {
            rsd[i] -= rdt * (sol[i] - prevSoln(i, 0)); // backward euler method (result will be appropriately extracted from the residual)
            timeintMask[i] = 1;                        // enable time stepping for this solution component (automatically resets each iteration)
        }
    }

private:
    NonLinSol_user *m_user;
    Cantera::Domain1D *m_left;
    Cantera::Domain1D *m_right;
    Cantera::Sim1D *m_sim;
    int m_neq;
};
