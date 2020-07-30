/// @file NonLinSol_user.h
/// An interface to enable the use of nonlinear algebraic solver applications.
/// -> Derive classes from `NonLinSol_user` to provide residual and initial guess functions to a nonlinsolver

#ifndef CT_NLS_H
#define CT_NLS_H

#include "cantera/base/ct_defs.h"
#include "cantera/base/ctexceptions.h"
#include "cantera/base/global.h"

class NonLinSol_user
{
public:
    NonLinSol_user(){};
    ~NonLinSol_user(){};

    // Specify the residual function for the system
    //  sol - iteration solution vector (input)
    //  rsd - residual vector (output)
    virtual void residFunction(double *sol, double *rsd) = 0;

    // Specify guesses for the initial values.
    // Note: called during Sim1D initialization
    virtual doublereal initialGuess(size_t i, size_t j) = 0;

    // Number of equations (state variables) for this reactor
    virtual size_t neq() = 0;
};
#endif
