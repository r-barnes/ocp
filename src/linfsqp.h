// Copyright (c) 2009,2010,2011,2012 Brian C. Fabien
// All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the license.txt file.

/* 
 * File:   LinfSQP.h
 * Author: fabien
 *
 * Created on July 10, 2010, 6:48 PM
 */

#ifndef _LinfSQP_H
#define	_LinfSQP_H

#ifdef	__cplusplus
extern "C" {
#endif

struct linfsqp_ {
    // Input
    double tolerance;       // Convergence tolerance (1.0e-6 default) 
    int maximumIterations;  // Maximum number of itersations (500 default) 
    int displayData;        // Display progress of the solver if YES 

    // Output
    Vector x;       // The optimal parameters: double[n] 
    Vector lambda;  // Lagrange multipliers associated with the equality constraints: double[m] 
    Vector mu;      // Lagrange multipliers associated with the inequality constraints: double[p] 
    double f;       // The optimal function value 
    double normT;   // Norm of the gradient of the Lagrangian 
    double normc;   // Norm of the active constraints 
    int numberOfIterations;     // Number of iterations performed 
    int numberOfFunctionCalls;  // Number of function calls 
    int numberOfGradientCalls;  // Number of gradient evaluations performed 
    int numberOfRefinements;    // Number of iterative refinements in the QP method 
    int numberOfQPSolves;       // Number of QP problems solved
    int numberOfSOC;            // Number of second-order corrections
    
    // Private
    NLP nlp;
    QP qp, rqp;
    SCQP scqp, rscqp;
    int n, m, p;
    double f0, theta0, Dtheta, alpha;
    Vector rqpx, rqplambda, rqpmu;
    Vector d, df, h, g, df0, h0, g0, kkt, Hd, d_soc;
    Matrix H, dh, dg, dh0, dg0;
    Vector x0;
    Vector _v, _q, _Bp, _r;
    double tol, nh, ng, nx, nd, normkkt, rho_hat;
    int nfun, ngrad, niter, nsoc, maxIter, nrelax;
    int code;
    int relaxed, do_soc;
    int nrefine;
    int done, display, fd_gradient;
};

typedef struct linfsqp_ * LinfSQP;

LinfSQP LinfSQPNew(NLP nlp);
int     LinfSQPSolve(LinfSQP a, Vector x0, Vector lambda0, Vector mu0);
void    LinfSQPDelete(LinfSQP a);

#ifdef	__cplusplus
}
#endif

#endif	/* _LinfSQP_H */
