// Copyright (c) 2009,2010,2011,2012 Brian C. Fabien
// All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the license.txt file.

#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <float.h>

#include "runtime.h"
#include "vector.h"
#include "matrix.h"
#include "cholesky.h"
#include "qp.h"
#include "scqp.h"
#include "equations.h"
#include "nlp.h"
#include "l1sqp.h"
#include "linfsqp.h"
#include "nlpsolve.h"

// Solve a NLP problem
//
// Function input:
//
// double (*fhg) (Vector x, Vector h, Vector g):
//      input: x: n-Vector  (unknown variables)
//      output: h: m-Vector (equality constraints)
//      output: g: p-Vector (inequality constraints)
// void (*Dfhg)(Vector x, Vector df, Matrix dh, Matrix dg): (optional, can be set to NULL):
//      input: x: n-Vector  (unknown variables)
//      output: Df: n-Vector (function gradient)
//      output: Dh: m-by-n Matrix (Jacobian equality constraints)
//      output: Dg: p-by-n Matrix (Jacobian inequality constraints)
// x: n-Vector
//      input: solution estimate, output: solution to the NLP
// lambda: m-Vector: Lagrange multipliers associated with the equality constraints
// mu: p-Vector: Lagrange multipliers associated with the inequality constraints
// options:
//      option.method = L1SQPmethod | LinfSQPmethod | Hybridmethod
//      option.tolerance = 1.0e-6;
//      option.maximumIterations = 1000;
//      option.displayData = NO;
//
// Function output:
//
//  stats: NLPStatistics data structure
//      stats.f: function value
//      stats.normc: norm of the active constraints
//      stats.normT: norm of the KKT conditions
//      stats.numberOfIterations;
//      stats.numberOfFunctionCalls;
//      stats.numberOfGradientCalls;
//      stats.numberOfQPSolves;
//      stats.numberOfRefinements;
//      stats.numberOfSOC;
//      stats.exitCode;
NLPStatistics
NLPsolve(double (*fhg) (Vector x, Vector h, Vector g),
         void (*Dfhg) (Vector x, Vector df, Matrix dh, Matrix dg),
         Vector x, Vector lambda, Vector mu,
         NLPOptions options) {

    NLPStatistics stats;
    int n = 0, m = 0, p = 0;
    stats.exitCode              = 0;
    stats.f                     = 0;
    stats.normc                 = 0;
    stats.normT                 = 0;
    stats.numberOfIterations    = 0;
    stats.numberOfFunctionCalls = 0;
    stats.numberOfGradientCalls = 0;
    stats.numberOfQPSolves      = 0;
    stats.numberOfRefinements   = 0;
    stats.numberOfSOC           = 0;
    
    if (x == NULL) {
        RuntimeWarning("NLPsolve: no initial solution estimate given");
        return stats;
    } else {
        n = VEC_LENGTH(x);
    }
    if (lambda != NULL) {
        m = VEC_LENGTH(lambda);
    }
    if (mu != NULL) {
        p = VEC_LENGTH(mu);
    }

    NLP nlp = NLPNew();
    nlp->n = n;
    nlp->m = m;
    nlp->p = p;
    nlp->fhg = fhg;
    nlp->Dfhg = Dfhg;

    stats.n = n;
    stats.m = m;
    stats.p = p;
    
    if (options.method == L1SQPmethod) {
        L1SQP sqp0 = L1SQPNew(nlp);
        sqp0->tolerance = options.tolerance;
        sqp0->maximumIterations = options.maximumIterations;
        sqp0->displayData = options.displayData;
        stats.exitCode = L1SQPSolve(sqp0, x, lambda, mu);
        stats.f                     = sqp0->f;
        stats.normc                 = sqp0->normc;
        stats.normT                 = sqp0->normT;
        stats.numberOfIterations    = sqp0->numberOfIterations;
        stats.numberOfFunctionCalls = sqp0->numberOfFunctionCalls;
        stats.numberOfGradientCalls = sqp0->numberOfGradientCalls;
        stats.numberOfQPSolves      = sqp0->numberOfQPSolves;
        stats.numberOfRefinements   = sqp0->numberOfRefinements;
        stats.numberOfSOC           = sqp0->numberOfSOC;
        sprintf(stats.methodName, "%s", "L1SQPmethod");
        L1SQPDelete(sqp0);
    } else if (options.method == LinfSQPmethod) {
        LinfSQP sqp1 = LinfSQPNew(nlp);
        sqp1->tolerance = options.tolerance;
        sqp1->maximumIterations = options.maximumIterations;
        sqp1->displayData = options.displayData;
        stats.exitCode = LinfSQPSolve(sqp1, x, lambda, mu);
        stats.f                     = sqp1->f;
        stats.normc                 = sqp1->normc;
        stats.normT                 = sqp1->normT;
        stats.numberOfIterations    = sqp1->numberOfIterations;
        stats.numberOfFunctionCalls = sqp1->numberOfFunctionCalls;
        stats.numberOfGradientCalls = sqp1->numberOfGradientCalls;
        stats.numberOfQPSolves      = sqp1->numberOfQPSolves;
        stats.numberOfRefinements   = sqp1->numberOfRefinements;
        stats.numberOfSOC           = sqp1->numberOfSOC;
        sprintf(stats.methodName, "%s", "LinfSQPmethod");
        LinfSQPDelete(sqp1);
    } /* else if (options.method == HybridSQPmethod) {
        L1SQP sqp2 = L1SQPNew(nlp);
        sqp2->tolerance = options.tolerance;
        sqp2->maximumIterations = options.maximumIterations;
        sqp2->displayData = options.displayData;
        stats.exitCode = L1SQPSolveX(sqp2, x, lambda, mu);
        stats.f                     = sqp2->f;
        stats.normc                 = sqp2->normc;
        stats.normT                 = sqp2->normT;
        stats.numberOfIterations    = sqp2->numberOfIterations;
        stats.numberOfFunctionCalls = sqp2->numberOfFunctionCalls;
        stats.numberOfGradientCalls = sqp2->numberOfGradientCalls;
        stats.numberOfQPSolves      = sqp2->numberOfQPSolves;
        stats.numberOfRefinements   = sqp2->numberOfRefinements;
        stats.numberOfSOC           = sqp2->numberOfSOC; 
        L1SQPDelete(sqp2);
    } else if (options.method == LagNewtmethod) {
        LagNewt sqp3 = LagNewtNew(nlp);
        sqp3->tolerance = options.tolerance;
        sqp3->maximumIterations = options.maximumIterations;
        sqp3->displayData = options.displayData;
        stats.exitCode = LagNewtSolve(sqp3, x, lambda, mu);
        stats.f                     = sqp3->f;
        stats.normc                 = sqp3->normc;
        stats.normT                 = sqp3->normT;
        stats.numberOfIterations    = sqp3->numberOfIterations;
        stats.numberOfFunctionCalls = sqp3->numberOfFunctionCalls;
        stats.numberOfGradientCalls = sqp3->numberOfGradientCalls;
        stats.numberOfQPSolves      = sqp3->numberOfQPSolves;
        stats.numberOfRefinements   = sqp3->numberOfRefinements;
        //stats.numberOfSOC           = sqp3->numberOfSOC; 
        LagNewtDelete(sqp3);
    }*/ else {
        RuntimeWarning("NLPsolve: invalid solution method specified");
    }
    
    NLPDelete(nlp);
    
    return stats;
}

void NLPPrintStatistics(NLPStatistics s) {
    printf("exitCode: %d\n",    s.exitCode);
    printf("f: %e\n",           s.f);
    printf("|c|: %e\n",         s.normc);
    printf("|T|: %e\n",         s.normT);
    printf("# iter: %d\n",      s.numberOfIterations);
    printf("# fun: %d\n",       s.numberOfFunctionCalls);
    printf("# grad: %d\n",      s.numberOfGradientCalls);
    printf("# QP: %d\n",        s.numberOfQPSolves);
    printf("# refine: %d\n",    s.numberOfRefinements);
    printf("# SOC: %d\n",       s.numberOfSOC);
}
