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

// line search parameters

#define SIGMA 0.1
#define BETA (1.0/3.0)
#define FAC 1.0e-15

// maximum line search iterations

#define MAX_LS_ITER 50

// second-order corrections

#define SOC 1

int  L1SQPcheckNLP(L1SQP a, Vector x0, Vector lambda0, Vector mu0);
void L1SQPshowInfo(L1SQP a);
void L1SQPcomputeSearchDirection(L1SQP a);
void L1SQPterminate(L1SQP a);
void L1SQPlineSearch(L1SQP a);
void L1SQPupdate(L1SQP a);
void L1SQPallocateStorage(L1SQP a);
void L1SQPdeleteStorage(L1SQP a);
int  L1SQPsolveQP(L1SQP a);
void L1SQPcomputeDtheta(L1SQP a);
int  L1SQPsolveRelaxedQP(L1SQP a);
void L1SQPcomputeDthetaRelaxed(L1SQP a);
void L1SQPsocStep(L1SQP a);
void L1SQPBFGSupdate(L1SQP a);

// create a L1SQP data structure
L1SQP L1SQPNew(NLP nlp) {
    L1SQP a;
    a = (L1SQP) malloc(sizeof (struct l1sqp_));

    // Input

    a->tolerance         = 1.0e-6; // Convergence tolerance (1.0e-6 default)
    a->maximumIterations = 500;    // Maximum number of itersations (500 default)
    a->displayData       = NO;     // Display progress of the solver if true

    // Input/Output

    a->x      = NULL;   // The optimal parameters: double[n]
    a->lambda = NULL;   // Lagrange multipliers associated with the equality constraints: double[m]
    a->mu     = NULL;   // Lagrange multipliers associated with the inequality constraints: double[p]

    // Output

    a->f                     = 0.0; // The optimal function value
    a->normT                 = 0.0; // Norm of the gradient of the Lagrangian
    a->normc                 = 0.0; // Norm of the active constraints
    a->numberOfIterations    = 0;   // Number of iterations performed
    a->numberOfFunctionCalls = 0;   // Number of function calls
    a->numberOfGradientCalls = 0;   // Number of gradient evaluations performed
    a->numberOfRefinements   = 0;   // Number of iterative refinements in the QP method
    a->numberOfQPSolves      = 0;   // Number of QP problems solved
    a->numberOfSOC           = 0;   // Number of second-order corrections

    // Private

    a->nlp         = nlp;    // the NLP problem
    a->qp          = NULL;   // the QP
    a->rqp         = NULL;   // the relaxed QP
    a->scqp        = NULL;   // the SCQP solver
    a->rscqp       = NULL;   // the relaxed SCQP solver
    a->n           = 0;      // # of parameters
    a->m           = 0;      // # equlity constraints
    a->p           = 0;      // # inequlity constraints
    a->f0          = 0.0;
    a->rho_hat     = 1.0;
    a->theta0      = 0.0;    // L1 penalty function
    a->Dtheta      = 0.0;    // 1st-order change in the penalty function
    a->alpha       = 0.0;    // search step
    a->rqpx        = NULL;   // x from the relaxed QP
    a->rqplambda   = NULL;   // lambda from the relaxed QP
    a->rqpmu       = NULL;   // mu from the relaxed QP
    a->d           = NULL;   // search direction
    a->df          = NULL;   // gradient of the function
    a->h           = NULL;   // equality constraints
    a->g           = NULL;   // inequality constraints
    a->df0         = NULL;   // gradient of the function
    a->h0          = NULL;   // equality constraints
    a->g0          = NULL;   // inequality constraints
    a->kkt         = NULL;   // the KKT vector
    a->Hd          = NULL;   // H*d
    a->d_soc       = NULL;   // second-order correction
    a->H           = NULL;   // Hessian
    a->dh          = NULL;   // Jacobian of the equality constraints
    a->dg          = NULL;   // Jacobian of the inequality constraints
    a->dh0         = NULL;   // Jacobian of the equality constraints
    a->dg0         = NULL;   // Jacobian of the inequality constraints
    a->x0          = NULL;   // current solution
    a->rho         = NULL;   // penalty weights
    a->_v          = NULL;   // BFGS update
    a->_q          = NULL;   // BFGS update
    a->_Bp         = NULL;   // BFGS update
    a->_r          = NULL;   // BFGS update
    a->tol         = 0.0;    // tolerance
    a->nh          = 0.0;    // norm h
    a->ng          = 0.0;    // norm g
    a->nx          = 0.0;    // norm x
    a->nd          = 0.0;    // norm d
    a->normkkt     = 0.0;    // norm of the KKT condition
    a->nfun        = 0;      // # func evaluations
    a->ngrad       = 0;      // # grad evaluations
    a->niter       = 0;      // # of iterations
    a->nsoc        = 0;      // # of second order corrections
    a->maxIter     = 0;
    a->nrelax      = 0;      // # of relaxed QP problems
    a->code        = 0;      // return code
    a->relaxed     = NO;     // is this a relaxed QP?
    a->do_soc      = NO;     // perform second order correction?
    a->nrefine     = 0;      // # of iterrative refinemens in the GI algorithm
    a->done        = NO;
    a->display     = YES;    // show progress of the algorithm?
    a->fd_gradient = YES;    // use forward difference gradients

    a->_y = NULL;
    a->_d = NULL;
    a->_L = NULL;

    return a;
}

int L1SQPcheckNLP(L1SQP a, Vector x0, Vector lambda0, Vector mu0) {
    // check nlp, nlp->n, nlp->m, nlp->p, nlp->fhg, nlp->Dfhg
    if (a->nlp == NULL) {
        RuntimeWarning("L1SQPcheckNLP: nlp == NULL");
        return 1;
    }
    if (a->nlp->n < 0) {
        RuntimeWarning("L1SQPcheckNLP: nlp->n < 0");
        return 1;
    }
    if (a->nlp->m < 0) {
        RuntimeWarning("L1SQPcheckNLP: nlp->m < 0");
        return 1;
    }
    if (a->nlp->m > a->nlp->n) {
        RuntimeWarning("L1SQPcheckNLP: nlp->m > nlp->n");
        return 1;
    }
    if (a->nlp->p < 0) {
        RuntimeWarning("L1SQPcheckNLP: nlp->p < 0");
        return 1;
    }
    if (a->nlp->fhg == NULL) {
        RuntimeWarning("L1SQPcheckNLP: nlp->fhg == NULL");
        return 1;
    }
    if (a->nlp->Dfhg != NULL) {
        a->fd_gradient = NO;
    }
    // check x0, lambda0, mu0
    if (x0 == NULL) {
        RuntimeWarning("L1SQPcheckNLP: x0 == NULL");
        return 1;
    }
    if (x0->r != a->nlp->n) {
        RuntimeWarning("L1SQPcheckNLP: invalid dimension for x0");
        return 1;
    }
    if ((a->nlp->m > 0) && (lambda0 == NULL)) {
        RuntimeWarning("L1SQPcheckNLP: lambda0 == NULL");
        return 1;
    }
    if ((a->nlp->m > 0) && (lambda0->r != a->nlp->m)) {
        RuntimeWarning("L1SQPcheckNLP: invalid dimension for lambda0");
        return 1;
    }
    if ((a->nlp->p > 0) && (mu0 == NULL)) {
        RuntimeWarning("L1SQPcheckNLP: mu0 == NULL");
        return 1;
    }
    if ((a->nlp->p > 0) && (mu0->r != a->nlp->p)) {
        RuntimeWarning("L1SQPcheckNLP: invalid dimension for mu0");
        return 1;
    }

    a->x = x0;
    a->lambda = lambda0;
    a->mu = mu0;

    if (a->tolerance < DBL_EPSILON) {
        a->tol = 1.0e-6;
        RuntimeWarning("tolerance < DBL_EPSILON, using 1.0e-6");
    } else {
        a->tol = a->tolerance;
    }
    if (a->maximumIterations < 1) {
        a->maxIter = 500;
        RuntimeWarning("maximumIterations < 0, using 500");
    } else {
        a->maxIter = a->maximumIterations;
    }
    if ((a->displayData != YES) && (a->displayData != NO)) {
        a->display = YES;
    } else {
        a->display = a->displayData;
    }

    L1SQPallocateStorage(a);

    return 0;
}

void L1SQPallocateStorage(L1SQP a) {
    // allocate storage
    int i;
    int n = a->nlp->n;
    int m = a->nlp->m;
    int p = a->nlp->p;

    a->n = n;
    a->m = m;
    a->p = p;

    a->x0 = VectorNew(n);           // current solution
    VectorSetEqual(a->x0, a->x);
    a->d = VectorNew(n);            // search direction
    a->d_soc = VectorNew(n);        // second-order correction
    a->kkt = VectorNew(n);          // the KKT vector

    a->df = VectorNew(n);           // gradient of the function
    a->df0 = VectorNew(n);          // gradient of the function

    a->H = MatrixNew(n, n);         // Hessian
    for (i = 0; i < n; i++) {
        a->H->e[i][i] = 1.0;
    }
    a->Hd = VectorNew(n);           // H*d

    if ((m + p) > 0) {
        a->rho = VectorNew(m + p); // penalty weights
        VectorSetAllTo(a->rho, 1.0);
    }
    if (m > 0) {
        a->h = VectorNew(m);        // equality constraints
        a->dh = MatrixNew(m, n);    // Jacobian of the equality constraints
        a->h0 = VectorNew(m);       // equality constraints
        a->dh0 = MatrixNew(m, n);   // Jacobian of the equality constraints
    }
    if (p > 0) {
        a->g = VectorNew(p);        // inequality constraints
        a->g0 = VectorNew(p);       // inequality constraints
        a->dg = MatrixNew(p, n);    // Jacobian of the inequality constraints
        a->dg0 = MatrixNew(p, n);   // Jacobian of the inequality constraints
    }
    a->qp = QPNew(); // the QP
    a->qp->H = a->H;
    a->qp->c = a->df;
    a->qp->A1 = a->dh;
    a->qp->b1 = a->h;
    a->qp->A2 = a->dg;
    a->qp->b2 = a->g;

    a->scqp = SCQPNew(a->qp); // the SCQP solver
    a->scqp->hessianType = FullHessian; //FullHessian;  InverseLowerCholeskyFactor;

    a->rqp = NULL;          // the relaxed QP
    a->rscqp = NULL;        // the relaxed SCQP solver
    a->rqpx = NULL;         // x from the relaxed QP
    a->rqplambda = NULL;    // lambda from the relaxed QP
    a->rqpmu = NULL;        // mu from the relaxed QP

    a->_v = VectorNew(n);   // BFGS update
    a->_q = VectorNew(n);   // BFGS update
    a->_Bp = VectorNew(n);  // BFGS update
    a->_r = VectorNew(n);   // BFGS update

    a->nfun = 0;
    a->ngrad = 0;
    a->nsoc = 0;
    a->nrelax = 0;
    a->nrefine = 0;

    a->f = a->nlp->fhg(a->x0, a->h, a->g);
    a->nfun++;
    if (a->fd_gradient == YES) {
        NLPDfhg(a->nlp, a->x0, a->df, a->dh, a->dg);
    } else {
        a->nlp->Dfhg(a->x0, a->df, a->dh, a->dg);
    }
    a-> ngrad++;

    VectorSetEqual(a->df0, a->df);
    a->nh = 0;
    a->ng = 0;
    if (m > 0) {
        VectorSetEqual(a->h0, a->h);
        MatrixSetEqual(a->dh0, a->dh);
        a->nh = VectorNorm(a->h);
    }
    if (p > 0) {
        VectorSetEqual(a->g0, a->g);
        MatrixSetEqual(a->dg0, a->dg);
        a->ng = fmax(0.0, VectorMax(a->g));
    }
    a->rho_hat = 1.0;
    a->theta0 = 0;
    a->Dtheta = 0;
    a->normc = fmax(a->nh, a->ng);
    a->done = NO;
    a->code = 0;

    a->_y = VectorNew(n);
    a->_d = VectorNew(n);
    a->_L = MatrixNew(n, n);
    a->_Li = MatrixNew(n, n);
    for (i = 0; i < n; i++) {
        a->_Li->e[i][i] = 1.0;
    }
}

void L1SQPDelete(L1SQP a) {
    if (a == NULL) {
        return;
    }
    free(a);
}

void L1SQPdeleteStorage(L1SQP a) {
    //int n = a->n;
    int m = a->m;
    int p = a->p;
    VectorDelete(a->x0);
    VectorDelete(a->d);
    VectorDelete(a->d_soc);
    VectorDelete(a->kkt);

    VectorDelete(a->df);
    VectorDelete(a->df0);

    MatrixDelete(a->H);

    VectorDelete(a->Hd);

    if ((m + p) > 0) {
        VectorDelete(a->rho);
    }
    if (m > 0) {
        VectorDelete(a->h);
        MatrixDelete(a->dh);
        VectorDelete(a->h0);
        MatrixDelete(a->dh0);
    }
    if (p > 0) {
        VectorDelete(a->g);
        VectorDelete(a->g0);
        MatrixDelete(a->dg);
        MatrixDelete(a->dg0);
    }
    QPDelete(a->qp);
    SCQPDelete(a->scqp);
    if (a->rqp != NULL) {
        // H, c, A1, b1, A2, b2, rqp
        if (a->rqp->H != NULL) MatrixDelete(a->rqp->H);
        if (a->rqp->A1 != NULL) MatrixDelete(a->rqp->A1);
        if (a->rqp->A2 != NULL) MatrixDelete(a->rqp->A2);
        if (a->rqp->c != NULL) VectorDelete(a->rqp->c);
        if (a->rqp->b1 != NULL) VectorDelete(a->rqp->b1);
        if (a->rqp->b2 != NULL) VectorDelete(a->rqp->b2);
        QPDelete(a->rqp);
    }
    if (a->rscqp != NULL) SCQPDelete(a->rscqp);
    if (a->rqpx != NULL) VectorDelete(a->rqpx);
    if (a->rqplambda != NULL) VectorDelete(a->rqplambda);
    if (a->rqpmu != NULL) VectorDelete(a->rqpmu);
    VectorDelete(a->_v);
    VectorDelete(a->_q);
    VectorDelete(a->_Bp);
    VectorDelete(a->_r);

    VectorDelete(a->_y);
    VectorDelete(a->_d);
    MatrixDelete(a->_L);
    MatrixDelete(a->_Li);
}

void L1SQPshowInfo(L1SQP a) {
    if (a->niter == 1) {
        printf("iter: %4d f: %8.4e |c|: %8.4e\n",
               a->niter, a->f, a->normc);
    } else {
        printf("iter: %4d f: %8.4e |c|: %8.4e |T|: %8.4e\n",
               a->niter, a->f, a->normc, a->normT);
    }
}

void L1SQPterminate(L1SQP a) {
    a->numberOfIterations = a->niter;
    a->numberOfFunctionCalls = a->nfun;
    a->numberOfGradientCalls = a->ngrad;
    a->numberOfRefinements = a->nrefine;
    a->numberOfSOC = a->nsoc;
    if (a->display) {
        printf("$ L1SQP\n");
        printf("$ n: %d\n", a->n);
        printf("$ m: %d\n", a->m);
        printf("$ p: %d\n", a->p);
        printf("$ f: %14.10e\n", a->f);
        printf("$ |c|: %g\n", a->normc);
        printf("$ |T|: %g\n", a->normT);
        printf("$ exit code: %d\n", a->code);
        printf("$ iterations: %d\n", a->niter);
        printf("$ function evaluations: %d\n", a->nfun);
        printf("$ gradient evaluations: %d\n", a->ngrad);
        printf("$ second-order corrections: %d\n", a->nsoc);
        printf("$ QP problems: %d\n", a->numberOfQPSolves);
        //printf("$ relaxed QP problems: %d\n", a->nrelax);
        //printf("$ iterative refinement: %d\n", a->nrefine);
    }
    L1SQPdeleteStorage(a);
}

void L1SQPcomputeSearchDirection(L1SQP a) {

    int i;
    int err = L1SQPsolveQP(a);

    if (a->scqp->refine) {
        a->nrefine++;
    }
    if (err == -1) { // reset the Hessian if is not positive definite
        MatrixSetAllTo(a->H, 0.0);
        for (i = 0; i < a->n; i++) {
            a->H->e[i][i] = 1.0;
        }
        err = L1SQPsolveQP(a);
        if (a->scqp->refine) {
            a->nrefine++;
        }
    }
    if (err != 0) {
    	RuntimeWarning("L1SQPcomputeSearchDirection(): inconsistent QP");
        err = L1SQPsolveRelaxedQP(a);
        if (err != 0) {
            a->done = YES;
            if (a->normc > a->tol) {
                a->code = -2;
            } else {
                a->code = 2;
            }
            return;
        }

        VectorSetEqual(a->kkt, a->df);
        if (a->m > 0) { // kkt += dh^T*lambda
            MatrixGemTv(a->kkt, a->dh, a->lambda, 1.0, 1.0);
        }
        if (a->p > 0) { // kkt += dg^T*mu
            MatrixGemTv(a->kkt, a->dg, a->mu, 1.0, 1.0);
        }

        if (a->scqp->hessianType == FullHessian) {
            MatrixMultVec(a->Hd, a->H, a->d);
        } else {
            for (i = 0; i < a->n; i++) {
                a->Hd->e[i] = -a->kkt->e[i];
            }
        }
        L1SQPcomputeDthetaRelaxed(a); // compute theta0, Dtheta and rho
        a->nd = VectorNorm(a->d);
        a->nx = VectorNorm(a->x);
        return;
    }

    // compute Hd = -kkt, normkkt, normT

    VectorSetEqual(a->kkt, a->df);
    if (a->m > 0) { // kkt += dh^T*lambda
        MatrixGemTv(a->kkt, a->dh, a->lambda, 1.0, 1.0);
    }
    if (a->p > 0) { // kkt += dg^T*mu
        MatrixGemTv(a->kkt, a->dg, a->mu, 1.0, 1.0);
    }
    if (a->scqp->hessianType == FullHessian) {
            MatrixMultVec(a->Hd, a->H, a->d);
    } else {
        for (i = 0; i < a->n; i++) {
            a->Hd->e[i] = -a->kkt->e[i];
        }
    }

    a->normkkt = VectorNormInf(a->kkt);
    a->normT = fmax(a->normkkt, a->normc);

    // compute theta0, Dtheta, rho

    L1SQPcomputeDtheta(a);

    // convergence tests

    if (a->normT <= a->tol) {
        a->done = YES;
        a->code = 0;
        return;
    }
    if ((fabs(a->Dtheta) <= FAC * (1.0 + fabs(a->theta0)))
            && (a->normc <= a->tol)) {
        a->done = YES;
        a->code = 1;
        return;
    }

    a->nd = VectorNorm(a->d);
    a->nx = VectorNorm(a->x);
    if ((a->nd <= FAC * (1.0 + a->nx)) && (a->normc <= a->tol)) {
        a->done = YES;
        a->code = 4;
    }
}

int L1SQPsolveQP(L1SQP a) {
    a->relaxed = NO;
    if (a->scqp->hessianType == FullHessian) {
        a->qp->H = a->H;
    } else if (a->scqp->hessianType == InverseLowerCholeskyFactor) {
        a->scqp->iterativeRefinement = NO;
        a->qp->H = a->_Li;
    }
    a->qp->c = a->df;
    if (a->m > 0) {
        a->qp->A1 = a->dh;
        a->qp->b1 = a->h;
    }
    if (a->p > 0) {
        a->qp->A2 = a->dg;
        a->qp->b2 = a->g;
    }
    a->numberOfQPSolves++;
    return SCQPSolve(a->scqp, a->d, a->lambda, a->mu);
}

void L1SQPcomputeDtheta(L1SQP a) {
    double Dt0, Dt1 = 0.0;
    int i;

    int m = a->m;
    int p = a->p;

    double *h = NULL;
    double *g = NULL;
    double *rho = NULL;
    double *lambda = NULL;
    double *mu = NULL;

    if (m > 0) {
        h = a->h->e;
        rho = a->rho->e;
        lambda = a->lambda->e;
    }
    if (p > 0) {
        g = a->g->e;
        rho = a->rho->e;
        mu = a->mu->e;
    }
    for (i = 0; i < m; i++) {
        rho[i] = fmax(0.5 * (rho[i] + fabs(lambda[i])),
                fabs(lambda[i]));
        Dt1 += rho[i] * fabs(h[i]);
    }
    for (i = 0; i < p; i++) {
        rho[i + m] = fmax(0.5 * (rho[i + m] + fabs(mu[i])), fabs(mu[i]));
        Dt1 += rho[i + m] * fmax(0.0, g[i]);
    }

    a->theta0 = a->f + Dt1;

    Dt1 = VectorDot(a->df, a->d) - Dt1;
    Dt0 = -VectorDot(a->d, a->Hd);
    a->Dtheta = fmin(Dt0, Dt1);
}

int L1SQPsolveRelaxedQP(L1SQP a) {
    double beta_hat = 1.0e6;
    double gamma;
    int i, j;
    int n = a->n, m = a->m, p = a->p;
    int n_ = n + m + p;
    int m_ = m;
    int p_ = 2 * p + m;

    a->relaxed = YES;
    a->nrelax++;
    if (a->rqp == NULL) {
        a->rqp = QPNew();

        a->rqp->H = MatrixNew(n_, n_);
        a->rqp->c = VectorNew(n_);

        a->rqplambda = NULL;

        if (m_ > 0) {
            a->rqp->A1 = MatrixNew(m_, n_);
            a->rqp->b1 = VectorNew(m_);
            a->rqplambda = VectorNew(m_);
        } else {
            a->rqp->A1 = NULL;
            a->rqp->b1 = NULL;
        }
        a->rqp->A2 = MatrixNew(p_, n_);
        a->rqp->b2 = VectorNew(p_);

        a->rqpx = VectorNew(n_);
        a->rqpmu = VectorNew(p_);

        a->rscqp = SCQPNew(a->rqp);
        a->rscqp->hessianType = a->scqp->hessianType;
    } else {
        MatrixSetAllTo(a->rqp->H, 0.0);
        for (i = 0; i < n_; i++) {
            a->rqp->c->e[i] = 0.0;
        }
        if (m_ > 0) {
            MatrixSetAllTo(a->rqp->A1, 0.0);
            for (i = 0; i < m_; i++) {
                a->rqp->b1->e[i] = 0.0;
            }
        }
        MatrixSetAllTo(a->rqp->A2, 0.0);
        for (i = 0; i < p_; i++) {
            a->rqp->b2->e[i] = 0.0;
        }
    }
    // rqp.c = [df; rho]
    for (i = 0; i < n; i++) {
        a->rqp->c->e[i] = a->df->e[i];
    }
    for (i = 0; i < (m + p); i++) {
        a->rqp->c->e[n + i] = a->rho->e[i];
    }
    // rqp.H = diag(H, beta_hat)
    if (a->scqp->hessianType == FullHessian) {
        for (i = 0; i < n; i++) {
            for (j = 0; j < n; j++) {
                a->rqp->H->e[i][j] = a->H->e[i][j];
            }
        }
        for (i = n; i < (n + m + p); i++) {
            a->rqp->H->e[i][i] = beta_hat;
        }
    } else if (a->scqp->hessianType == InverseLowerCholeskyFactor) {
        a->rscqp->iterativeRefinement = NO;
        for (i = 0; i < n; i++) {
            for (j = 0; j < n; j++) {
                a->rqp->H->e[i][j] = a->_Li->e[i][j];
            }
        }
        for (i = n; i < (n + m + p); i++) {
            a->rqp->H->e[i][i] = 1.0 / sqrt(beta_hat);
        }
    }

    for (i = 0; i < m; i++) {
        gamma = -1.0;
        if (a->h->e[i] < 0.0) {
            gamma = 1.0;
        }
        for (j = 0; j < n; j++) {
            a->rqp->A1->e[i][j] = a->dh->e[i][j];
        }
        a->rqp->A1->e[i][n + i] = gamma;
        a->rqp->b1->e[i] = a->h->e[i];
        a->rqp->A2->e[p + i][n + i] = -1.0;
    }
    for (i = 0; i < p; i++) {
        for (j = 0; j < n; j++) {
            a->rqp->A2->e[i][j] = a->dg->e[i][j];
        }
        a->rqp->A2->e[i][n + m + i] = -1.0;
        a->rqp->b2->e[i] = a->g->e[i];
        a->rqp->A2->e[p + m + i][n + m + i] = -1.0;
    }
    a->numberOfQPSolves++;
    return SCQPSolve(a->rscqp, a->rqpx, a->rqplambda, a->rqpmu);
}

void L1SQPcomputeDthetaRelaxed(L1SQP a) {
    int i;
    double u, w, Dt1 = 0.0;
    int n = a->n, m = a->m, p = a->p;

    double *d = a->d->e;

    double *rho = NULL;
    double *h = NULL;
    double *g = NULL;
    double *lambda = NULL;
    double *mu = NULL;

    double *rqpx = a->rqpx->e;
    double *rqplambda = NULL;
    double *rqpmu = NULL;

    if (m > 0) {
        rho = a->rho->e;
        h = a->h->e;
        lambda = a->lambda->e;
        rqplambda = a->rqplambda->e;
    }
    if (p > 0) {
        rho = a->rho->e;
        g = a->g->e;
        mu = a->mu->e;
        rqpmu = a->rqpmu->e;
    }
    // copy d, lambda, mu
    for (i = 0; i < n; i++) {
        d[i] = rqpx[i];
    }
    for (i = 0; i < m; i++) {
        lambda[i] = rqplambda[i];
    }
    for (i = 0; i < p; i++) {
        mu[i] = rqpmu[i];
    }
    for (i = 0; i < m; i++) {
        u = rqpx[i + n];
        if (fabs(h[i]) <= DBL_EPSILON) {
            if (u > 0) {
                rho[i] = fabs(lambda[i]);
            } else {
                rho[i] = fmax(0.5 * (rho[i]
                        + fabs(lambda[i])), fabs(lambda[i]));
            }
        } else if (h[i] > DBL_EPSILON) {
            if (u > h[i]) {
                rho[i] = fabs(lambda[i]);
            } else {
                rho[i] = fmax(0.5 * (rho[i]
                        + fabs(lambda[i])), fabs(lambda[i]));
            }
        } else {
            if (u > -h[i]) {
                rho[i] = fabs(lambda[i]);
            } else {
                rho[i] = fmax(0.5 * (rho[i] + fabs(lambda[i])),
                        fabs(lambda[i]));
            }
        }
        Dt1 += rho[i] * fabs(h[i]);
    }
    for (i = 0; i < p; i++) {
        w = rqpx[i + n + m];
        if (fabs(g[i]) <= DBL_EPSILON) {
            if (w > 0) {
                rho[m + i] = fabs(mu[i]);
            } else {
                rho[m + i] = fmax(0.5 * (rho[m + i] + fabs(mu[i])),
                        fabs(mu[i]));
            }
        } else if (g[i] > DBL_EPSILON) {
            if (w > g[i]) {
                rho[m + i] = fabs(mu[i]);
            } else {
                rho[m + i] = fmax(0.5 * (rho[m + i] + fabs(mu[i])),
                        fabs(mu[i]));
            }
        }
        Dt1 += rho[i + m] * fmax(0.0, g[i]);
    }
    a->theta0 = a->f + Dt1;
    a->Dtheta = -VectorDot(a->d, a->Hd);
}

void L1SQPlineSearch(L1SQP a) {
    int m = a->m, p = a->p;
    double nc0 = a->normc;
    double theta1, fac0, fac1, fac3, n_soc;
    int i, ls_iter;

    a->alpha = 1.0;
    a->do_soc = NO;

    for (ls_iter = 0; ls_iter < MAX_LS_ITER; ls_iter++) {
        VectorSetEqual(a->x, a->x0);
        VectorAxpy(a->x, a->alpha, a->d);
        if (a->do_soc && (ls_iter > 0)) {
            // x = x0 + alpha * d + alpha^2 * d_soc
            VectorAxpy(a->x, a->alpha * a->alpha, a->d_soc);
        }

        a->f = a->nlp->fhg(a->x, a->h, a->g);
        a->nfun++;
        theta1 = a->f;
        for (i = 0; i < m; i++) {
            theta1 += a->rho->e[i] * fabs(a->h->e[i]);
        }
        for (i = 0; i < p; i++) {
            theta1 += a->rho->e[m + i] * fmax(a->g->e[i], 0.0);
        }
        a->nh = 0;
        a->ng = 0;
        if (m > 0) {
            a->nh = VectorNormInf(a->h);
        }
        if (p > 0) {
            a->ng = fmax(VectorMax(a->g), 0.0);
        }
        a->normc = fmax(a->nh, a->ng);

        fac0 = a->alpha * a->nd;
        fac1 = a->theta0 + a->alpha * SIGMA * a->Dtheta;
        fac3 = fmin(FAC * (1.0 + a->nx), 1.0e-15);

        if (fac0 <= fac3) {
            a->done = YES;
            if (a->normc <= a->tol) {
                a->code = 3;
            } else {
                a->code = -3; // step size too small
            }
            return;
        } else if (theta1 <= fac1) {
            return;
        } else if (SOC && (ls_iter == 0) && (a->relaxed == NO)
                && ((a->normc > nc0) || (nc0 <= a->tol))) {

            L1SQPsocStep(a);

            n_soc = VectorNorm(a->d_soc);
            if ((n_soc > 10.0 * DBL_EPSILON) && (n_soc < a->nd)) {
                a->do_soc = YES;
                a->nsoc++;
                VectorAxpy(a->x, 1, a->d_soc);
                a->f = a->nlp->fhg(a->x, a->h, a->g);
                a->nfun++;
                theta1 = a->f;
                for (i = 0; i < m; i++) {
                    theta1 += a->rho->e[i] * fabs(a->h->e[i]);
                }
                for (i = 0; i < p; i++) {
                    theta1 += a->rho->e[m + i] * fmax(a->g->e[i], 0.0);
                }
                a->nh = 0;
                a->ng = 0;
                if (m > 0) {
                    a->nh = VectorNormInf(a->h);
                }
                if (p > 0) {
                    a->ng = fmax(VectorMax(a->g), 0.0);
                }
                a->normc = fmax(a->nh, a->ng);
                if (theta1 <= fac1) {
                    return;
                }
            } else {
                a->do_soc = NO;
            }
        }
        a->alpha *= BETA;
    }
    a->done = YES;
    a->code = -6; // too many line searches
}

void L1SQPsocStep(L1SQP a) {
    int n = a->n, m = a->m;
    int i, j, mhat, i1, i2, j1, j2;
    Vector b = NULL, c = NULL, v = NULL, w = NULL, wv = NULL;
    Matrix Rtranspose = NULL, J2 = NULL;
    for (i = 0; i < n; i++) {
        a->d_soc->e[i] = 0.0;
    }

    mhat = m + a->scqp->n_active;

    if (mhat < 1) {
        return;
    }

    // form b
    b = VectorNew(mhat);
    for (i = 0; i < m; i++) {
        b->e[i] = a->h->e[i];
    }
    for (i = 0; i < a->scqp->n_active; i++) {
        b->e[i + m] = a->g->e[a->scqp->which_constraint[i]];
    }

    // form w
    w = VectorNew(mhat);
    Rtranspose = MatrixNew(mhat, mhat);
    for (i = 0; i < mhat; i++) {
        for (j = i; j < mhat; j++) {
            Rtranspose->e[j][i] = a->scqp->R->e[i][j];
        }
    }
    MatrixSolveLower(Rtranspose, b, w);

    // form c
    c = VectorNew(n);
    if (a->scqp->hessianType == FullHessian) {
        MatrixMultVec(c, a->H, a->d);
    } else if (a->scqp->hessianType == InverseLowerCholeskyFactor) {
        VectorAxpy(c, 1.0, a->kkt);
    }
    VectorAxpy(c, 1.0, a->df);

    // form v
    v = NULL;
    if (mhat < n) {
        v = VectorNew(n - mhat);
        i1 = 0;
        i2 = n - 1;
        j1 = mhat;
        j2 = n-1;
        J2 = MatrixNew(i2 - i1 + 1, j2 - j1 + 1);
        for (i = i1; i <= i2; i++) {
            for (j = j1; j <= j2; j++) {
                J2->e[i - i1][j - j1] = a->scqp->J->e[i][j];
            }
        }
        MatrixMultTVec(v, J2, c);
    }

    // form -wv
    wv = VectorNew(n);
    for (i = 0; i < mhat; i++) {
        wv->e[i] = -w->e[i];
    }
    for (i = mhat; i < n; i++) {
        wv->e[i] = -v->e[i - mhat];
    }

    // form d_soc
    MatrixMultVec(a->d_soc, a->scqp->J, wv);

    if (b != NULL) VectorDelete(b);
    if (c != NULL) VectorDelete(c);
    if (v != NULL) VectorDelete(v);
    if (w != NULL) VectorDelete(w);
    if (wv != NULL) VectorDelete(wv);
    if (Rtranspose != NULL) MatrixDelete(Rtranspose);
    if (J2 != NULL) MatrixDelete(J2);
}

void L1SQPupdate(L1SQP a) {
    double deltaf = fabs(a->f - a->f0);
    //printf("deltaf: %e\n", deltaf);
    if (deltaf <= FAC * (1.0 + fabs(a->f)) && (a->normc <= a->tol)) {
        a->done = YES;
        a->code = 5;
        return;
    }
    if (a->fd_gradient == YES) {
        NLPDfhg(a->nlp, a->x, a->df, a->dh, a->dg);
    } else {
        a->nlp->Dfhg(a->x, a->df, a->dh, a->dg);
    }

    a->ngrad++;

    // Do NOT update the Hessian if,
    // (i) the QP problem is inconsistent (i.e., relaxed == true), or
    // (ii) a second order correction is used (i.e., do_soc == true)
    if (a->relaxed == NO) {
        L1SQPBFGSupdate(a);
        // compute L, Linverse
        if (a->scqp->hessianType == InverseLowerCholeskyFactor) {
            int i, j, err = 0;
            if (CholeskyFactor(a->H, a->_L) != 0) {
                //RuntimeWarning("SCQPcomputeLinverse: unable to compute Cholesky factor of H");
                //return 1;
                err = 1;
            }
            MatrixSetAllTo(a->_Li, 0.0);
            VectorSetAllTo(a->_y, 0.0);
            VectorSetAllTo(a->_d, 0.0);
            for (i = 0; i < a->n; i++) {
                a->_y->e[i] = 1.0;
                err = MatrixSolveLower(a->_L, a->_y, a->_d); // d = L.inverse * y
                if (err != 0) {
                    break;
                }
                for (j = 0; j < a->n; j++) {
                    a->_Li->e[j][i] = a->_d->e[j];
                }
                a->_y->e[i] = 0.0;
            }
            if (err != 0) {
                //printf("*** Hessian == identity ***\n");
                for (j = 0; j < a->n; j++) {
                    a->_Li->e[j][j] = 1.0;
                }
            }
        }
    }

    a->f0 = a->f;
    VectorSetEqual(a->x0, a->x);
    VectorSetEqual(a->df0, a->df);
    if (a->m > 0) {
        MatrixSetEqual(a->dh0, a->dh);
    }
    if (a->p > 0) {
        MatrixSetEqual(a->dg0, a->dg);
    }
}

void L1SQPBFGSupdate(L1SQP a) {
    int n = a->n, m = a->m, p = a->p;
    double *v = a->_v->e;
    double *q = a->_q->e;
    double *Bp = a->_Bp->e;
    double *r = a->_r->e;
    double *x = a->x->e;
    double *x0 = a->x0->e;
    double *df = a->df->e;
    double *df0 = a->df0->e;
    double **H = a->H->e;
    double *lambda = NULL;
    double *mu = NULL;
    double **dh = NULL;
    double **dh0 = NULL;
    double **dg = NULL;
    double **dg0 = NULL;
    double *Hi = NULL;

    double dL0, dL1;
    double s0, s1;
    double _a, _b, pr, gamma1, gamma, pBp, pq;
    int i, j;

    if (m > 0) {
        lambda = a->lambda->e;
        dh = a->dh->e;
        dh0 = a->dh0->e;
    }

    if (p > 0) {
        mu = a->mu->e;
        dg = a->dg->e;
        dg0 = a->dg0->e;
    }
    // v = x - x0
    // dL1 = df + dh^T * lambda + dg^T * mu
    // dL0 = df0 +  dh0^T * lambda + dg0^T * mu
    // q = dL1 - dL0
    // pq = v.dot(q)

    pq = 0.0;
    for (i = 0; i < n; i++) {
        v[i] = x[i] - x0[i];
        dL1 = df[i];
        dL0 = df0[i];
        if (m > 0) {
            s1 = 0.0;
            s0 = 0.0;
            for (j = 0; j < m; j++) {
                s1 += dh[j][i] * lambda[j];
                s0 += dh0[j][i] * lambda[j];
            }
            dL1 += s1;
            dL0 += s0;
        }
        if (p > 0) {
            s1 = 0.0;
            s0 = 0.0;
            for (j = 0; j < p; j++) {
                s1 += dg[j][i] * mu[j];
                s0 += dg0[j][i] * mu[j];
            }
            dL1 += s1;
            dL0 += s0;
        }
        q[i] = dL1 - dL0;
        pq += v[i] * q[i]; // = v.dot(q)
    }

    pBp = 0.0;
    for (i = 0; i < n; i++) { // Bp = H * v
        s0 = 0.0;
        Hi = H[i];
        for (j = 0; j < n; j++) {
            s0 += Hi[j] * v[j];
        }
        Bp[i] = s0;
        pBp += s0 * v[i]; // = v.dot(Bp)
    }

    if (pBp < DBL_EPSILON) {
        return;
    }

    gamma = 1.0;
    if (pq < (0.2 * pBp)) {
        gamma = (0.8 * pBp) / (pBp - pq);
    }

    gamma1 = 1.0 - gamma;
    pr = 0.0;
    for (i = 0; i < n; i++) {
        s0 = gamma * q[i] + gamma1 * Bp[i];
        r[i] = s0;
        pr += s0 * v[i]; // = v.dot(r)
    }

    if (pr < DBL_EPSILON) {
        return;
    }

    _a = -(1.0 / pBp);
    _b = (1.0 / pr);
    for (i = 0; i < n; i++) {
        for (j = i; j < n; j++) {
            Hi = H[i];
            s0 = _a * Bp[i] * Bp[j] + _b * r[i] * r[j];
            Hi[j] += s0;
            if (i != j) {
                H[j][i] += s0;
            }
        }
    }
}

int L1SQPSolve(L1SQP a, Vector x0, Vector lambda0, Vector mu0) {
    // set tol, maxIter, f, g, h, etc.
    if (L1SQPcheckNLP(a, x0, lambda0, mu0) != 0) {
        RuntimeWarning("L1SQPSolve: invalid problem specification");
        return -1;
    }

    for (a->niter = 1; a->niter < a->maxIter; a->niter++) {

        if (a->display) {
            L1SQPshowInfo(a);
        }
        L1SQPcomputeSearchDirection(a);
        if (a->done) {
            L1SQPterminate(a);
            return a->code;
        }
        L1SQPlineSearch(a);
        if (a->done) {
            L1SQPterminate(a);
            return a->code;
        }
        L1SQPupdate(a);
        if (a->done) {
            L1SQPterminate(a);
            return a->code;
        }
    }
    a->code = -4;
    L1SQPterminate(a);
    return a->code;
}
