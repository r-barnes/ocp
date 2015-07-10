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

void SCQPkkt(SCQP s);
void SCQPJTmultA(SCQP s, double *a, double *y);
void SCQPcheck(SCQP s, int *a, Vector _x, Vector _lambda, Vector _mu);
void SCQPcomputeDeltaEta(SCQP s, double *a, int mhat, double *y);
int  SCQPcomputeLinverse(SCQP s);
void SCQPcsEval(SCQP scqp, double x_i, double x_j);
int  SCQPdropConstraint(SCQP scqp, int n_active);
int  SCQPequality(SCQP s);
void SCQPgetStorage(SCQP s);
int  SCQPinequality(SCQP s);
void SCQPrefinement(SCQP s);
int  SCQPunconstrained(SCQP s);
void SCQPupdateJR(SCQP scqp, int mhat, int add, double *y);
int  SCQPviolatedConstraint(SCQP scqp);

void _printMat(double **M, int n);
void _printIntVec(int *a, int n);

// create a new SCQP data structure
SCQP SCQPNew(QP qp) {
    SCQP s;

    s = (SCQP) malloc(sizeof (struct scqp_));
    s->qp = qp; // the QP problem to solve (input)
    s->iterativeRefinement = YES; // if YES perform iterative refinement if necessary (input)
    s->refine = NO; // if YES iterative refinement was needed (output)
    s->hessianType = FullHessian;

    // private
    s->n = 0;
    s->m = 0;
    s->p = 0;

    s->x = NULL;
    s->lambda = NULL;
    s->mu = NULL;

    s->getStorage = YES;

    s->y = NULL;
    s->d = NULL;
    s->delta = NULL;
    s->eta = NULL;
    s->kappa = NULL;

    s->R = NULL;
    s->L = NULL;
    s->Li = NULL;
    s->J = NULL;

    s->active = NULL;
    s->which_constraint = NULL;

    s->n_active = 0;
    s->n_dropped = 0;

    s->norm_delta = 0.0;
    s->eta_j = 0.0;
    s->t2 = 0.0;

    s->_b = NULL;
    s->_c = NULL;
    s->_e = NULL;
    s->_f = NULL;
    s->_g = NULL;
    s->_h = NULL;
    s->m_hat = 0;

    s->s_max = 0.0;
    s->removed_constraint = NO;
    
    return s;
}

// delete the SCQP data structure
void SCQPDelete(SCQP a) {
    if (a == NULL) {
        return;
    }
    if (a->y != NULL) VectorDelete(a->y);
    if (a->d != NULL) VectorDelete(a->d);
    if (a->delta != NULL) VectorDelete(a->delta);
    if (a->eta != NULL) VectorDelete(a->eta);
    if (a->kappa != NULL) VectorDelete(a->kappa);
    if (a->R != NULL) MatrixDelete(a->R);
    if (a->L != NULL) MatrixDelete(a->L);
    if (a->Li != NULL) MatrixDelete(a->Li);
    if (a->J != NULL) MatrixDelete(a->J);
    if (a->active != NULL) free(a->active);
    if (a->which_constraint != NULL) free(a->which_constraint);
    if (a->_b != NULL) VectorDelete(a->_b);
    if (a->_c != NULL) VectorDelete(a->_c);
    if (a->_e != NULL) VectorDelete(a->_e);
    if (a->_f != NULL) VectorDelete(a->_f);
    if (a->_g != NULL) VectorDelete(a->_g);
    if (a->_h != NULL) VectorDelete(a->_h);
    free(a);
}

// check for consistent matrix/vector dimensions
// input:
// s :- the SCQP data structure
// a :- int[4]
//   on output
//   a[0] = 0 if the dimensions are valid, a[0] != 0 if the dimensions are invalid
//   a[1] = n, the number of unknown variables
//   a[2] = m, the number of equality constraints
//   a[3] = p, the number of inequality constraints
// _x :- the vector of unknowns
// _lambda :- the vector of Lagrange multipliers associated with the equality constraints
// _mu :- the Lagrange multipliers associated with the inequality constraints
void SCQPcheck(SCQP s, int *a, Vector _x, Vector _lambda, Vector _mu) {
    QPcheck(s->qp, a);
    if (a[0] != 0) {
        return;
    }
    if (_x == NULL) {
        a[0] = 1;
        RuntimeWarning("SCQPcheck: x == NULL");
        return;
    }
    if (_x->r != a[1]) {
        a[0] = 1;
        RuntimeWarning("SCQPcheck: x->r does not match the dimensions of the qp problem");
        return;
    }
    if (a[2] > 0) {
        if (_lambda == NULL) {
            a[0] = 1;
            RuntimeWarning("SCQPcheck: lambda == NULL");
            return;
        }
        if (a[2] != _lambda->r) {
            a[0] = 1;
            RuntimeWarning("SCQPcheck: lambda->r does not match the dimensions of the qp problem");
            return;
        }
    }
    if (a[3] > 0) {
        if (_mu == NULL) {
            a[0] = 1;
            RuntimeWarning("SCQPcheck: mu == NULL");
            return;
        }
        if (a[3] != _mu->r) {
            a[0] = 1;
            RuntimeWarning("SCQPcheck: mu->r does not match the dimensions of the qp problem");
            return;
        }
    }
    return;
}

// allocate storage for the SCQP solver
void SCQPgetStorage(SCQP s) {
    int n = s->n;
    int m = s->m;
    int p = s->p;
    int nmp;
    
    s->getStorage = NO;

    if ((m == 0) && (p == 0)) {
        s->d = VectorNew(n);
        return;
    }
    s->y = VectorNew(n);
    s->d = VectorNew(n);
    s->delta = VectorNew(n);
    nmp = m + p;
    nmp = nmp > n ? n : nmp;
    s->eta = VectorNew(nmp);
    s->kappa = VectorNew(nmp);
    s->L = MatrixNew(n, n);
    s->Li = MatrixNew(n, n);
    s->J = MatrixNew(n, n);
    s->R = MatrixNew(n, n);
    if (p > 0) {
        s->active = (int *) calloc(p, sizeof (int)); //new int[p];
        s->which_constraint = (int *) calloc(p, sizeof (int)); //new int[p];
    }
    s->_b = VectorNew(n);
    s->_c = VectorNew(n);
    s->_e = VectorNew(n);
    s->_f = VectorNew(n);
    s->_g = VectorNew(n);
    s->_h = VectorNew(n);
}

// solve the unconstrained SCQP problem
int SCQPunconstrained(SCQP s) {
    int error = 0, i;
    if (s->hessianType == FullHessian) {
        error = CholeskySolveLinearSystem(s->qp->H, s->qp->c, s->d);
        if (error != 0) {
            RuntimeWarning("SCQPunconstrained: H not positive definite");
            return -1;
        }
        for (i = 0; i < s->n; i++) {
            s->x->e[i] = -s->d->e[i];
        }
    } else if (s->hessianType == InverseLowerCholeskyFactor) {
        // d = -Linverse * c
        // x = Linverse.Transpose() * d
        MatrixGemv(s->d, s->qp->H, s->qp->c, -1.0, 0.0);
        MatrixGemTv(s->x, s->qp->H, s->d, 1.0, 0.0);
    }
    return 0;
}

// solve the strictly convex QP problem
// return code:
//  0 :- feasible solution found
// -1 :- H is not positive definite
// -2 :- equality constraints are inconsistent
// -3 :- inequality constraints are inconsistent
// -4 :- too many iterations
int SCQPSolve(SCQP s, Vector x, Vector lambda, Vector mu) {
    int a[4];
    int error = 0;

    s->refine = NO;

    if (s->qp == NULL) {
        return -5;
    }

    SCQPcheck(s, a, x, lambda, mu);

    if (a[0] != 0) {
        return -5;
    }
    s->n = a[1];
    s->m = a[2];
    s->p = a[3];
    s->x = x;
    s->lambda = lambda;
    s->mu = mu;
    // allocate memory
    if (s->getStorage == YES) {
        SCQPgetStorage(s);
    }
    // no constraints
    if ((s->m == 0) && (s->p == 0)) {
        return SCQPunconstrained(s);
    }
    // only equality constraints
    if ((s->m > 0) && (s->p == 0)) {
        return SCQPequality(s);
    }
    // both equality and inequality constraints
    if (s->m > 0) {
        error = SCQPequality(s);
        if (error != 0) {
            return error;
        }
    }
    return SCQPinequality(s);
}

// compute L the Cholesky factor of the Hessian, and L inverse
// return code:
// 0 :- inverse computed
// 1 :- unable to compute L and L inverse 
int SCQPcomputeLinverse(SCQP s) {
    int err = 0;
    int i, j;
    Vector y = s->y;
    Vector d = s->d;
    
    if (s->hessianType == FullHessian) {
        if (CholeskyFactor(s->qp->H, s->L) != 0) {
            //RuntimeWarning("SCQPcomputeLinverse: unable to compute Cholesky factor of H");
            return 1;
        }
        VectorSetAllTo(y, 0.0);
        VectorSetAllTo(d, 0.0);
        for (i = 0; i < s->n; i++) {
            y->e[i] = 1.0;
            err = MatrixSolveLower(s->L, y, d); // d = L.inverse * y
            for (j = 0; j < s->n; j++) {
                s->Li->e[j][i] = d->e[j];
            }
            y->e[i] = 0.0;
        }
    } else if (s->hessianType == InverseLowerCholeskyFactor) {
        for (i = 0; i < s->n; i++) {
            for (j = 0; j <= i; j++) {
                s->Li->e[i][j] = s->qp->H->e[i][j];
            }
        }
    }
    return err;
}

// solve the SCQP problem with only equality constraints
int SCQPequality(SCQP scqp) {
    Vector y = scqp->y;
    Vector d = scqp->d;
    Vector delta = scqp->delta;
    Vector eta = scqp->eta;
    Vector kappa = scqp->kappa;
    Vector x = scqp->x;
    Vector lambda = scqp->lambda;
    Matrix L = scqp->L;
    Matrix Li = scqp->Li;
    Matrix J = scqp->J;
    Matrix R = scqp->R;
    int n = scqp->n;
    int m = scqp->m;
    int p = scqp->p;

    int error;
    int i, j, _j;
    double t1, ad, s, *_v, _vj;

    error = SCQPcomputeLinverse(scqp);
    if (error != 0) {
        return -1;
    }
    if (scqp->hessianType == FullHessian) {
        error = MatrixSolveLower(L, scqp->qp->c, y);
        error = MatrixSolveLowerT(L, y, d);
        for (i = 0; i < n; i++) {
            d->e[i] *= -1.0;
        }
    } else if (scqp->hessianType == InverseLowerCholeskyFactor) {
        MatrixGemv(y, Li, scqp->qp->c, -1.0, 0.0);
        MatrixGemTv(d, Li, y, 1.0, 0.0);
    }

    MatrixTranspose(J, Li);

    for (i = 0; i < m; i++) {
        // computes y, delta, eta, norm_delta, eta_j
        SCQPcomputeDeltaEta(scqp, scqp->qp->A1->e[i], i, y->e); // y = J^T * A1[i,:] is needed later in updatJR

        if (scqp->norm_delta < (10.0 * DBL_EPSILON)) {
            //RuntimeWarning("SCQPequality: inconsistent equality constraints (delta)");
            return -2;
        }
        //ad = Matrix.dot(qp.A1[i], delta);
        //s = Matrix.dot(qp.A1[i], d) + qp.b1[i];
        ad = 0.0;
        s = scqp->qp->b1->e[i];
        {
            _v = scqp->qp->A1->e[i]; // the i-th row of A1
            for (_j = 0; _j < n; _j++) {
                _vj = _v[_j];
                ad += _vj * delta->e[_j];
                s += _vj * d->e[_j];
            }
        }
        t1 = -s / ad;
        VectorAxpy(d, t1, delta); // d += t1 * delta
        for (j = 0; j < i; j++) {
            kappa->e[j] -= t1 * eta->e[j];
        }
        kappa->e[i] = t1;

        SCQPupdateJR(scqp, i, 1, y->e);

        if (fabs(R->e[i][i]) < (10.0 * DBL_EPSILON)) {
            //RuntimeWarning("SCQPequality: inconsistent equality constraints (R)");
            return -2;
        }
    }
    if (p < 1) {
        VectorSetEqual(x, d);
        for (i = 0; i < m; i++) {
            lambda->e[i] = kappa->e[i];
        }
        if (scqp->iterativeRefinement == YES) {
            SCQPrefinement(scqp);
        }
    }
    return 0;
}

// y = J.transpose() * a
void SCQPJTmultA(SCQP s, double *a, double *y) {
    int i, j;
    int n = s->n;
    double **J = s->J->e;
    for (j = 0; j < n; j++) {
        y[j] = 0.0;
    }
    for (i = 0; i < n; i++) {
        double *Ji = J[i];
        double  ai = a[i];
        for (j = 0; j < n; j++) {
            y[j] += Ji[j] * ai;
        }
    }
}

// computes y, delta, eta, norm_delta, eta_j
void SCQPcomputeDeltaEta(SCQP scqp, double *a, int mhat, double *y) {
    int n = scqp->n;
    int m = scqp->m;
    double **J = scqp->J->e;
    double **R = scqp->R->e;
    double *delta = scqp->delta->e;
    double *eta = scqp->eta->e;

    int i, j;
    double sum, abs_sum, n_j, etai;

    // y = J^T * a
    SCQPJTmultA(scqp, a, y);
    scqp->norm_delta = 0.0;
    // delta = -J[:, (m_hat+1) : n] * y[(m_hat+1) : n]
    for (i = 0; i < n; i++) {
        sum = 0.0;
        double *Ji = J[i];
        for (j = mhat; j < n; j++) {
            sum += Ji[j] * y[j];
        }
        delta[i] = -sum;

        abs_sum = fabs(sum);
        if (abs_sum > scqp->norm_delta) {
            scqp->norm_delta = abs_sum;
        }
    }
    // eta = R[1:m_hat, 1:m_hat] \ y[1:m_hat], {R is upper triangular}
    n_j = 0;
    for (i = mhat - 1; i >= 0; i--) {
        sum = 0.0;
        double *Ri = R[i];
        for (j = i + 1; j < mhat; j++) {
            sum += Ri[j] * eta[j];
        }
        etai = (y[i] - sum) / Ri[i];
        eta[i] = etai;
        if (i >= m) {
            if (etai > n_j) {
                n_j = etai;
            }
        }
    }
    scqp->eta_j = n_j;
}

// update the matrix factors J and R
void SCQPupdateJR(SCQP scqp, int mhat, int add, double *y) {
    int n = scqp->n;
    double **R = scqp->R->e;
    double **J = scqp->J->e;
    double *Ri, *Ri1;

    int i, j, k;
    double c, s, a1, a2, a1p, a2p;

    if (add > 0) {
        for (i = 0; i < n; i++) {
            R[i][mhat] = y[i];
        }
        for (j = n - 1; j > mhat; j--) {
            SCQPcsEval(scqp, R[j - 1][mhat], R[j][mhat]);
            c = scqp->cs[0];
            s = scqp->cs[1];
            R[j - 1][mhat] = c * R[j - 1][mhat] + s * R[j][mhat];
            R[j][mhat] = 0.0;
            for (i = 0; i < n; i++) { // do this in parallel
                a1 = J[i][j - 1];
                a2 = J[i][j];
                a1p = a1 * c + a2 * s;
                a2p = -a1 * s + a2 * c;
                J[i][j - 1] = a1p;
                J[i][j] = a2p;
            }
        }
    } else {
        k = -add - 1;
        // Shift the columns of R
        for (i = k; i < mhat; i++) {
            for (j = 0; j <= (i + 1); j++) {
                R[j][i] = R[j][i + 1];
            }
        }
        for (i = k; i < mhat; i++) {
            SCQPcsEval(scqp, R[i][i], R[i + 1][i]);
            c = scqp->cs[0];
            s = scqp->cs[1];
            Ri = R[i];
            Ri1 = R[i+1];
            for (j = i; j < mhat; j++) {
                a1 = Ri[j];
                a2 = Ri1[j];
                a1p = a1 * c + a2*s;
                a2p = -a1 * s + a2*c;
                Ri[j] = a1p;
                Ri1[j] = a2p;
            }
            for (j = 0; j < n; j++) {
                a1 = J[j][i];
                a2 = J[j][i + 1];
                a1p = a1 * c + a2 * s;
                a2p = -a1 * s + a2 * c;
                J[j][i] = a1p;
                J[j][i + 1] = a2p;
            }
        }
    }
}

// compute the Givens matrix coefficients
void SCQPcsEval(SCQP scqp, double x_i, double x_j) {
    double c, s, b;
    if (x_j == 0.0) {
        c = 1.0;
        s = 0.0;
    } else if (fabs(x_j) > fabs(x_i)) {
        b = x_i / x_j;
        s = 1.0 / sqrt(1.0 + b * b);
        c = b * s;
    } else {
        b = x_j / x_i;
        c = 1.0 / sqrt(1.0 + b * b);
        s = b * c;
    }
    scqp->cs[0] = c;
    scqp->cs[1] = s;
}

// find the most 'violated' inequality constraint
// return: constraint number
int SCQPviolatedConstraint(SCQP scqp) {
    int n = scqp->n;
    int p = scqp->p;
    double *d = scqp->d->e;

    double *v;
    int i, j, r;
    double s, s_max;

    r = -1;
    s_max = 0;

    for (i = 0; i < p; i++) {
        if (scqp->active[i] == 0) {
            //s = Matrix.dot(qp.A2[i], d) + qp.b2[i];
            s = scqp->qp->b2->e[i];
            v = scqp->qp->A2->e[i];
            for (j = 0; j < n; j++) {
                s += v[j] * d[j];
            }
            if (s > s_max) {
                s_max = s;
                r = i;
            }
        }
    }
    scqp->s_max = s_max;
    if (s_max <= (10.0 * DBL_EPSILON)) {
        return -1;
    }
    return r;
}

// determine which inequality constraint should be removed from the active set
int SCQPdropConstraint(SCQP scqp, int n_active) {

    double *eta = scqp->eta->e;
    double *kappa = scqp->kappa->e;
    int m = scqp->m;

    double t, etaj;
    int j_drop, j;

    scqp->t2 = 0.0;
    j_drop = -1;
    for (j = m; j < (m + n_active); j++) {
        etaj = eta[j];
        if (etaj > (10.0 * DBL_EPSILON)) {
            t = kappa[j] / etaj;
            if (scqp->t2 == 0.0) {
                scqp->t2 = t;
                j_drop = j;
            } else {
                if (t < scqp->t2) {
                    scqp->t2 = t;
                    j_drop = j;
                }
            }
        }
    }
    return j_drop;
}

// solve the SCQP problem that includes inequality constraints
int SCQPinequality(SCQP scqp) {
    int n = scqp->n;
    int m = scqp->m;
    int p = scqp->p;
    Vector x = scqp->x;
    Vector y = scqp->y;
    Vector d = scqp->d;
    Matrix L = scqp->L;
    Matrix Li = scqp->Li;
    Matrix J = scqp->J;
    Matrix R = scqp->R;
    Vector delta = scqp->delta;
    Vector eta = scqp->eta;
    Vector kappa = scqp->kappa;
    int *active = scqp->active;
    int *which_constraint = scqp->which_constraint;

    int error;
    int i, j, r, mhat, MAX_ITER, iter, j_drop, _j;
    double t1, ad, s, nu_j, nu_bar, *_v, _vj;

    if (m == 0) {
        // compute L and Linverse
        error = SCQPcomputeLinverse(scqp); 
        if (error != 0) {
            return -1;
        }
        // the unconstrained minimum
        if (scqp->hessianType == FullHessian) {
            // FIX: use Li instead?
            error = MatrixSolveLower(L, scqp->qp->c, y); 
            error = MatrixSolveLowerT(L, y, d);
            for (i = 0; i < n; i++) {
                d->e[i] *= -1.0;
            }
        } else if (scqp->hessianType == InverseLowerCholeskyFactor) {
            MatrixGemv(y, Li, scqp->qp->c, -1.0, 0.0);
            MatrixGemTv(d, Li, y, 1.0, 0.0);
        }
        
        MatrixTranspose(J, Li);
    }
    for (i = 0; i < p; i++) {
        active[i] = 0;
        which_constraint[i] = -1;
    }
    scqp->n_active = 0;
    scqp->n_dropped = 0;
    MAX_ITER = 5 * (n + m + p);

    for (iter = 1; iter <= MAX_ITER; iter++) {
        mhat = m + scqp->n_active;
        nu_bar = 0.0;
        nu_j = 0.0;
        t1 = 0.0;

        r = SCQPviolatedConstraint(scqp);

        if (r < 0) {
            break;
        }
        // computes y, delta, eta, norm_delta, eta_j
        SCQPcomputeDeltaEta(scqp, scqp->qp->A2->e[r], mhat, y->e);
        
        if ((scqp->norm_delta < (10.0 * DBL_EPSILON))
                && (scqp->eta_j < (10.0 * DBL_EPSILON))) {
            //RuntimeWarning("inconsistent constraints (inequality delta)");
            return -3;
        }
        if (scqp->norm_delta >= (10.0 * DBL_EPSILON)) {
            //ad = Matrix.dot(qp.A2[r], delta);
            //s = Matrix.dot(qp.A2[r], d) + qp.b2[r];
            ad = 0.0;
            s = scqp->qp->b2->e[r];
            {
                _v = scqp->qp->A2->e[r];
                //double _vj;
                for (_j = 0; _j < n; _j++) {
                    _vj = _v[_j];
                    ad += _vj * delta->e[_j];
                    s += _vj * d->e[_j];
                }
            }
            t1 = -s / ad;
            nu_j = 0.0;
            for (j = m; j < (m + scqp->n_active); j++) {
                nu_j = kappa->e[j] - t1 * eta->e[j];
                if (nu_j < 0.0) {
                    break;
                }
            }
        }
        while ((scqp->norm_delta < (10.0 * DBL_EPSILON))
                || (nu_j < 0.0)) { // drop a constraint
            // compute j_drop and t2
            j_drop = SCQPdropConstraint(scqp, scqp->n_active);
            if (j_drop < m) {
                break;
            }
            scqp->n_dropped++;

            VectorAxpy(d, scqp->t2, delta);
            for (j = 0; j < (m + scqp->n_active); j++) {
                kappa->e[j] -= scqp->t2 * eta->e[j];
            }
            nu_bar += scqp->t2;

            // shift kappa
            if (j_drop < (m + scqp->n_active - 1)) {
                for (i = j_drop; i < (m + scqp->n_active - 1); i++) {
                    kappa->e[i] = kappa->e[i + 1];
                }
            }
            // de-activate the dropped constraint
            active[which_constraint[j_drop - m]] = 0;
            // shift which_constraint
            if (j_drop < (m + scqp->n_active - 1)) {
                for (i = (j_drop - m); i < (scqp->n_active - 1); i++) {
                    which_constraint[i] = which_constraint[i + 1];
                }
            }
            scqp->n_active--;
            mhat = m + scqp->n_active;

            if ((m + scqp->n_active) > 0) {
                SCQPupdateJR(scqp, mhat, -(j_drop + 1), NULL);
            } else {
                MatrixSetAllTo(R, 0.0);
                for (i = 0; i < kappa->r; i++) {
                    kappa->e[i] = 0.0;
                }
                MatrixTranspose(J, Li);
            }

            SCQPcomputeDeltaEta(scqp, scqp->qp->A2->e[r], mhat, y->e);

            t1 = 0.0;
            nu_j = 0.0;

            if ((scqp->norm_delta < (10.0 * DBL_EPSILON))
                    && (scqp->eta_j < (10.0 * DBL_EPSILON))) {
                //RuntimeWarning("inconsistent constraints (inequality delta) 2");
                return -3;
            }
            if (scqp->norm_delta >= (10.0 * DBL_EPSILON)) {
                //ad = Matrix.dot(qp.A2[r], delta);
                //s = Matrix.dot(qp.A2[r], d) + qp.b2[r];
                ad = 0.0;
                s = scqp->qp->b2->e[r];
                {
                    _v = scqp->qp->A2->e[r];
                    // double _vj;
                    for (_j = 0; _j < n; _j++) {
                        _vj = _v[_j];
                        ad += _vj * delta->e[_j];
                        s += _vj * d->e[_j];
                    }
                }
                t1 = -s / ad;
                nu_j = 0.0;
                for (j = m; j < (m + scqp->n_active); j++) {
                    nu_j = kappa->e[j] - t1 * eta->e[j];
                    if (nu_j < 0.0) {
                        break;
                    }
                }
            }
        }
        VectorAxpy(d, t1, delta);
        for (j = 0; j < (m + scqp->n_active); j++) {
            kappa->e[j] -= t1 * eta->e[j];
        }
        kappa->e[m + scqp->n_active] = t1 + nu_bar;
        active[r] = 1;
        which_constraint[scqp->n_active] = r;

        SCQPupdateJR(scqp, m + scqp->n_active, 1, y->e);

        scqp->n_active += 1;
    }

    if (iter >= MAX_ITER) {
        //RuntimeWarning("SCQPinequality: too many iterations");
        return -4;
    }

    VectorSetEqual(x, d);

    if (m > 0) {
        for (i = 0; i < m; i++) {
            scqp->lambda->e[i] = kappa->e[i];
        }
    }
    for (i = 0; i < p; i++) {
        scqp->mu->e[i] = 0.0;
    }
    if (scqp->n_active > 0) {
        for (i = 0; i < scqp->n_active; i++) {
            j = scqp->which_constraint[i];
            scqp->mu->e[j] = fmax(kappa->e[m + i], 0.0);
        }
    }
    if (scqp->iterativeRefinement == YES) {
        SCQPrefinement(scqp);
    }

    return 0;
}

// perform iterative refinement of the solution
void SCQPrefinement(SCQP s) {
    int n = s->n;
    int m = s->m;
    int p = s->p;
    int n_active = s->n_active;

    Vector w = NULL, ed = NULL, ek = NULL, uv = NULL, wv = NULL,
            c0 = NULL, c1 = NULL, c2 = NULL, rd = NULL, rk = NULL;
    Matrix Rmhat;
    int mhat = m + n_active;
    int i, j, k;
    double norm_rd, norm_rk, z;

    rd = VectorNew(n);
    if (mhat > 0) {
        rk = VectorNew(mhat);
    }
    // form the residual rd = H*d + A1^T*lambda + A2^T*mu + c
    c0 = VectorNew(n);
    MatrixMultVec(c0, s->qp->H, s->x); // c0 = H*x
    for (i = 0; i < n; i++) {
        rd->e[i] = s->qp->c->e[i] + c0->e[i];
    }
    if (m > 0) {
        c1 = VectorNew(n);
        MatrixMultTVec(c1, s->qp->A1, s->lambda); // c1 = A1^T*lambda
        for (i = 0; i < n; i++) {
            rd->e[i] += c1->e[i];
        }
    }
    if (p > 0) {
        c2 = VectorNew(n);
        MatrixMultTVec(c2, s->qp->A2, s->mu); // c2 = A2^T*mu
        for (i = 0; i < n; i++) {
            rd->e[i] += c2->e[i];
        }
    }
    norm_rd = VectorNorm(rd);

    // form the constraint residual
    for (i = 0; i < m; i++) {
        //rk->e[i] = MatrixDot(qp->A1[i], x) + s->qp->b1->e[i];
        z = s->qp->b1->e[i];
        for (j = 0; j < n; j++) {
            z += s->qp->A1->e[i][j] * s->x->e[j];
        }
        rk->e[i] = z;
    }
    for (i = 0; i < n_active; i++) {
        j = s->which_constraint[i];
        //rk[m + i] = MatrixDot(qp->A2[j], x) + s->qp->b2[j];
        z = s->qp->b2->e[j];
        for (k = 0; k < n; k++) {
            z += s->qp->A2->e[j][k] * s->x->e[k];
        }
        rk->e[m + i] = z; //fmax(z, 0.0); // Use z here???
    }
    norm_rk = 0.0;
    if (mhat > 0) {
        norm_rk = VectorNorm(rk);
    }
    if ((norm_rd < SQRT_EPS) && (norm_rk < SQRT_EPS)) {
        VectorDelete(c0);
        VectorDelete(rd);
        if (m > 0) {
            VectorDelete(c1);
        }
        if (p > 0) {
            VectorDelete(c2);
        }
        if (mhat > 0) {
            VectorDelete(rk);
        }
        return;
    }
    s->refine = YES;
    uv = VectorNew(n);
    wv = VectorNew(n);
    MatrixMultTVec(uv, s->J, rd);
    for (i = 0; i < n; i++) {
        wv->e[i] = uv->e[i];
    }

    w = NULL;
    Rmhat = NULL;

    if (mhat > 0) {
        w = VectorNew(mhat);
        Rmhat = MatrixNew(mhat, mhat);
        MatrixGetSlice(s->R, 0, mhat-1, 0, mhat-1, Rmhat);
        MatrixSolveUpperT(Rmhat, rk, w);
        for (i = 0; i < mhat; i++) {
            wv->e[i] = w->e[i];
        }
    }
    ed = VectorNew(n);
    MatrixMultVec(ed, s->J, wv);

    for (i = 0; i < n; i++) {
        s->x->e[i] -= ed->e[i];
    }
    if (mhat > 0) {
        for (i = 0; i < mhat; i++) {
            w->e[i] -= uv->e[i];
        }
        ek = VectorNew(mhat);
        MatrixSolveUpper(Rmhat, w, ek);
        for (i = 0; i < m; i++) {
            s->lambda->e[i] += ek->e[i];
        }
        for (i = 0; i < n_active; i++) {
            j = s->which_constraint[i];
            s->mu->e[j] = fmax(s->mu->e[j] + ek->e[m + i], 0.0);
        }
    }
    if (uv != NULL) VectorDelete(uv);
    if (wv != NULL) VectorDelete(wv);
    if (ed != NULL) VectorDelete(ed);
    if (ek != NULL) VectorDelete(ek);
    if (w != NULL) VectorDelete(w);
    if (c0 != NULL) VectorDelete(c0);
    if (c1 != NULL) VectorDelete(c1);
    if (c2 != NULL) VectorDelete(c2);
    if (rd != NULL) VectorDelete(rd);
    if (rk != NULL) VectorDelete(rk);
    if (Rmhat != NULL) MatrixDelete(Rmhat);
}

// compute the KKT conditions for the QP problem
void SCQPkkt(SCQP s) {
    // r0 = H*x + A1^T*lambda + A2^T*mu + c
    // r1 = A1*x + b1
    // r2 = mu*(A2^T*x + b2)
    // r3 = max(0, A2^T*x + b2)
    int i, n, m, p;
    Vector r3 = NULL, r2 = NULL, r1 = NULL, r0 = NULL, c0 = NULL, c1 = NULL, c2 = NULL;

    QP qp = s->qp;
    Vector qpx = s->x;
    Vector qplambda = s->lambda;
    Vector qpmu = s->mu;

    n = s->n;
    m = s->m;
    p = s->p;

    r0 = VectorNew(n);
    c0 = VectorNew(n);

    if (p > 0) {
        c2 = VectorNew(n);
        r2 = VectorNew(p);
        r3 = VectorNew(p);
        MatrixMultVec(r2, qp->A2, qpx);
        for (i = 0; i < p; i++) {
            r2->e[i] += qp->b2->e[i];
            r3->e[i] = fmax(0.0, r2->e[i]);
            r2->e[i] *= qpmu->e[i];
        }
        printf("||r2|| = %g\n", VectorNorm(r2));
        printf("||r3|| = %g\n", VectorNorm(r3));
        MatrixMultTVec(c2, qp->A2, qpmu);
        for (i = 0; i < n; i++) {
            r0->e[i] += c2->e[i];
        }
    }
    if (m > 0) {
        c1 = VectorNew(n);
        r1 = VectorNew(m);
        MatrixMultVec(r1, qp->A1, qpx);
        for (i = 0; i < m; i++)
            r1->e[i] += qp->b1->e[i];
        printf("||r1|| = %g\n", VectorNorm(r1));
        MatrixMultTVec(c1, qp->A1, qplambda);
        for (i = 0; i < n; i++) {
            r0->e[i] += c1->e[i];
        }
    }
    MatrixMultVec(c0, qp->H, qpx);
    for (i = 0; i < n; i++) {
        r0->e[i] += qp->c->e[i] + c0->e[i];
    }
    printf("||r0|| = %g\n", VectorNorm(r0));

    VectorDelete(r0);
    VectorDelete(c0);
    if (m > 0) {
        VectorDelete(c1);
        VectorDelete(r1);
    }
    if (p > 0) {
        VectorDelete(c2);
        VectorDelete(r2);
        VectorDelete(r3);
    }
}
