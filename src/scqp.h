// Copyright (c) 2009,2010,2011,2012 Brian C. Fabien
// All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the license.txt file.

/* 
 * File:   scqp.h
 * Author: fabien
 *
 * Created on July 9, 2010, 6:28 PM
 */

#ifndef _SCQP_H
#define	_SCQP_H

#ifdef	__cplusplus
extern "C" {
#endif
    
enum HessianType { FullHessian, LowerCholeskyFactor, 
                   InverseLowerCholeskyFactor };

struct scqp_ {
    QP qp;
    int refine;
    int iterativeRefinement;
    enum HessianType hessianType;

    // private
    
    int n, m, p;
    Vector x, lambda, mu;

    Vector y, d, delta, eta, kappa;
    Matrix R, L, Li, J;

    int *active, *which_constraint;
    int n_active, n_dropped, m_hat;

    double norm_delta, eta_j, t2;
    double cs[2];
    
    int getStorage;

    Vector _b, _c, _e, _f, _g, _h;

    double s_max;
    int removed_constraint;
};

typedef struct scqp_ * SCQP;

SCQP SCQPNew(QP qp);
void SCQPDelete(SCQP a);
int  SCQPSolve(SCQP a, Vector x, Vector lambda, Vector mu);

//int  SCQPSolve1(SCQP a, Vector x, Vector lambda, Vector mu);
//void SCQPsolveKKT1(SCQP scqp, double *rx, double *ry);

//int  SCQPSolve2(SCQP a, Vector x, Vector lambda, Vector mu);
//void SCQPsolveKKT2(SCQP scqp, double *rx, double *ry);

//int  SCQPSolve3(SCQP a, Vector x, Vector lambda, Vector mu);
//void SCQPsolveKKT3(SCQP scqp, double *rx, double *ry);
#ifdef	__cplusplus
}
#endif

#endif	/* _SCQP_H */
