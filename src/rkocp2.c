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
#include "lu.h"
#include "equations.h"
#include "nlpsolve.h"
#include "ocp.h"
#include "rkocp.h"

void _getRKcoeffRadauIIA3(RKOCP r) {
    // RadauIIA order 3
    r->rkMethodName = "RadauIIA3";
    r->rk_stages = 2;
    r->rk_order = 3;
    r->rk_a = MatrixNew(r->rk_stages, r->rk_stages);
    r->rk_w = MatrixNew(r->rk_stages, r->rk_stages);
    r->rk_b = VectorNew(r->rk_stages);
    r->rk_c = VectorNew(r->rk_stages);
    r->rk_e = VectorNew(r->rk_stages);

    double **rk_a, **rk_w, *rk_b, *rk_c;

    rk_a = r->rk_a->e;
    rk_w = r->rk_w->e;
    rk_b = r->rk_b->e;
    rk_c = r->rk_c->e;
    
    rk_a[0][0] = 5.0/12.0;
    rk_a[0][1] = -1.0/12.0;
    
    rk_a[1][0] = 0.75;
    rk_a[1][1] = 0.25;
    
    rk_w[0][0] = 1.5;
    rk_w[0][1] = 0.5;

    rk_w[1][0] = -4.5;
    rk_w[1][1] = 2.5;

    rk_b[0] = 0.75;
    rk_b[1] = 1.0;
    
    rk_c[0] = 1.0/3.0;
    rk_c[1] = 1.0;
}

void _getRKcoeffRadauIIA5(RKOCP r) {
    // RadauIIA order 5
    r->rkMethodName = "RadauIIA5";
    r->rk_stages = 3;
    r->rk_order = 5;
    r->rk_a = MatrixNew(r->rk_stages, r->rk_stages);
    r->rk_w = MatrixNew(r->rk_stages, r->rk_stages);
    r->rk_b = VectorNew(r->rk_stages);
    r->rk_c = VectorNew(r->rk_stages);
    r->rk_e = VectorNew(r->rk_stages);

    double **rk_a, **rk_w, *rk_b, *rk_c;

    double s6 = sqrt(6.0);

    rk_a = r->rk_a->e;
    rk_w = r->rk_w->e;
    rk_b = r->rk_b->e;
    rk_c = r->rk_c->e;
    
    rk_a[0][0] = (88.0-7.0*s6)/360.0;
    rk_a[0][1] = (296.0-169.0*s6)/1800.0;
    rk_a[0][2] = (-2.0+3.0*s6)/225.0;
    
    rk_a[1][0] = (296.0+169.0*s6)/1800.0;
    rk_a[1][1] = (88.0+7.0*s6)/360.0;
    rk_a[1][2] = (-2.0-3.0*s6)/225.0;
    
    rk_a[2][0] = (16.0-s6)/36.0;
    rk_a[2][1] = (16.0+s6)/36.0;
    rk_a[2][2] = 1.0/9.0;

    rk_w[0][0] = 3.2247448713915885e+00;
    rk_w[0][1] = 1.1678400846904053e+00;
    rk_w[0][2] = -2.5319726474218074e-01;
    
    rk_w[1][0] = -3.5678400846904053e+00;
    rk_w[1][1] = 7.7525512860841150e-01;
    rk_w[1][2] = 1.0531972647421806e+00;
    
    rk_w[2][0] = 5.5319726474218109e+00;
    rk_w[2][1] = -7.5319726474218101e+00;
    rk_w[2][2] = 5.0000000000000000e+00;

    rk_b[0] = (16.0-s6)/36.0;
    rk_b[1] = (16.0+s6)/36.0;
    rk_b[2] = 1.0/9.0;
    
    rk_c[0] = (4.0-s6)/10.0;
    rk_c[1] = (4.0+s6)/10.0;
    rk_c[2] = 1.0;
}

// get the stage values Y_i = y + h *sum_{l = 1, s} a_{i,l} K[l]
void _formIRKstage_values(RKOCP r, Vector y, double h, Vector *Kstage, Matrix Ystage) {
    int n_states = r->n_states;
    int rk_stages = r->rk_stages;
    double **rk_a = r->rk_a->e;
    int i, j, k;
    for (i = 0; i < rk_stages; i++) {
        for (j = 0; j < n_states+1; j++) {
            Ystage->e[i][j] = y->e[j];
        }
        for (k = 0; k < rk_stages; k++) {
            for (j = 0; j < n_states+1; j++) {
                Ystage->e[i][j] += h * rk_a[i][k] * Kstage[k]->e[j];
            }
        }
    }
}

// form H = (1/h) W (x) I - I (x) J
// get the L, U, p factors of H
int _formIRKmatrixZ(RKOCP r, Vector y, Vector u, Vector p, double t, double h, Matrix L, Matrix U, int *P) {
    int n = r->n_states+1;
    int s = r->rk_stages;
    double **W = r->rk_w->e;
    int err = 0;
    int i, j, k, l;
    Matrix H = MatrixNew(n*s, n*s);
    Matrix J = MatrixNew(n, n);
    
    // get Fy, Fu, Fp, Ly, Lu, Lp
    if (r->ocp->Ddifferential_equations == NULL) {
        OCPFDfL(r->ocp, y, u, p, t, r->Fy, r->Fu, r->Fp, r->Ly, r->Lu, r->Lp);
    } else {
        r->ocp->Ddifferential_equations(y, u, p, t, r->Fy, r->Fu, r->Fp, r->Ly, r->Lu, r->Lp);
    }
    
    for (i = 0; i < r->n_states; i ++) {
        for (j = 0; j < r->n_states; j++) {
            J->e[i][j] = r->Fy->e[i][j];
        }
        J->e[r->n_states][i] = r->Ly->e[i];
    }
    
    // H = (1/h) W (x) I - I (x) J
    for (i = 0; i < s; i++) {
        for (j = 0; j < s; j++) {
            double wij_h = W[i][j] / h;
            for (k = 0; k < n; k++) {
                H->e[i*n + k][j*n + k] += wij_h;
            }
        }
        for (k = 0; k < n; k++) {
            for (l = 0; l < n; l++) {
                H->e[i*n + k][i*n + l] -= J->e[k][l];
            }
        }
    }

    err = LUFactor(H, L, U, P);
    
    MatrixDelete(J);
    MatrixDelete(H);
    return err;
}

// given y0, Z0, y1, tau = h1/h0, estimate Z1
// using a Lagrange interpolation formula
void _estimateZ(RKOCP r, Vector y0, Vector *Z0, Vector y1, double tau, Vector *Z1) {
    int s = r->rk_stages;
    int n = r->n_states+1;
    int i, j;
    if (s != 3) {
        for (i = 0; i < s-1; i++) {
            for (j = 0; j < n; j++) {
                Z1[i]->e[j] = Z0[s-1]->e[j];
            }
        }
        return;
    }
    double t0 = 0;
    double t1 = r->rk_c->e[0];
    double t2 = r->rk_c->e[1];
    double t3 = r->rk_c->e[2];
    double L1, L2, L3, t;

    for (i = 0; i < s; i++) {
        t = 1.0 + r->rk_c->e[i]*tau;
        L1 =  (t - t0)*(t - t2)*(t - t3)/((-t0 + t1)*(t1 - t2)*(t1 - t3));
        L2 =  (t - t0)*(t - t1)*(t - t3)/((-t0 + t2)*(-t1 + t2)*(t2 - t3));
        L3 =  (t - t0)*(t - t1)*(t - t2)/((-t0 + t3)*(-t1 + t3)*(-t2 + t3));
        for (j = 0; j < n; j++) {
            Z1[i]->e[j] = y0->e[j] + L1*Z0[0]->e[j] + L2*Z0[1]->e[j] + L3*Z0[2]->e[j] - y1->e[j];
        }
    }
}
 
void _updateZx_stage(RKOCP r, Matrix *Zx_stage, Matrix dZx) {
    int i, j, k;
    for (i = 0; i < r->rk_stages; i++) {
        for (j = 0; j < r->n_states+1; j++) {
            //for (k = 0; k < r->nlp_n; k++) { // FIX: make this i_node dependent
            for (k = 0; k < r->column_index; k++) {
                Zx_stage[i]->e[j][k] += dZx->e[i*(r->n_states+1) + j][k];
            }
        }
    }
}

void _setYx_dot_stage_implicit(RKOCP r, Matrix *Zx_stage, double h, int l, Matrix Yx_dot) {
    int i, j, k;
    double h1 = 1.0/h;
    MatrixSetAllTo(Yx_dot, 0);
    for (i = 0; i < r->rk_stages; i++) {
        for (j = 0; j < r->n_states+1; j++) {
            //for (k = 0; k < r->nlp_n; k++) { // FIX: make this i_node dependent r->column_index
            for (k = 0; k < r->column_index; k++) {
                Yx_dot->e[j][k] += h1*r->rk_w->e[l][i]*Zx_stage[i]->e[j][k];
            }
        }
    }
}

//
void _setYx_stage_implicit(RKOCP r, Matrix Yx, Matrix *Zx_stage, Matrix *Yx_stage) {
    int i, j, k;
    for (k = 0; k < r->rk_stages; k++) {
        for (i = 0; i < r->n_states+1; i++) {
            //for (j = 0; j < r->nlp_n; j++) { // FIX: make this i_node dependent
            for (j = 0; j < r->column_index; j++) {
                Yx_stage[k]->e[i][j] = Yx->e[i][j] + Zx_stage[k]->e[i][j];
            }
        }
    }
}

// FIX: extrapolate to get an estimate of the stage increments
void _setZx_stage(RKOCP r, int i_node, Matrix *Zx_stage) {
    int i, j, k;
    if (i_node == 0) {
        for (i = 0; i < r->rk_stages; i++) {
            for (j = 0; j < r->n_states+1; j++) {
                for (k = 0; k < r->nlp_n; k++) {
                    Zx_stage[i]->e[j][k] = 0;
                }
            }
        }
    } else {
        for (i = 0; i < r->rk_stages-1; i++) {
            for (j = 0; j < r->n_states+1; j++) {
                //for (k = 0; k < r->nlp_n; k++) { // FIX: make this i_node dependent
                for (k = 0; k < r->column_index; k++) {
                    Zx_stage[i]->e[j][k] = Zx_stage[r->rk_stages-1]->e[j][k];
                }
            }
        }                
    }
}

// res = -(Yx_dot_i - (Fy*Yx_i + Fu*Ux_i + Fp*Px_i))
void _setResidual_implicit(RKOCP r, int i_node, int l, Matrix Yx_dot_l, Matrix Yx_stage_l, Vector ys, Vector us, Vector ps, double t, Matrix Kx_stage_l, Matrix res) {
    // K = [Fy * Yx + Fu * Ux + Fp * Px; 
    //      Ly * Yx + Lu * Ux + Lp * Px]
    _setKx_stage( r, Kx_stage_l, Yx_stage_l, i_node, l, ys, us, ps, t);
    int j, k;
    for (j = 0; j < r->n_states+1; j++) {
        //for (k = 0; k < r->nlp_n; k++) { // FIX: make this i_node dependent
        for (k = 0; k < r->column_index; k++) {
            res->e[l*(r->n_states+1) + j][k] = Kx_stage_l->e[j][k] - Yx_dot_l->e[j][k];
        }
    }
}
