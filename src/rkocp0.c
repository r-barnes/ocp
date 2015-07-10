// Copyright (c) 2009,2010,2011,2012 Brian C. Fabien
// All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the license.txt file.

// TODO: provide rkocp.Nnodes?
// TODO: readData, writeData methods?
// TODO: use Butcher transformation for the implicit Runge-Kutta methods
// TODO: error estimator for the implicit methods

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

static RKOCP thisRKOCP;

// create a RKOCP data structure
RKOCP RKOCPNew(OCP o) {
    RKOCP r;
    r = (RKOCP) malloc(sizeof (struct rkocp_));
    r->ocp = o;
    r->T = NULL;
    r->P = NULL;
    r->Y = NULL;
    r->U = NULL;
    r->Yw = NULL;
    r->display_data = NO;
    r->nlp_n = 0;
    r->nlp_m = 0;
    r->nlp_p = 0;
    r->prob = NULL;
    r->n_nodes = 0;
    r->n_states = 0;
    r->n_parameters = 0;
    r->n_controls = 0;
    r->n_initial = 0;
    r->n_terminal = 0;
    r->n_inequality = 0;
    r->c_type = 0;
    r->yc = NULL;
    r->uc = NULL;
    r->pc = NULL;
    r->Gamma = NULL;
    r->Psi = NULL;
    r->F = NULL;
    r->G = NULL;
    r->Fy = NULL;
    r->Fu = NULL;
    r->Fp = NULL;
    r->Gy = NULL;
    r->Gu = NULL;
    r->Gp = NULL;
    r->Gammay = NULL;
    r->Gammap = NULL;
    r->Psiy = NULL;
    r->Psip = NULL;
    r->Ly = NULL;
    r->Lu = NULL;
    r->Lp = NULL;
    r->phiy = NULL;
    r->phip = NULL;
    
    // Integration data
    // -: f, h, g calculation
    r->Kstage = NULL; // stage derivatives
    r->Ydata = NULL; // stage values for the state
    r->Udata = NULL; // stage values for the controls
    r->lte = NULL;
    r->max_lte = 0.0;
    r->local_error = NULL;
    r->only_compute_lte = NO;
    
    // -: sensitivity calculation, df, dh, dg
    r->Yx = NULL; // state sensitivity
    r->Yx_stage = NULL; // sensitivity stage values
    r->Kx_stage = NULL; // sensitivity stage derivatives
    
    // Runge-Kutta
    r->rk_method = 0;
    r->rk_stages = 0;
    r->rk_order = 0;
    r->rk_b = NULL;
    r->rk_c = NULL;
    r->rk_e = NULL;
    r->rk_a = NULL;
    r->rk_w = NULL; // w = inv(a)
    
    // SQP
    r->x0 = NULL;
    r->lambda0 = NULL;
    r->mu0 = NULL;
    r->tol = 0.0;
    r->max_iter = 0;
    
    // control interpoaltion coefficients
    r->Ucubic = NULL;

    // Implicit Runge-Kutta methods
    r->Lnode = NULL;
    r->Unode = NULL;
    r->pLUnode = NULL;
    r->integration_fatal_error = NO;
    r->column_index = 0;
    
    r->_memory_allocated = NO;

    // Input
    r->control_type = LINEAR;       // Control type: LINEAR (default)
    r->runge_kutta_method = DormandPrince8_7; // Runge-Kutta method: DormandPrince8_7 (default) 
    r->controlTypeName = "LINEAR";
    r->rkMethodName = "Dormand-Prince 8(7)";
    r->tolerance = 1.0e-6;          // Convergence toerance 1.0e-6 (default) 
    r->maximum_iterations = 5000;   // Maximum number of iterations: 5000 (default) 
    r->nlp_method = L1SQPmethod;    // NLP solver method: L1SQPmethod (default)
    
    thisRKOCP = r;

    return r;
}

double _RKOCPfhg(Vector nlp_x, Vector nlp_h, Vector nlp_g) {
    RKOCP r = thisRKOCP;
    double phi = 0;
    double t;
    int i_node, l;
    int i, j;
    int m = r->nlp_m;
    int p = r->nlp_p;
    int n_nodes = r->n_nodes;
    int n_parameters = r->n_parameters;
    int n_states = r->n_states;
    int n_controls = r->n_controls;
    int n_inequality = r->ocp->n_inequality;
    int n_initial = r->ocp->n_initial;
    int n_terminal = r->ocp->n_terminal;
    int only_compute_lte = r->only_compute_lte;

    int rk_stages = r->rk_stages;
    double **rk_a = r->rk_a->e;
    double *rk_b = r->rk_b->e;
    double *rk_c = r->rk_c->e;
    double *rk_e = r->rk_e->e;

    Vector T = r->T;
    Matrix Yw = r->Yw;
    Vector yc = r->yc;
    Vector uc = r->uc;
    Vector pc = r->pc;
    Vector local_error = r->local_error;
    Vector lte = r->lte;

    Vector Gamma = r->Gamma;
    Vector G = r->G;
    Vector Psi = r->Psi;

    Matrix *Ydata = r->Ydata;
    Matrix *Udata = r->Udata;
    Vector *Kstage = r->Kstage;

    double _y;

    if (only_compute_lte == NO) {
        if (m > 0) {
            VectorSetAllTo(nlp_h, 0.0);
        }
        if (p > 0) {
            VectorSetAllTo(nlp_g, 0.0);
        }
    }
    for (i = 0; i < n_parameters; i++) {
        pc->e[i] = nlp_x->e[i];

    }
    for (j = 0; j < n_states; j++) {
        Yw->e[0][j] = nlp_x->e[n_parameters+j];
        //MAT(Yw,0,j) = VEC(nlp_x,n_parameters+j);
    }
    Yw->e[0][n_states] = 0;

    if (n_controls > 0) {
        _setUdata(r, nlp_x);
    }

    // integrate the ODEs using an explicit Runge-Kutta method

    Matrix Ystage;
    double delta, ti;

    for (i_node = 0; i_node < n_nodes-1; i_node++) {
        Ystage = Ydata[i_node]; // save the stage values for use in Dfhg

        ti = T->e[i_node];
        delta = T->e[i_node+1] - ti;
        t = ti;

        //- save the control values for use in Dfhg
        //- getControl(i_node, t, Udata[i_node][0]);
        for (i = 0; i < n_controls; i++) 
            uc->e[i] = Udata[i_node]->e[0][i];

        for (i = 0; i < n_states + 1; i++) {
            _y = Yw->e[i_node][i];
            Ystage->e[0][i] = _y;
            yc->e[i] = _y;
        }

        if (only_compute_lte == NO) {
            // form Gamma(y(t_initial), p)
            if ((i_node == 0) && (n_initial > 0)) {
                r->ocp->initial_constraints(yc, pc, Gamma);
                for (i = 0; i < n_initial; i++) {
                    nlp_h->e[i] = Gamma->e[i];
                }
            }
            // form g(y,u,p,t)
            if (n_inequality > 0) {
                r->ocp->inequality_constraints(yc, uc, pc, t, G);
                for (i = 0; i < n_inequality; i++) {
                    nlp_g->e[i + i_node * n_inequality] = G->e[i];
                }
            }
        }

        // K = F(y,u,p,t)
        Kstage[0]->e[n_states] = r->ocp->differential_equations(yc, uc, pc, t, Kstage[0]);

        for (l = 1; l < rk_stages; l++) {
            t = ti + rk_c[l]*delta;

            //getControl(i_node, t, Udata[i_node][l]);
            for (i = 0; i < n_controls; i++) 
                uc->e[i] = Udata[i_node]->e[l][i];

            //Y_l = yc + delta * sum_{j=0 to j=l-1} rk_a[l][j] * K_j
            for (i = 0; i < n_states + 1; i++) 
                Ystage->e[l][i] = Yw->e[i_node][i];
            for (j = 0; j < l; j++) {
                for (i = 0; i < n_states+1; i++) {
                    Ystage->e[l][i] += delta*rk_a[l][j]*Kstage[j]->e[i];
                }
            }

            for (i = 0; i < n_states+1; i++) 
                yc->e[i] = Ystage->e[l][i];
                
            Kstage[l]->e[n_states] = r->ocp->differential_equations(yc, uc, pc, t, Kstage[l]);
        }

        if (only_compute_lte == YES) {
            VectorSetAllTo(local_error, 0.0);
        }

        // Yw[i_node+1] = Yw[i_node] + delta * sum_{l=0 to stages-1} b[l] * K[l]
        for (i = 0; i < n_states+1; i++) {
            Yw->e[i_node+1][i] = Yw->e[i_node][i];
        }
        for (l = 0; l < rk_stages; l++) {
            for (i = 0; i < n_states + 1; i++) {
                Yw->e[i_node+1][i] += delta*rk_b[l]*Kstage[l]->e[i];
                if (only_compute_lte) {
                    local_error->e[i] += delta*rk_e[l]*Kstage[l]->e[i];
                }
            }
        }
        if (only_compute_lte == YES) {
            lte->e[i_node + 1] = _compute_LTE(r, local_error->e, Yw->e[i_node], Yw->e[i_node + 1]);
        }
    }

    if (only_compute_lte == YES) {
        return 0;
    }

    // form g at time t_final
    t = T->e[n_nodes-1];
    if (n_controls > 0) {
        //getControl(n_nodes - 2, t, Udata[n_nodes - 1][0]);
        for (i = 0; i < n_controls; i++) 
            uc->e[i] = Udata[n_nodes-2]->e[rk_stages-1][i];
    }

    for (i = 0; i < n_states; i++) 
        yc->e[i] = Yw->e[n_nodes-1][i];
    
    if (n_inequality > 0) {
        r->ocp->inequality_constraints(yc, uc, pc, t, G);
        for (i = 0; i < n_inequality; i++) {
            nlp_g->e[i + (n_nodes-1)*n_inequality] = G->e[i];
        }
    }
    
    // form phi, Psi
    phi = r->ocp->terminal_constraints(yc, pc, Psi);
    for (i = 0; i < n_terminal; i++) {
        nlp_h->e[i + n_initial] = Psi->e[i];
    }

    return phi + Yw->e[n_nodes-1][n_states];
}

void _RKOCPDfhg(Vector nlp_x, Vector nlp_df, Matrix nlp_dh, Matrix nlp_dg) {
/*
 TODO: Instead of copying data to ys and us
 set ys->r = n_states+1, ys->e = state data
 set us->r = n_controls, us->e = control data
*/
    // nlp_x: vec[P, Y(0), U(0, :), U(1, :), ..., U(n_nodes-1, :)]
    // Yx: state+cost sensitivity with respect to nlp_x
    RKOCP r = thisRKOCP;
    int i_node, i, l;
    double t, ti, delta;
    int m = r->nlp_m;
    int p = r->nlp_p;
    int n_nodes = r->n_nodes;
    int n_parameters = r->n_parameters;
    int n_states = r->n_states;
    int n_controls = r->n_controls;
    int n_inequality = r->ocp->n_inequality;
    int n_initial = r->ocp->n_initial;
    int n_terminal = r->ocp->n_terminal;

    int rk_stages = r->rk_stages;
    double *rk_c = r->rk_c->e;

    Vector T = r->T;
    Matrix Yw = r->Yw;
    Matrix Yx = r->Yx;
    Matrix *Ydata = r->Ydata;
    Matrix *Udata = r->Udata;
    Matrix *Yx_stage = r->Yx_stage;
    Matrix *Kx_stage = r->Kx_stage;
    Vector ys = r->yc;
    Vector us = r->uc;
    Vector ps = r->pc;

    VectorSetAllTo(nlp_df, 0.0);

    if (m > 0) {
        MatrixSetAllTo(nlp_dh, 0.0);
    }

    if (p > 0) {
        MatrixSetAllTo(nlp_dg, 0.0);
    }

    if (n_controls > 0) {
        us = VectorNew(n_controls);
    }

    // parameters: ps
    if (n_parameters > 0) {
        ps = VectorNew(n_parameters);
        for (i = 0; i < n_parameters; i++) 
            ps->e[i] = nlp_x->e[i];
    }

    // sensitivity: Yx(t_initial)
    MatrixSetAllTo(Yx, 0.0);
    for (i = 0; i < n_states; i++) {
        Yx->e[i][i+n_parameters] = 1.0;
    }

    for (i_node = 0; i_node < n_nodes-1; i_node++) {
        ti = T->e[i_node];
        delta = T->e[i_node+1] - ti;
        t = ti;

        // states: ys(t)
        for (i = 0; i < n_states+1; i++) 
            ys->e[i] = Yw->e[i_node][i];

        // controls: us(t)
        for (i = 0; i < n_controls; i++) 
            us->e[i] = Udata[i_node]->e[0][i];
        

        // dGamma -> nlp_dh
        if ((i_node == 0) && (n_initial > 0)) {
            _setDhGamma(r, ys, ps, nlp_dh);
        }

        // dg(t) -> nlp_dg
        if (n_inequality > 0) {
            _setDg(r, Yx, ys, us, ps, t, i_node, nlp_dg);
        }

        /*
        Integrate the sensitivities using an explicit Runge-Kutta method

        for l = 0 to (stages-1):
        Yx_stage[l] = Yx[i_node] + delta_{i_node} * Sum_{r = 0 to l-1} a(l,r) * Kx_stage[r]
        Kx_stage[l] = [Fy*Yx_stage[l] + Fu*u_x + Fp*P_x; Ly*Yx_stage[l] + Lu*u_x + Lp*P_x]
        Yx[i_node+1] = Yx[i_node] + delta_{i_node} * Sum_{r = 0 to stages-1} b(r) * Kx_stage[r]
         */

        l = 0;

        for (i = 0; i < n_states+1; i++) 
            ys->e[i] = Ydata[i_node]->e[l][i];
        MatrixSetEqual(Yx_stage[l], Yx);
        _setKx_stage(r, Kx_stage[l], Yx_stage[l], i_node, l, ys, us, ps, t);
        for (l = 1; l < rk_stages; l++) {
            t = ti + rk_c[l]*delta;
            _setYx_stage(r, Yx_stage[l], Yx, i_node, l, delta, Kx_stage);
            for (i = 0; i < n_controls; i++) 
                us->e[i] = Udata[i_node]->e[l][i];
            for (i = 0; i < n_states+1; i++) 
                ys->e[i] = Ydata[i_node]->e[l][i];
            _setKx_stage(r, Kx_stage[l], Yx_stage[l], i_node, l, ys, us, ps, t);
        }
        // => Yx += delta_{i_node} * Sum_{r = 0 to stages-1} b(r) * Kx_stage[r]
        _setYx(r, Yx, Kx_stage, i_node, delta);
    }

    i_node = n_nodes-1;
    t = T->e[i_node];

    // states: ys(t)
    for (i = 0; i < n_states+1; i++) 
        ys->e[i] = Yw->e[i_node][i];
    //ys = Yw->e[i_node];

    // controls: us(t)
    for (i = 0; i < n_controls; i++) 
        us->e[i] = Udata[n_nodes-2]->e[rk_stages-1][i];
    
    // dg(t_{i_node}) -> nlp_dg
    if (n_inequality > 0) {
        _setDg(r, Yx, ys, us, ps, t, i_node, nlp_dg);
    }

    // Terminal conditions Psi and penalty phi
    if (r->ocp->Dterminal_constraints == NULL) {
        OCPFDphiPsi(r->ocp, ys, ps, r->phiy, r->phip, r->Psiy, r->Psip);
    } else {
        r->ocp->Dterminal_constraints(ys, ps, r->phiy, r->phip, r->Psiy, r->Psip);
    }
    // dPsi -> nlp_dh
    if (n_terminal > 0) {
        _setDhPsi(r, Yx, nlp_dh);
    }

    // df = integral of L
    for (i = 0; i < r->nlp_n; i++) 
        nlp_df->e[i] = Yx->e[n_states][i];

    // df += dphi
    _addDphi(r, Yx, nlp_df);
}

double _RKOCPfhg_implicitZ(Vector nlp_x, Vector nlp_h, Vector nlp_g) {
    RKOCP r = thisRKOCP;
    double phi = 0;
    int i_node;
    int i, j, k;
    int m = r->nlp_m;
    int p = r->nlp_p;
    int n_nodes = r->n_nodes;
    int n_parameters = r->n_parameters;
    int n_states = r->n_states;
    int n_controls = r->n_controls;
    int n_inequality = r->ocp->n_inequality;
    int n_initial = r->ocp->n_initial;
    int n_terminal = r->ocp->n_terminal;
    int only_compute_lte = r->only_compute_lte;

    int rk_stages = r->rk_stages;
    double **W = r->rk_w->e;
    //double **rk_a = r->rk_a->e;
    //double *rk_b = r->rk_b->e;
    double *rk_c = r->rk_c->e;
    //double *rk_e = r->rk_e->e;

    Vector T = r->T;
    Matrix Yw = r->Yw;
    Vector yc = r->yc;
    Vector uc = r->uc;
    Vector pc = r->pc;
    //Vector local_error = r->local_error;
    //Vector lte = r->lte;

    Vector Gamma = r->Gamma;
    Vector G = r->G;
    Vector Psi = r->Psi;

    Matrix *Ydata = r->Ydata;  // Ydata[i] => Ystage
    Matrix *Udata = r->Udata;
    Vector *Zstage = r->Kstage;
    Matrix Ystage;
    double ti;

    // FIX: preallocate the variables
    Vector y_old = VectorNew(n_states+1);
    Vector *Z_old = NULL;
    Z_old = VectorArrayNew(r->rk_stages, n_states+1);
    Vector f = VectorNew(n_states+1);
    Vector res = VectorNew(rk_stages*(n_states+1));
    Vector dZ = VectorNew(rk_stages*(n_states+1));
    Vector y = VectorNew(n_states+1);
    Matrix L;// = MatrixNew(rk_stages*(n_states+1), rk_stages*(n_states+1));
    Matrix U;// = MatrixNew(rk_stages*(n_states+1), rk_stages*(n_states+1));
    int *P;// = (int *)malloc(rk_stages*(n_states+1)*sizeof(int));
    
    if (only_compute_lte == NO) {
        if (m > 0) {
            VectorSetAllTo(nlp_h, 0.0);
        }
        if (p > 0) {
            VectorSetAllTo(nlp_g, 0.0);
        }
    }
    for (i = 0; i < n_parameters; i++) {
        pc->e[i] = nlp_x->e[i];
    }
    for (j = 0; j < n_states; j++) {
        Yw->e[0][j] = nlp_x->e[n_parameters+j];
        y->e[j] = Yw->e[0][j];
    }
    Yw->e[0][n_states] = 0;
    y->e[n_states] = 0;

    if (n_controls > 0) {
        _setUdata(r, nlp_x);
        for (i = 0; i < n_controls; i++) {
            uc->e[i] = nlp_x->e[n_parameters+n_states+i];
        }
    }

    if (only_compute_lte == NO) {
        // form Gamma(y(t_initial), p)
        if (n_initial > 0) {
            r->ocp->initial_constraints(y, pc, Gamma);
            for (i = 0; i < n_initial; i++) {
                nlp_h->e[i] = Gamma->e[i];
            }
        }
        // form g(y,u,p,t_initial)
        if (n_inequality > 0) {
            r->ocp->inequality_constraints(y, uc, pc, T->e[0], G);
            for (i = 0; i < n_inequality; i++) {
                nlp_g->e[i] = G->e[i];
            }
        }
    }
    
    // integrate the ODEs using an implicit Runge-Kutta method

    int MAX_NEWTON_ITER = 15;
    int iter;
    int err;
    
    for (i_node = 0; i_node < (n_nodes-1); i_node++) {
        L = r->Lnode[i_node];
        U = r->Unode[i_node];
        P = r->pLUnode[i_node];
        // form H = (1/h) W (x) I - I (x) J
        // get the L, U, p factors of H
        double h = T->e[i_node+1] - T->e[i_node];
        err = _formIRKmatrixZ(r, y, uc, pc, T->e[i_node], h, L, U, P);
        if (err != 0) {
printf("IRK factor error: t = %e\n", T->e[i_node]);
//pause();
            r->integration_fatal_error = YES;
            if (m > 0) VectorSetAllTo(nlp_h, 0);
            if (p > 0) VectorSetAllTo(nlp_g, 0);
            return 0;
        }
        // initialize Kstage
        // FIX: extrapolate previous result
        if (i_node == 0 ) {
            for (i = 0; i < rk_stages; i++) {
                for (j = 0; j < n_states+1; j++) {
                    Zstage[i]->e[j] = 0;
                }
            }
        } else {
            // given y0, Z0, y1, tau = h1/h0, estimate Z1
            // using a Lagrange interpolation formula
            double tau;
            tau = h / (T->e[i_node] - T->e[i_node-1]);
            _estimateZ(r, y_old, Z_old, y, tau, Zstage);
        }
        Ystage = Ydata[i_node];
        
        // Y_i = y + Z_i
        // Ydot_i = (1/h) sum_{j = 1, s} w_{i,j} * Z_j
        // set res_i = -Ydot_i
        for (i = 0; i < rk_stages; i++) {
            for (j = 0; j < n_states+1; j++) {
                Ystage->e[i][j] = y->e[j] + Zstage[i]->e[j];
            }
            for (k = 0; k < n_states+1; k++) {
                //Ydot[i]->e[k] = 0.0;
                res->e[i*(n_states+1) + k] = 0;
            }
            for (j = 0; j < rk_stages; j++) {
                for (k = 0; k < n_states+1; k++) {
                    //Ydot[i]->e[k] += W[i][j] * Zstage[j]->e[k] / h;
                    res->e[i*(n_states+1) + k] -= W[i][j] * Zstage[j]->e[k] / h;
                }
            }
        }

        for (iter = 0; iter < MAX_NEWTON_ITER; iter++) {
            for (i = 0; i < rk_stages; i++) {
                ti = T->e[i_node] + h*rk_c[i];
                // get uc = u(t_i)
                for (j = 0; j < n_controls; j++) {
                    uc->e[j] = Udata[i_node]->e[i][j];
                }
                // get yc = Y_i
                for (j = 0; j < n_states+1; j++) {
                    yc->e[j] = Ystage->e[i][j];
                }
                // form f(Y_i, u_i, p, t_i)
                f->e[n_states] = r->ocp->differential_equations(yc, uc, pc, ti, f);
                // form res_i = Ydot_i - f(Y_i, u_i, p, t_i)
                for (j = 0; j < n_states+1; j++) {
                    res->e[i*(n_states+1) + j] += f->e[j];
                }
            }
            // solve H*dZ = -res = -(ydot - f(y,u,p,t))
            err = LUSolve(L, U, P, res, dZ);
            // Z = Z + dZ
            for (i = 0; i < rk_stages; i++) {
                for (j = 0; j < n_states+1; j++) {
                    Zstage[i]->e[j] += dZ->e[i*(n_states+1) + j];
                }
            }
            // Y_i = y + Z_i
            // Ydot_i = (1/h) sum_{j = 1, s} w_{i,j} * Z_j
            // set res_i = -Ydot_i
            for (i = 0; i < rk_stages; i++) {
                for (j = 0; j < n_states+1; j++) {
                    Ystage->e[i][j] = y->e[j] + Zstage[i]->e[j];
                }
                for (k = 0; k < n_states+1; k++) {
                    //Ydot[i]->e[k] = 0.0;
                     res->e[i*(n_states+1) + k] = 0;
                }
                for (j = 0; j < rk_stages; j++) {
                    for (k = 0; k < n_states+1; k++) {
                        //Ydot[i]->e[k] += W[i][j] * Zstage[j]->e[k] / h;
                        res->e[i*(n_states+1) + k] -= W[i][j] * Zstage[j]->e[k] / h;
                    }
                }
            }
            // test for convergence
            double norm_dZ = VectorNorm(dZ);
            if (norm_dZ <= 0.01*r->tolerance) {
                break;
            }
        }
//printf("i_node: %d, iter: %d\n", i_node, iter);
        // y_old = y
        // Z_old = Z
        for (j = 0; j < n_states+1; j++) {
            for (i = 0; i < rk_stages; i++) {
                Z_old[i]->e[j] = Zstage[i]->e[j];
            }
            y_old->e[j] = y->e[j];
        }
        
        // Assumes R-K methods with c[s] = 1
        // form Yw[i_node+1] = Ystage_s
        // update y = Yw[i_node+1]
        for (j = 0; j < n_states+1; j++) {
            Yw->e[i_node+1][j] = Ystage->e[rk_stages-1][j];
            y->e[j] = Yw->e[i_node+1][j];
        }
        
        // update u = u[i_node+1];
        int i_node1 = i_node+1;
        if ((r->c_type == CONSTANT) && (i_node1 == (r->n_nodes-1))) {
            i_node1 -= 1;
        }
        
        for (i = 0; i < n_controls; i++) {
            uc->e[i] = nlp_x->e[n_parameters + n_states + i_node1*n_controls + i];
        }
        
        // compute g(y(t+h),u(t+h),p,t+h)
        if (only_compute_lte == NO) {
            if (n_inequality > 0) {
                r->ocp->inequality_constraints(y, uc, pc, T->e[i_node+1], G);
                for (i = 0; i < n_inequality; i++) {
                    nlp_g->e[i + (i_node+1) * n_inequality] = G->e[i];
                }
            }
        }
        // FIX:  compute LTE 
    }

    if (only_compute_lte == YES) {
        return 0;
    }

    // form phi, Psi
    phi = r->ocp->terminal_constraints(yc, pc, Psi);
    for (i = 0; i < n_terminal; i++) {
        nlp_h->e[i + n_initial] = Psi->e[i];
    }

    VectorDelete(y_old);
    VectorArrayDelete(Z_old, r->rk_stages);
    VectorDelete(f);
    VectorDelete(res);
    VectorDelete(y);
    VectorDelete(dZ);
    //MatrixDelete(L);
    //MatrixDelete(U);
    //free(P);
//pause();
    return phi + Yw->e[n_nodes-1][n_states];
}

void _RKOCPDfhg_implicitZ(Vector nlp_x, Vector nlp_df, Matrix nlp_dh, Matrix nlp_dg) {
/*
 TODO: Instead of copying data to ys and us
 set ys->r = n_states+1, ys->e = state data
 set us->r = n_controls, us->e = control data
*/
    // nlp_x: vec[P, Y(0), U(0, :), U(1, :), ..., U(n_nodes-1, :)]
    // Yx: state+cost sensitivity with respect to nlp_x
    RKOCP r = thisRKOCP;
    int i_node, i, j;
    double t, ti, h;
    int m = r->nlp_m;
    int p = r->nlp_p;
    int n_nodes = r->n_nodes;
    int n_parameters = r->n_parameters;
    int n_states = r->n_states;
    int n_controls = r->n_controls;
    int n_inequality = r->ocp->n_inequality;
    int n_initial = r->ocp->n_initial;
    int n_terminal = r->ocp->n_terminal;

    int rk_stages = r->rk_stages;
    double *rk_c = r->rk_c->e;

    Vector T = r->T;
    Matrix Yw = r->Yw;
    Matrix Yx = r->Yx;
    Matrix *Ydata = r->Ydata;
    Matrix *Udata = r->Udata;
    Matrix *Yx_stage = r->Yx_stage;
    Matrix *Zx_stage = r->Kx_stage;
    Vector ys = r->yc;
    Vector us = r->uc;
    Vector ps = r->pc;

    // FIX: preallocate the following data
    Matrix Yx_dot = MatrixNew(n_states+1, r->nlp_n);
    Matrix res = MatrixNew(rk_stages*(n_states+1), r->nlp_n);
    Matrix dZx = MatrixNew(rk_stages*(n_states+1), r->nlp_n);
    Matrix Kx = MatrixNew(n_states+1, r->nlp_n);
    
    VectorSetAllTo(nlp_df, 0.0);

    if (m > 0) {
        MatrixSetAllTo(nlp_dh, 0.0);
    }

    if (p > 0) {
        MatrixSetAllTo(nlp_dg, 0.0);
    }

    if (r->integration_fatal_error == YES) {
        return;
    }
    
    if (n_controls > 0) {
        us = VectorNew(n_controls);
    }

    // parameters: ps
    if (n_parameters > 0) {
        ps = VectorNew(n_parameters);
        for (i = 0; i < n_parameters; i++) 
            ps->e[i] = nlp_x->e[i];
    }

    // sensitivity: Yx(t_initial)
    MatrixSetAllTo(Yx, 0.0);
    for (i = 0; i < n_states; i++) {
        Yx->e[i][i+n_parameters] = 1.0;
    }

    for (i_node = 0; i_node < n_nodes-1; i_node++) {
        ti = T->e[i_node];
        h = T->e[i_node+1] - ti;
        t = ti;

        // states: ys(t)
        for (i = 0; i < n_states+1; i++) 
            ys->e[i] = Yw->e[i_node][i];

        // controls: us(t)
        for (i = 0; i < n_controls; i++) 
            us->e[i] = nlp_x->e[n_parameters + n_states + i_node*n_controls + i];
        
        // dGamma -> nlp_dh
        if ((i_node == 0) && (n_initial > 0)) {
            _setDhGamma(r, ys, ps, nlp_dh);
        }

        // dg(t) -> nlp_dg
        if (n_inequality > 0) {
            _setDg(r, Yx, ys, us, ps, t, i_node, nlp_dg);
        }
        
        /*
        Integrate the sensitivity equations using an implicit Runge-Kuta method
         
        Yx_dot = Fy*Yx + Fu*ux + Fp*px
        
        Fy = [fy; Ly], Fu = [fu; Lu], Fp = [fp; Lp]
        Yx      : (n_states+1)-by-nlp_n matrix
        Yx_dot  : (n_states+1)-by-nlp_n matrix
        Ux      : n_controls-by-nlp_n matrix
        Ux      : n_parameters-by-nlp_n matrix
        
        Algorithm:
        
        recover the LUp factors of DPhi = (1/h) w (x) I - I (x) J
        
        estimate Zx[i], i = 1, 2, .. stages
        
        for iter = 0, 1, ..., MAX_NEWTON_ITER
            for i = 1,2,...,stages
                t_i = t + c_{i}*h
                Yx_i = Yx + Zx[i]
                Yx_dot_i = (1/h)*sum_{j=1,s} w_{i,j} * Zx[j]
                get Ux_i
                get Px_i
                compute Phi_i = Yx_dot_i - (Fy*Yx_i + Fu*Ux_i + Fp*Px_i)
            end
            
            solve DPhi*dZx = -Phi
            Zx += dZx
            Yx_i = Yx + Zx_i, i=1,2,...,s    
            check for convergence, i.e, ||dZx|| <= 0.01*tol
        end
        
        update Yx = Yx_s  (Assuming that the Runge-Kutta method has c[s] = 1
        
        Since these equations are linear we should expect convergence
        in a few iterations.
        */
        
        int iter, MAX_NEWTON = 15;
        // recover L, U, P
        Matrix L = r->Lnode[i_node];
        Matrix U = r->Unode[i_node];
        int *P = r->pLUnode[i_node];
        
        // set the column index for the sensitivity unknowns at this node
        int c0, c1;
        r->column_index = 0;
        switch (r->c_type) {
            case CONSTANT:
                r->column_index = n_parameters + n_states + (i_node + 1) * n_controls;
                break;
            case LINEAR:
                r->column_index = n_parameters + n_states + (i_node + 2) * n_controls;
                break;
            case CUBIC:
                c0 = (i_node + 6) * n_controls;
                c1 = n_nodes * n_controls;
                c0 = c0 < c1 ? c0 : c1;
                r->column_index = n_parameters + n_states + c0;
                break;
        }

        // get Zx_stage estimate
        _setZx_stage(r, i_node, Zx_stage);
        for (iter = 0; iter < MAX_NEWTON; iter++) {
            _setYx_stage_implicit(r, Yx, Zx_stage, Yx_stage);
            for (i = 0; i < rk_stages; i++) {
                ti = t + h*rk_c[i];
                _setYx_dot_stage_implicit(r, Zx_stage, h, i, Yx_dot);
                for (j = 0; j < n_controls; j++) us->e[j] = Udata[i_node]->e[i][j];
                for (j = 0; j < n_states+1; j++) ys->e[j] = Ydata[i_node]->e[i][j];
                _setResidual_implicit(r, i_node, i, Yx_dot, Yx_stage[i], ys, us, ps, t, Kx, res);
            }
            LUSolveM(L, U, P, res, dZx);
            // update Zx_stage
            _updateZx_stage(r, Zx_stage, dZx);
            // update Yx_stage
            _setYx_stage_implicit(r, Yx, Zx_stage, Yx_stage);
            // check for convergence
            double normdZx = MatrixNorm(dZx);
            if (normdZx <= 0.01*r->tolerance) {
                break;
            }
        }
        // set Yx = Yx_stage[rk_stage-1]
        for (i = 0; i < n_states+1; i++) {
            //for (j = 0; j < r->nlp_n; j++) { // FIX: make this i_node dependent
            for (j = 0; j < r->column_index; j++) {
                Yx->e[i][j] = Yx_stage[rk_stages-1]->e[i][j];
            }
        }
    }

    i_node = n_nodes-1;
    t = T->e[i_node];

    // states: ys(t)
    for (i = 0; i < n_states+1; i++) 
        ys->e[i] = Yw->e[i_node][i];

    // controls: us(t)
    for (i = 0; i < n_controls; i++) 
        us->e[i] = Udata[n_nodes-2]->e[rk_stages-1][i];
    
    // dg(t_{i_node}) -> nlp_dg
    if (n_inequality > 0) {
        _setDg(r, Yx, ys, us, ps, t, i_node, nlp_dg);
    }

    // Terminal conditions Psi and penalty phi
    if (r->ocp->Dterminal_constraints == NULL) {
        OCPFDphiPsi(r->ocp, ys, ps, r->phiy, r->phip, r->Psiy, r->Psip);
    } else {
        r->ocp->Dterminal_constraints(ys, ps, r->phiy, r->phip, r->Psiy, r->Psip);
    }
    // dPsi -> nlp_dh
    if (n_terminal > 0) {
        _setDhPsi(r, Yx, nlp_dh);
    }

    // df = integral of L
    for (i = 0; i < r->nlp_n; i++) 
        nlp_df->e[i] = Yx->e[n_states][i];

    // df += dphi
    _addDphi(r, Yx, nlp_df);

    MatrixDelete(Yx_dot);
    MatrixDelete(res);
    MatrixDelete(dZx);
    MatrixDelete(Kx);  
}

double _RKOCPfhgE(Vector nlp_x, Vector nlp_h, Vector nlp_g) {
    RKOCP r = thisRKOCP;
    double phi = 0;
    double t;
    int i_node, l;
    int i, j;
    int m = r->nlp_m;
    int p = r->nlp_p;
    int n_nodes = r->n_nodes;
    int n_parameters = r->n_parameters;
    int n_states = r->n_states;
    int n_controls = r->n_controls;
    int n_inequality = r->ocp->n_inequality;
    int n_initial = r->ocp->n_initial;
    int n_terminal = r->ocp->n_terminal;

    int rk_stages = r->rk_stages;
    double **rk_a = r->rk_a->e;
    double *rk_b = r->rk_b->e;
    double *rk_c = r->rk_c->e;
    double *rk_e = r->rk_e->e;

    Vector T = r->T;
    Matrix Yw = r->Yw;
    Vector yc = r->yc;
    Vector uc = r->uc;
    Vector pc = r->pc;
    Vector local_error = r->local_error;
    Vector lte = r->lte;

    Vector Gamma = r->Gamma;
    Vector G = r->G;
    Vector Psi = r->Psi;

    Matrix *Ydata = r->Ydata;
    Matrix *Udata = r->Udata;
    Vector *Kstage = r->Kstage;

    double _y;

    double fac0 = 0.2, fac1 = 5.0, beta = 0.9;
    int step, max_steps = 100;
    Vector Y0, Y1;
    double T0, h = 1.0e-8;
    Y0 = VectorNew(n_states+1);
    Y1 = VectorNew(n_states+1);

    
    if (m > 0) {
        VectorSetAllTo(nlp_h, 0.0);
    }
    
    if (p > 0) {
        VectorSetAllTo(nlp_g, 0.0);
    }
    
    for (i = 0; i < n_parameters; i++) {
        pc->e[i] = nlp_x->e[i];
    }
    
    for (j = 0; j < n_states; j++) {
        Yw->e[0][j] = nlp_x->e[n_parameters+j];
    }
    Yw->e[0][n_states] = 0;

    if (n_controls > 0) {
        _setUdata(r, nlp_x);
    }

    // integrate the ODEs using an explicit Runge-Kutta method

    Matrix Ystage;
    double delta, ti;

    for (i_node = 0; i_node < n_nodes-1; i_node++) {
        Ystage = Ydata[i_node]; // save the stage values for use in Dfhg

        ti = T->e[i_node];
        delta = T->e[i_node+1] - ti;
        t = ti;

        //- save the control values for use in Dfhg
        //- getControl(i_node, t, Udata[i_node][0]);
        for (i = 0; i < n_controls; i++) 
            uc->e[i] = Udata[i_node]->e[0][i];

        for (i = 0; i < n_states+1; i++) {
            _y = Yw->e[i_node][i];
            Ystage->e[0][i] = _y;
            yc->e[i] = _y;
            Y0->e[i] = _y;
        }

        // form Gamma(y(t_initial), p)
        if ((i_node == 0) && (n_initial > 0)) {
            r->ocp->initial_constraints(yc, pc, Gamma);
            for (i = 0; i < n_initial; i++) {
                nlp_h->e[i] = Gamma->e[i];
            }
        }
        // form g(y,u,p,t)
        if (n_inequality > 0) {
            r->ocp->inequality_constraints(yc, uc, pc, t, G);
            for (i = 0; i < n_inequality; i++) {
                nlp_g->e[i + i_node * n_inequality] = G->e[i];
            }
        }
        
        // K = F(y,u,p,t)
        Kstage[0]->e[n_states] = r->ocp->differential_equations(yc, uc, pc, t, Kstage[0]);
        
        T0 = t;
        if (i_node == 0)
            h = delta;
                
        for (step = 0; step <= max_steps; step++) {
            
            for (l = 1; l < rk_stages; l++) {
                t = T0 + rk_c[l]*h;
                //getControl(i_node, t, Udata[i_node][l]);
                //for (i = 0; i < n_controls; i++) 
                //    uc->e[i] = Udata[i_node]->e[0][i]; // constant control
    
                //Y_l = yc + h * sum_{j=0 to j=l-1} rk_a[l][j] * K_j
                for (i = 0; i < n_states + 1; i++) 
                    Ystage->e[l][i] = Y0->e[i];//Yw->e[i_node][i];
                for (j = 0; j < l; j++) {
                    for (i = 0; i < n_states+1; i++) {
                        Ystage->e[l][i] += h*rk_a[l][j]*Kstage[j]->e[i];
                    }
                }
    
                for (i = 0; i < n_states+1; i++) 
                    yc->e[i] = Ystage->e[l][i];
                    
                Kstage[l]->e[n_states] = r->ocp->differential_equations(yc, uc, pc, t, Kstage[l]);
            }
    
            // Y1 = Y0 + h * sum_{l=0 to stages-1} b[l] * K[l]
            for (i = 0; i < n_states+1; i++) {
                Y1->e[i] = Y0->e[i];
                local_error->e[i] = 0;
            }
            for (l = 0; l < rk_stages; l++) {
                for (i = 0; i < n_states+1; i++) {
                    Y1->e[i] += h*rk_b[l]*Kstage[l]->e[i];
                    local_error->e[i] += h*rk_e[l]*Kstage[l]->e[i];
                }
            }
            lte->e[i_node+1] = _compute_LTE(r, local_error->e, Y0->e, Y1->e);
            
            double sigma = pow(lte->e[i_node+1], -1.0/r->rk_order);
            double hnew = h * fmin(fac1, fmax(fac0, beta*sigma));
            //double _dt = fabs(T0 + h - T->e[i_node+1]);
            if (lte->e[i_node+1] < 1.0) {
                T0 = T0 + h;
                if (fabs(T0-T->e[i_node+1]) <= 100.0*DBL_EPSILON) {
                    for (i = 0; i < n_states+1; i++)
                        Yw->e[i_node+1][i] = Y1->e[i];
                    break;
                }
                for (i = 0; i < n_states+1; i++) {
                    Y0->e[i] = Y1->e[i];
                    yc->e[i] = Y1->e[i];
                }
                Kstage[0]->e[n_states] = r->ocp->differential_equations(yc, uc, pc, T0, Kstage[0]);
            }

            if (step >= max_steps) {
                for (i = 0; i < n_states+1; i++)
                        Yw->e[i_node+1][i] = Y1->e[i];
/*printf("local error:\n");
VectorPrint(local_error);
printf("Y1:\n");
VectorPrint(Y1);
printf("step: %d, lte: %e, T0: %e, T1: %e, h: %e, hnew: %e, sigma: %e, _dt: %e\n", step, lte->e[i_node+1], T0, T->e[i_node+1], h, hnew, sigma, _dt);
pause();*/
                break;
            }

/*printf("control:\n");
VectorPrint(uc);
printf("local error:\n");
VectorPrint(local_error);
printf("step: %d, lte: %e, T0: %e, T1: %e, h: %e, hnew: %e, sigma: %e, _dt: %e\n", step, lte->e[i_node+1], T0, T->e[i_node+1], h, hnew, sigma, _dt);
pause();*/
          
            h = hnew;

            if (T0+h > T->e[i_node+1])
                h = T->e[i_node+1] - T0;
//
        }
//pause();
    }

    // form g at time t_final
    t = T->e[n_nodes-1];
    if (n_controls > 0) {
        //getControl(n_nodes - 2, t, Udata[n_nodes - 1][0]);
        for (i = 0; i < n_controls; i++) 
            uc->e[i] = Udata[n_nodes-2]->e[rk_stages-1][i];
    }

    for (i = 0; i < n_states; i++) 
        yc->e[i] = Yw->e[n_nodes-1][i];
    
    if (n_inequality > 0) {
        r->ocp->inequality_constraints(yc, uc, pc, t, G);
        for (i = 0; i < n_inequality; i++) {
            nlp_g->e[i + (n_nodes-1)*n_inequality] = G->e[i];
        }
    }
    
    // form phi, Psi
    phi = r->ocp->terminal_constraints(yc, pc, Psi);
    for (i = 0; i < n_terminal; i++) {
        nlp_h->e[i + n_initial] = Psi->e[i];
    }

    VectorDelete(Y0);
    VectorDelete(Y1);
    
    return phi + Yw->e[n_nodes-1][n_states];
}
