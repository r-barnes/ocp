// Copyright (c) 2009,2010,2011,2012 Brian C. Fabien
// All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the license.txt file.

// TODO: remove remesh code

#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <float.h>

#include "runtime.h"
#include "vector.h"
#include "matrix.h"
#include "equations.h"
#include "nlpsolve.h"
#include "ocp.h"
#include "rkocp.h"

int RKOCPSolve(RKOCP r, Vector t0, Matrix y0, Matrix u0, Vector p0) {
    if (_RKOCPcheck(r, t0, y0, u0, p0) != 0) {
        return -1;
    }
    return _RKOCPsolveNLP(r);
}

void RKOCPDelete(RKOCP r) {
    if (r == NULL) {
        return;
    }
    if (r->rk_a != NULL) MatrixDelete(r->rk_a);
    if (r->rk_b != NULL) VectorDelete(r->rk_b);
    if (r->rk_c != NULL) VectorDelete(r->rk_c);
    if (r->rk_e != NULL) VectorDelete(r->rk_e);
    if (r->rk_w != NULL) MatrixDelete(r->rk_w);
    
    if (r->Yw != NULL) MatrixDelete(r->Yw);
    if (r->Ucubic != NULL) MatrixArrayDelete(r->Ucubic, r->n_nodes-1);
    
    if (r->x0 != NULL) VectorDelete(r->x0);
    if (r->lambda0 != NULL) VectorDelete(r->lambda0);
    if (r->mu0 != NULL) VectorDelete(r->mu0);
    
    if (r->yc != NULL) VectorDelete(r->yc);
    if (r->pc != NULL) VectorDelete(r->pc);
    if (r->uc != NULL) VectorDelete(r->uc);

    if (r->F != NULL) VectorDelete(r->F);
    if (r->Fy != NULL) MatrixDelete(r->Fy);
    if (r->Fu != NULL) MatrixDelete(r->Fu);
    if (r->Fp != NULL) MatrixDelete(r->Fp);
    
    if (r->Ly != NULL) VectorDelete(r->Ly);
    if (r->Lu != NULL) VectorDelete(r->Lu);
    if (r->Lp != NULL) VectorDelete(r->Lp);

    if (r->G != NULL) VectorDelete(r->G);
    if (r->Gy != NULL) MatrixDelete(r->Gy);
    if (r->Gu != NULL) MatrixDelete(r->Gu);
    if (r->Gp != NULL) MatrixDelete(r->Gp);    
    
    if (r->Gamma != NULL) VectorDelete(r->Gamma);
    if (r->Gammay != NULL) MatrixDelete(r->Gammay);
    if (r->Gammap != NULL) MatrixDelete(r->Gammap);

    if (r->phiy != NULL) VectorDelete(r->phiy);
    if (r->phip != NULL) VectorDelete(r->phip);
    
    if (r->Psi != NULL) VectorDelete(r->Psi);
    if (r->Psiy != NULL) MatrixDelete(r->Psiy);
    if (r->Psip != NULL) MatrixDelete(r->Psip);
    
    if (r->Ydata != NULL) MatrixArrayDelete(r->Ydata, r->n_nodes-1);
    if (r->Kstage != NULL) VectorArrayDelete(r->Kstage, r->rk_stages);
    if (r->Udata != NULL) MatrixArrayDelete(r->Udata, r->n_nodes-1);
    if (r->Yx != NULL) MatrixDelete(r->Yx);
    if (r->Yx_stage != NULL) MatrixArrayDelete(r->Yx_stage, r->rk_stages);
    if (r->Kx_stage != NULL) MatrixArrayDelete(r->Kx_stage, r->rk_stages);
    
    if (r->lte != NULL) VectorDelete(r->lte);
    if (r->local_error != NULL) VectorDelete(r->local_error);
    
    if (r->Lnode != NULL) MatrixArrayDelete(r->Lnode, r->n_nodes);
    if (r->Unode != NULL) MatrixArrayDelete(r->Unode, r->n_nodes);
    if (r->pLUnode != NULL) {
        int i;
        for (i = 0; i < r->n_nodes; i++) {
            free(r->pLUnode[i]);
        }
        free(r->pLUnode);
    }    
    
    free(r);
}

int _RKOCPcheck(RKOCP r, Vector t0, Matrix y0, Matrix u0, Vector p0) {
    //t0: n_nodes
    //y0: n_nodes, n_states
    //u0: n_nodes, n_controls
    //p0: n_parameters
    // set U = u0, Y = y0, T = y0, P = p0
    int ns = 0, nc = 0, np = 0, nn = 0;
    int i, j;

    // check the OCP
    int err = OCPCheck(r->ocp);
    if (err != 0) {
        RuntimeWarning("_RKOCPcheck: invalid OCP specified");
        return -1;
    }
    // check the inputs to RKOCP
    if ((t0 == NULL) || (y0 == NULL)) {
        RuntimeWarning("_RKOCPcheck: NULL pointer for t0 or y0");
        return -1;
    }

    nn = t0->r;
    if (nn < 2) {
        RuntimeWarning("_RKOCPcheck: t0->r < 2");
        return -1;
    }

    // check that the time nodes are in the correct sequence, and > epsilon
    double dt0 = t0->e[1] - t0->e[0];
    for (i = 0; i < nn-1; i++) {
        double dt = t0->e[i+1] - t0->e[i];
        if (fabs(dt) <= DBL_EPSILON) {
            RuntimeWarning("_RKOCPcheck: time node step size too small");
            return -1;
        }
        if (dt*dt0 < 0) {
            RuntimeWarning("_RKOCPcheck: time nodes not in sequence");
            return -1;
        }
        dt0 = dt;
    }

    ns = y0->c;

    if (ns != r->ocp->n_states) {
        RuntimeWarning("_RKOCPcheck: y0 does not have the correct number of columns");
        return -1;
    }

    if (nn != y0->r) {
        RuntimeWarning("_RKOCPcheck: dimensions of t0 and y0 are incompatible");
        return -1;
    }

    // check control type

    r->c_type = r->control_type;
    switch (r->c_type) {
        case CONSTANT:
            r->controlTypeName = "CONSTANT";
            break;
        case LINEAR:
            r->controlTypeName = "LINEAR";
            break;
        case CUBIC:
            r->controlTypeName = "CUBIC";
            if (nn < 6) {
            	RuntimeWarning("_RKOCPcheck: CUBIC control type requires at least 6 nodes");
            	return -1;
            }
            break;
        default:
            RuntimeWarning("_RKOCPcheck: invalid control_type specified.");
            return -1;
    }

    if ((r->ocp->n_controls > 0) && (u0 == NULL)) {
        RuntimeWarning("_RKOCPcheck: ocp.n_controls > 0 but u0 == NULL");
        return -1;
    }

    if (u0 != NULL) {
        nc = u0->c;
        if (nc != r->ocp->n_controls) {
            RuntimeWarning("_RKOCPcheck: u0 does not have the correct number of rows");
            return -1;
        }
        if (nn != u0->r) {
            RuntimeWarning("_RKOCPcheck: u0 does not have the correct number of columns");
            return -1;
        }
    }

    // parameters
    if ((r->ocp->n_parameters > 0) && (p0 == NULL)) {
        RuntimeWarning("_RKOCPcheck: ocp.n_parameters > 0 but p0 == NULL");
        return -1;
    }

    if (p0 != NULL) {
        np = p0->r;
        if (np != r->ocp->n_parameters) {
            RuntimeWarning("_RKOCPcheck: p0 does not have the correct number of elements (n_parameters)");
            return -1;
        }
    }

    // check tolerance
    if (r->tolerance < (100.0 * DBL_EPSILON)) {
        RuntimeWarning("_RKOCPcheck: convergence tolerance set too small.\nUsing tolerance = 1.0e-6");
        r->tol = 1.0e-6;
    } else {
        r->tol = r->tolerance;
    }

    // check max iterations
    if (r->maximum_iterations < 1) {
        RuntimeWarning("_RKOCPcheck: maximum_iterations < 1\n.Using maximum_iterations = 5000");
        r->max_iter = 5000;
    } else {
        r->max_iter = r->maximum_iterations;
    }

    // check max nodes
/*
    if (r->maximum_nodes < nn) {
        RuntimeWarning("_RKOCPcheck: maximum_nodes set too small.\nUsing maximum_nodes = n_nodes");
        r->max_nodes = nn;
    } else {
        r->max_nodes = r->maximum_nodes;
    }
*/
    // check RK method
/*
    if ((r->runge_kutta_method < 0) || (r->runge_kutta_method > 10)) {
        RuntimeWarning("_RKOCPcheck: invalid runge_kutta_method specified.\nUsing runge_kutta_method = DormandPrince8_7");
        r->rk_method = DormandPrince8_7;
    } else {
        r->rk_method = r->runge_kutta_method;
    }
*/
    r->rk_method = r->runge_kutta_method;
    switch (r->rk_method) {
        case BogackiShampine3_2:
            _getRKcoefficients0(r);
            break;
        case Classical4_3:
            _getRKcoefficients1(r);
            break;
        case DormandPrince5_4:
            _getRKcoefficients2(r);
            break;
        case Verner6_5:
            _getRKcoefficients3(r);
            break;
        case Verner7_6:
            _getRKcoefficients4(r);
            break;
        case Verner8_7:
            _getRKcoefficients5(r);
            break;
        case Verner9_8:
            _getRKcoefficients6(r);
            break;
        case DormandPrince8_7:
            _getRKcoefficients7(r);
            break;
        case DVERK6_5:
            _getRKcoefficients8(r);
            break;
        case RungeKuttaFehlberg7_8:
            _getRKcoefficients9(r);
            break;
        case Merson4_5:
            _getRKcoefficients10(r);
            break;
        case RadauIIA3:
            _getRKcoeffRadauIIA3(r);
            break;
        case RadauIIA5:
            _getRKcoeffRadauIIA5(r);
            break;
        default:
            RuntimeWarning("_RKOCPcheck: invalid runge_kutta_method specified.");
            return -1;
    }

    r->T = t0;
    r->Y = y0;
    r->Yw = MatrixNew(nn, ns+1);
    for (i = 0; i < nn; i++) {
        for (j = 0; j < ns; j++) {
            r->Yw->e[i][j] = y0->e[i][j];
        }
    }
    if (nc > 0) {
        r->U = u0;
    }
    if (np > 0) {
        r->P = p0;
    }
    r->n_states = ns;
    r->n_controls = nc;
    r->n_parameters = np;
    r->n_nodes = nn;

    if (r->c_type == CUBIC) {
        r->Ucubic = MatrixArrayNew(r->n_nodes-1, r->rk_stages, 6);
        _RKOCPformUcubic(r);
    }
    return 0;
}

double ***_RKOCPa3d_allo_dbl(int dim1, int dim2, int dim3) {
    double ***array;
    int i, j;

    if (dim1 < 1 || dim2 < 1 || dim3 < 1) {
        array = NULL;
        return (array);
    }

    array = (double ***) malloc(dim1 * sizeof (double **));

    if (array == NULL) {
        return (array);
    }
    for (i = 0; i < dim1; i++) {
        array[i] = _RKOCPa2d_allo_dbl(dim2, dim3);
        if (array[i] == NULL) {
            for (j = i - 1; j >= 0; j--) {
                _RKOCPa2d_free_dbl(array[j], dim2);
            }
            free(array);
            array = NULL;
            return (array);
        }
    }
    return (array);
}

void _RKOCPa3d_free_dbl(double ***array, int dim1, int dim2) {
    int i;

    if (array == NULL) {
        return;
    }
    for (i = dim1 - 1; i >= 0; i--) {
        if (array[i] != NULL) {
            _RKOCPa2d_free_dbl(array[i], dim2);
        }
    }
    if (array != NULL) {
        free(array);
    }
}

double **_RKOCPa2d_allo_dbl(int dim1, int dim2) {
    double **array;
    int i, j;

    if (dim1 < 1 || dim2 < 1) {
        array = NULL;
        return (array);
    }

    array = (double **) malloc(dim1 * sizeof (double *));

    if (array == NULL) {
        return ( array);
    }
    for (i = 0; i < dim1; i++) {
        array[i] = (double *) malloc(dim2 * sizeof (double));
        if (array[i] == NULL) {
            for (j = i-1; j >= 0; j--)
                free(array[j]);
            free(array);
            array = NULL;
            return (array);
        }
    }
    return (array);
}

void _RKOCPa2d_free_dbl(double **array, int dim1) {
    int i;

    if (array == NULL) {
        return;
    }
    for (i = dim1-1; i >= 0; i--) {
        if (array[i] != NULL) {
            free(array[i]);
        }
    }
    if (array != NULL) {
        free(array);
    }
}

void _RKOCPformUcubic(RKOCP r) {
    int n_nodes = r->n_nodes;
    int rk_stages = r->rk_stages;
    Matrix *Ucubic = r->Ucubic;
    double **h = _RKOCPa2d_allo_dbl(r->rk_stages, 4);
    double *T = r->T->e;
    int j, k;

    for (k = 0; k < r->rk_stages; k++) {
        _RKOCPHermite(h[k], r->rk_c->e[k]);
    }
    double t0, t1, t2, t3, t4, t5;
    double h00, h10, h01, h11;
    for (j = 0; j < (n_nodes - 1); j++) {
        if ((j > 1) && (j < n_nodes - 3)) {
            t0 = T[j - 2];
            t1 = T[j - 1];
            t2 = T[j];
            t3 = T[j + 1];
            t4 = T[j + 2];
            t5 = T[j + 3];
            for (k = 0; k < rk_stages; k++) {
                h00 = h[k][0];
                h10 = h[k][1];
                h01 = h[k][2];
                h11 = h[k][3];
                Ucubic[j]->e[k][0] = h10 * (t2 - t1) * (t2 - t3) * (t2 - t4) * (t3 - t2) / ((t0 - t1) * (t0 - t2) * (t0 - t3) * (t0 - t4));
                Ucubic[j]->e[k][1] = h11 * pow((t3 - t2), 2) * (t3 - t4) * (t3 - t5) / ((t1 - t2) * (t1 - t3) * (t1 - t4) * (t1 - t5)) + h10 * (t2 - t0) * (t2 - t3) * (t2 - t4) * (t3 - t2) / ((t1 - t0) * (t1 - t2) * (t1 - t3) * (t1 - t4));
                Ucubic[j]->e[k][2] = h00 + h10 * (t3 - t2) * (1.0 / (t2 - t0) + 1.0 / (t2 - t1) + 1.0 / (t2 - t3) + 1.0 / (t2 - t4)) + h11 * (t3 - t1) * (t3 - t2) * (t3 - t4) * (t3 - t5) / ((t2 - t1) * (t2 - t3) * (t2 - t4) * (t2 - t5));
                Ucubic[j]->e[k][3] = h01 + h11 * (t3 - t2) * (1.0 / (t3 - t1) + 1.0 / (t3 - t2) + 1.0 / (t3 - t4) + 1.0 / (t3 - t5)) + h10 * (t2 - t0) * (t2 - t1) * (t2 - t4) / ((t3 - t0) * (t3 - t1) * (t3 - t4));
                Ucubic[j]->e[k][4] = h11 * pow((t3 - t2), 2) * (t3 - t1) * (t3 - t5) / ((t4 - t1) * (t4 - t2) * (t4 - t3) * (t4 - t5)) + h10 * (t2 - t0) * (t2 - t1) * (t2 - t3) * (t3 - t2) / ((t4 - t0) * (t4 - t1) * (t4 - t2) * (t4 - t3));
                Ucubic[j]->e[k][5] = h11 * pow((t3 - t2), 2) * (t3 - t1) * (t3 - t4) / ((t5 - t1) * (t5 - t2) * (t5 - t3) * (t5 - t4));
            }
        } else if (j == 0) {
            t0 = T[j];
            t1 = T[j + 1];
            t2 = T[j + 2];
            t3 = T[j + 3];
            t4 = T[j + 4];
            for (k = 0; k < rk_stages; k++) {
                h00 = h[k][0];
                h10 = h[k][1];
                h01 = h[k][2];
                h11 = h[k][3];
                Ucubic[j]->e[k][0] = h00 + h10 * (t1 - t0) * (1.0 / (t0 - t1) + 1.0 / (t0 - t2) + 1.0 / (t0 - t3) + 1.0 / (t0 - t4)) + h11 * (t1 - t0) * (t1 - t2) * (t1 - t3) * (t1 - t4) / ((t0 - t1) * (t0 - t2) * (t0 - t3) * (t0 - t4));
                Ucubic[j]->e[k][1] = h01 + h11 * (t1 - t0) * (1.0 / (t1 - t0) + 1.0 / (t1 - t2) + 1.0 / (t1 - t3) + 1.0 / (t1 - t4)) + h10 * (t0 - t2) * (t0 - t3) * (t0 - t4) / ((t1 - t2) * (t1 - t3) * (t1 - t4));
                Ucubic[j]->e[k][2] = h11 * pow((t1 - t0), 2) * (t1 - t3) * (t1 - t4) / ((t2 - t0) * (t2 - t1) * (t2 - t3) * (t2 - t4)) + h10 * (t0 - t1) * (t0 - t3) * (t0 - t4) * (t1 - t0) / ((t2 - t0) * (t2 - t1) * (t2 - t3) * (t2 - t4));
                Ucubic[j]->e[k][3] = h11 * pow((t1 - t0), 2) * (t1 - t2) * (t1 - t4) / ((t3 - t0) * (t3 - t1) * (t3 - t2) * (t3 - t4)) + h10 * (t0 - t1) * (t0 - t2) * (t0 - t4) * (t1 - t0) / ((t3 - t0) * (t3 - t1) * (t3 - t2) * (t3 - t4));
                Ucubic[j]->e[k][4] = h11 * pow((t1 - t0), 2) * (t1 - t2) * (t1 - t3) / ((t4 - t0) * (t4 - t1) * (t4 - t2) * (t4 - t3)) + h10 * (t0 - t1) * (t0 - t2) * (t0 - t3) * (t1 - t0) / ((t4 - t0) * (t4 - t1) * (t4 - t2) * (t4 - t3));
                Ucubic[j]->e[k][5] = 0.0;
            }
        } else if (j == 1) {
            t0 = T[j - 1];
            t1 = T[j];
            t2 = T[j + 1];
            t3 = T[j + 2];
            t4 = T[j + 3];
            for (k = 0; k < rk_stages; k++) {
                h00 = h[k][0];
                h10 = h[k][1];
                h01 = h[k][2];
                h11 = h[k][3];
                Ucubic[j]->e[k][0] = h11 * pow((t2 - t1), 2) * (t2 - t3) * (t2 - t4) / ((t0 - t1) * (t0 - t2) * (t0 - t3) * (t0 - t4)) + h10 * (t1 - t2) * (t1 - t3) * (t1 - t4) * (t2 - t1) / ((t0 - t1) * (t0 - t2) * (t0 - t3) * (t0 - t4));
                Ucubic[j]->e[k][1] = h00 + h10 * (t2 - t1) * (1.0 / (t1 - t0) + 1.0 / (t1 - t2) + 1.0 / (t1 - t3) + 1.0 / (t1 - t4)) + h11 * (t2 - t0) * (t2 - t1) * (t2 - t3) * (t2 - t4) / ((t1 - t0) * (t1 - t2) * (t1 - t3) * (t1 - t4));
                Ucubic[j]->e[k][2] = h01 + h11 * (t2 - t1) * (1.0 / (t2 - t0) + 1.0 / (t2 - t1) + 1.0 / (t2 - t3) + 1.0 / (t2 - t4)) + h10 * (t1 - t0) * (t1 - t3) * (t1 - t4) / ((t2 - t0) * (t2 - t3) * (t2 - t4));
                Ucubic[j]->e[k][3] = h11 * pow((t2 - t1), 2) * (t2 - t0) * (t2 - t4) / ((t3 - t0) * (t3 - t1) * (t3 - t2) * (t3 - t4)) + h10 * (t1 - t0) * (t1 - t2) * (t1 - t4) * (t2 - t1) / ((t3 - t0) * (t3 - t1) * (t3 - t2) * (t3 - t4));
                Ucubic[j]->e[k][4] = h11 * pow((t2 - t1), 2) * (t2 - t0) * (t2 - t3) / ((t4 - t0) * (t4 - t1) * (t4 - t2) * (t4 - t3)) + h10 * (t1 - t0) * (t1 - t2) * (t1 - t3) * (t2 - t1) / ((t4 - t0) * (t4 - t1) * (t4 - t2) * (t4 - t3));
                Ucubic[j]->e[k][5] = 0.0;
            }
        } else if (j == n_nodes - 3) {
            t0 = T[j - 2];
            t1 = T[j - 1];
            t2 = T[j];
            t3 = T[j + 1];
            t4 = T[j + 2];
            for (k = 0; k < rk_stages; k++) {
                h00 = h[k][0];
                h10 = h[k][1];
                h01 = h[k][2];
                h11 = h[k][3];
                Ucubic[j]->e[k][0] = 0.0;
                Ucubic[j]->e[k][1] = h11 * pow((t3 - t2), 2) * (t3 - t1) * (t3 - t4) / ((t0 - t1) * (t0 - t2) * (t0 - t3) * (t0 - t4)) + h10 * (t2 - t1) * (t2 - t3) * (t2 - t4) * (t3 - t2) / ((t0 - t1) * (t0 - t2) * (t0 - t3) * (t0 - t4));
                Ucubic[j]->e[k][2] = h11 * pow((t3 - t2), 2) * (t3 - t0) * (t3 - t4) / ((t1 - t0) * (t1 - t2) * (t1 - t3) * (t1 - t4)) + h10 * (t2 - t0) * (t2 - t3) * (t2 - t4) * (t3 - t2) / ((t1 - t0) * (t1 - t2) * (t1 - t3) * (t1 - t4));
                Ucubic[j]->e[k][3] = h00 + h10 * (t3 - t2) * (1.0 / (t2 - t0) + 1.0 / (t2 - t1) + 1.0 / (t2 - t3) + 1.0 / (t2 - t4)) + h11 * (t3 - t0) * (t3 - t1) * (t3 - t2) * (t3 - t4) / ((t2 - t0) * (t2 - t1) * (t2 - t3) * (t2 - t4));
                Ucubic[j]->e[k][4] = h01 + h11 * (t3 - t2) * (1.0 / (t3 - t0) + 1.0 / (t3 - t1) + 1.0 / (t3 - t2) + 1.0 / (t3 - t4)) + h10 * (t2 - t0) * (t2 - t1) * (t2 - t4) / ((t3 - t0) * (t3 - t1) * (t3 - t4));
                Ucubic[j]->e[k][5] = h11 * pow((t3 - t2), 2) * (t3 - t0) * (t3 - t1) / ((t4 - t0) * (t4 - t1) * (t4 - t2) * (t4 - t3)) + h10 * (t2 - t0) * (t2 - t1) * (t2 - t3) * (t3 - t2) / ((t4 - t0) * (t4 - t1) * (t4 - t2) * (t4 - t3));
            }

        } else if (j == n_nodes - 2) {
            t0 = T[j - 3];
            t1 = T[j - 2];
            t2 = T[j - 1];
            t3 = T[j];
            t4 = T[j + 1];
            for (k = 0; k < rk_stages; k++) {
                h00 = h[k][0];
                h10 = h[k][1];
                h01 = h[k][2];
                h11 = h[k][3];
                Ucubic[j]->e[k][0] = 0.0;
                Ucubic[j]->e[k][1] = h11 * pow((t4 - t3), 2) * (t4 - t1) * (t4 - t2) / ((t0 - t1) * (t0 - t2) * (t0 - t3) * (t0 - t4)) + h10 * (t3 - t1) * (t3 - t2) * (t3 - t4) * (t4 - t3) / ((t0 - t1) * (t0 - t2) * (t0 - t3) * (t0 - t4));
                Ucubic[j]->e[k][2] = h11 * pow((t4 - t3), 2) * (t4 - t0) * (t4 - t2) / ((t1 - t0) * (t1 - t2) * (t1 - t3) * (t1 - t4)) + h10 * (t3 - t0) * (t3 - t2) * (t3 - t4) * (t4 - t3) / ((t1 - t0) * (t1 - t2) * (t1 - t3) * (t1 - t4));
                Ucubic[j]->e[k][3] = h11 * pow((t4 - t3), 2) * (t4 - t0) * (t4 - t1) / ((t2 - t0) * (t2 - t1) * (t2 - t3) * (t2 - t4)) + h10 * (t3 - t0) * (t3 - t1) * (t3 - t4) * (t4 - t3) / ((t2 - t0) * (t2 - t1) * (t2 - t3) * (t2 - t4));
                Ucubic[j]->e[k][4] = h00 + h10 * (t4 - t3) * (1.0 / (t3 - t0) + 1.0 / (t3 - t1) + 1.0 / (t3 - t2) + 1.0 / (t3 - t4)) + h11 * (t4 - t0) * (t4 - t1) * (t4 - t2) * (t4 - t3) / ((t3 - t0) * (t3 - t1) * (t3 - t2) * (t3 - t4));
                Ucubic[j]->e[k][5] = h01 + h11 * (t4 - t3) * (1.0 / (t4 - t0) + 1.0 / (t4 - t1) + 1.0 / (t4 - t2) + 1.0 / (t4 - t3)) + h10 * (t3 - t0) * (t3 - t1) * (t3 - t2) / ((t4 - t0) * (t4 - t1) * (t4 - t2));
            }
        }
    }
    _RKOCPa2d_free_dbl(h, r->rk_stages);
}

void _RKOCPHermite(double *L, double t) {
    double t2 = t * t;
    double t3 = t2 * t;
    L[0] = 2.0 * t3 - 3.0 * t2 + 1.0;
    L[1] = t3 - 2.0 * t2 + t;
    L[2] = -2.0 * t3 + 3.0 * t2;
    L[3] = t3 - t2;
}

int _RKOCPsolveNLP(RKOCP r) {
    int err = 0;
    NLPOptions opt;
    opt.method = r->nlp_method;
    opt.tolerance = r->tol;
    opt.maximumIterations = r->max_iter;
    opt.displayData = r->display_data;
    
    _formInitialEstimate(r);

    if (r->rk_method > 20) {
        //r->nlp_stats = NLPsolve(_RKOCPfhg_implicit, NULL, r->x0, r->lambda0, r->mu0, opt);
        //r->nlp_stats = NLPsolve(_RKOCPfhg_implicitZ, NULL, r->x0, r->lambda0, r->mu0, opt);
        r->nlp_stats = NLPsolve(_RKOCPfhg_implicitZ, _RKOCPDfhg_implicitZ, r->x0, r->lambda0, r->mu0, opt);
        //r->only_compute_lte = YES;
        //_RKOCPfhg(r->x0, NULL, NULL);
        //r->max_lte = VectorMax(r->lte);

        err = r->nlp_stats.exitCode;
    
        _getSolution(r);
        
        if (r->integration_fatal_error == YES) {
            VectorSetAllTo(r->lte, 1.0);
            r->max_lte = 1.0;
            return -99;
        }
    } else {
        r->nlp_stats = NLPsolve(_RKOCPfhg, _RKOCPDfhg, r->x0, r->lambda0, r->mu0, opt);
        
        err = r->nlp_stats.exitCode;
    
        _getSolution(r);
        
        if (r->integration_fatal_error == YES) {
            VectorSetAllTo(r->lte, 1.0);
            r->max_lte = 1.0;
            return -99;
        }
        
        r->only_compute_lte = YES;
        _RKOCPfhg(r->x0, NULL, NULL);
        r->max_lte = VectorMax(r->lte);
    }
    
    //r->nlp_stats = NLPsolve(_RKOCPfhg, NULL, r->x0, r->lambda0, r->mu0, opt);
    //r->nlp_stats = NLPsolve(_RKOCPfhgE, NULL, r->x0, r->lambda0, r->mu0, opt);
    //r->max_lte = VectorMax(r->lte);
    
    return err;
}

void _formInitialEstimate(RKOCP r) {
    // set nlp.n, nlp.m, nlp.p, nlp.tolerance, nlp.maximum_iterations
    int i, j, k0, m, p, n = 0;
    int n_nodes = r->n_nodes;
    int n_states = r->n_states;
    int n_controls = r->n_controls;
    int n_parameters = r->n_parameters;

    switch (r->c_type) {
        case CONSTANT:
            n = (n_nodes-1) * n_controls + n_parameters + n_states;
            break;
        case LINEAR:
            n = n_nodes * n_controls + n_parameters + n_states;
            break;
        case CUBIC:
            n = n_nodes * n_controls + n_parameters + n_states;
            break;
    }

    // set nlp.m (number of equality constrints in the NLP problem)
    m = r->ocp->n_initial + r->ocp->n_terminal;
    r->n_initial = r->ocp->n_initial;
    r->n_terminal = r->ocp->n_terminal;

    // set nlp.p (number of inequality constrints in the NLP problem)
    p = n_nodes * r->ocp->n_inequality;
    r->n_inequality = r->ocp->n_inequality;

    // set x0, sqp.x, sqp.lambda, sqp.mu
    r->x0 = VectorNew(n);
    if (m > 0) {
        r->lambda0 = VectorNew(m);
    }
    if (p > 0) {
        r->mu0 = VectorNew(p);
    }

    // x = [P; Yi; U[0, :]; U[1: ]; ... U[n_nodes-1, :]]

    if (n_parameters > 0) {
        for (i = 0; i < n_parameters; i++) {
            r->x0->e[i] = r->P->e[i];
        }
    }

    for (i = 0; i < n_states; i++) {
        r->x0->e[n_parameters + i] = r->Y->e[0][i]; // Y(n_nodes, n_states)
    }

    if (n_controls > 0) {
        k0 = n_parameters + n_states;
        switch (r->c_type) {
            case CONSTANT:
                for (i = 0; i < n_nodes - 1; i++) {
                    for (j = 0; j < n_controls; j++) {
                        r->x0->e[k0 + j] = r->U->e[i][j];
                    }
                    k0 += n_controls;
                }
                break;
            case LINEAR:
                for (i = 0; i < n_nodes; i++) {
                    for (j = 0; j < n_controls; j++) {
                        r->x0->e[k0 + j] = r->U->e[i][j];
                    }
                    k0 += n_controls;
                }
                break;
            case CUBIC:
                for (i = 0; i < n_nodes; i++) {
                    for (j = 0; j < n_controls; j++) {
                        r->x0->e[k0 + j] = r->U->e[i][j];
                    }
                    k0 += n_controls;
                }
                break;
        }
    }

    r->yc = VectorNew(n_states+1);
    if (n_parameters > 0) {
        r->pc = VectorNew(n_parameters);
    }
    if (n_controls > 0) {
        r->uc = VectorNew(n_controls);
    }

    r->F = VectorNew(n_states);
    r->Fy = MatrixNew(n_states, n_states);
    if (n_controls > 0) {
        r->Fu = MatrixNew(n_states, n_controls);
    }
    if (n_parameters > 0) {
        r->Fp = MatrixNew(n_states, n_parameters);
    }

    r->Ly = VectorNew(n_states);
    if (n_controls > 0) {
        r->Lu = VectorNew(n_controls);
    }
    if (n_parameters > 0) {
        r->Lp = VectorNew(n_parameters);
    }

    if (r->n_inequality > 0) {
        r->G = VectorNew(r->n_inequality);
        r->Gy = MatrixNew(r->n_inequality, n_states);
        if (n_controls > 0) {
            r->Gu = MatrixNew(r->n_inequality, n_controls);
        }
        if (n_parameters > 0) {
            r->Gp = MatrixNew(r->n_inequality, n_parameters);
        }
    }

    if (r->n_initial > 0) {
        r->Gamma = VectorNew(r->n_initial);
        r->Gammay = MatrixNew(r->n_initial, n_states);
        if (n_parameters > 0) {
            r->Gammap = MatrixNew(r->n_initial, n_parameters);
        }
    }

    r->phiy = VectorNew(n_states);
    if (n_parameters > 0) {
        r->phip = VectorNew(n_parameters);
    }

    if (r->n_terminal > 0) {
        r->Psi = VectorNew(r->n_terminal);
        r->Psiy = MatrixNew(r->n_terminal, n_states);
        if (n_parameters > 0) {
            r->Psip = MatrixNew(r->n_terminal, n_parameters);
        }
    }

    r->Ydata = MatrixArrayNew(n_nodes-1, r->rk_stages, n_states+1);
    r->Kstage = VectorArrayNew(r->rk_stages, n_states+1);
    if (n_controls > 0) {
        r->Udata = MatrixArrayNew(n_nodes-1, r->rk_stages, n_controls);
    }

    r->Yx = MatrixNew(n_states+1, n);
    r->Yx_stage = MatrixArrayNew(r->rk_stages, n_states+1, n);
    r->Kx_stage = MatrixArrayNew(r->rk_stages, n_states+1, n);

    r->lte = VectorNew(n_nodes);
    r->local_error = VectorNew(n_states+1);
    r->only_compute_lte = NO;

    r->nlp_n = n;
    r->nlp_m = m;
    r->nlp_p = p;
    
    // implicit Runge-Kutta data
    if (r->runge_kutta_method > 20) {
        int n_rk = r->rk_stages*(n_states+1);
        r->Lnode = MatrixArrayNew(n_nodes, n_rk, n_rk);
        r->Unode = MatrixArrayNew(n_nodes, n_rk, n_rk);
        r->pLUnode = (int **) malloc(n_nodes * sizeof (int *));
        for (i = 0; i < n_nodes; i++) {
            r->pLUnode[i] = (int *) malloc(n_rk * sizeof (int));
        }
    }
}

void _getSolution(RKOCP r) {
    // put nlp.x -> P, U; Yw ->Y
    int i, j;
    int rk_stages = r->rk_stages;
    int n_nodes = r->n_nodes;
    int n_states = r->n_states;
    int n_controls = r->n_controls;
    int n_parameters = r->n_parameters;

    if (n_parameters > 0) {
        for (i = 0; i < n_parameters; i++) {
            r->P->e[i] = r->x0->e[i];
        }
    }
    for (i = 0; i < n_nodes; i++) {
        for (j = 0; j < n_states; j++) {
            r->Y->e[i][j] = r->Yw->e[i][j];
        }
    }
    if (n_controls > 0) {
        for (i = 0; i < n_nodes - 1; i++) {
            for (j = 0; j < n_controls; j++) {
                r->U->e[i][j] = r->Udata[i]->e[0][j];
            }
        }
        for (j = 0; j < n_controls; j++) {
            r->U->e[n_nodes - 1][j] = r->Udata[n_nodes - 2]->e[rk_stages - 1][j];
        }
    }
}

int _remesh(RKOCP r) {
    r->only_compute_lte = YES;
    _RKOCPfhg(r->x0, NULL, NULL);
    r->max_lte = VectorMax(r->lte);
    return NO;
}

// IMPORTANT: We need rk_c[rk_stages-1] = 1.0
void _getRKcoefficients0(RKOCP r) {
    // set rk_a, rk_b, rk_c, rk_e, rk_stages
    // Bogacki-Shampine 4 stage 3(2)
    // FSAL 4-th stage added for error estimation
    // (2rd-oder error estimator)
    int i;
    r->rkMethodName = "Bogacki-Shampine-3(2)";
    r->rk_stages = 4;
    r->rk_order = 3;
    r->rk_a = MatrixNew(r->rk_stages, r->rk_stages);
    r->rk_b = VectorNew(r->rk_stages);
    r->rk_c = VectorNew(r->rk_stages);
    r->rk_e = VectorNew(r->rk_stages);
    double **rk_a, *rk_b, *rk_c, *rk_e;
    rk_a = r->rk_a->e;
    rk_b = r->rk_b->e;
    rk_c = r->rk_c->e;
    rk_e = r->rk_e->e;

    rk_a[1][ 0] = 1.0 / 2.0;

    rk_a[2][ 1] = 3.0 / 4.0;

    rk_a[3][ 0] = 2.0 / 9.0;
    rk_a[3][ 1] = 1.0 / 3.0;
    rk_a[3][ 2] = 4.0 / 9.0;

    rk_c[0] = 0.0;
    rk_c[1] = 0.5;
    rk_c[2] = 0.75;
    rk_c[3] = 1.0;

    rk_b[0] = 2.0 / 9.0;
    rk_b[1] = 1.0 / 3.0;
    rk_b[2] = 4.0 / 9.0;

    rk_e[0] = 7.0 / 24.0;
    rk_e[1] = 1.0 / 4.0;
    rk_e[2] = 1.0 / 3.0;
    rk_e[3] = 1.0 / 8.0;
    for (i = 0; i < r->rk_stages; i++) {
        rk_e[i] -= rk_b[i];
    }
}

void _getRKcoefficients1(RKOCP r) {
    // set rk_a, rk_b, rk_c, rk_e, r->rk_stages
    // Classical 4-th order Runge-Kutta method
    // FSAL 5-th stage added for error estimation
    // (3rd-oder error estimator)

    double **rk_a, *rk_b, *rk_c, *rk_e;
    r->rkMethodName = "Classical-4-th-order-Runge-Kutta-4(3)";
    r->rk_stages = 5;
    r->rk_order = 4;
    r->rk_a = MatrixNew(r->rk_stages, r->rk_stages);
    r->rk_b = VectorNew(r->rk_stages);
    r->rk_c = VectorNew(r->rk_stages);
    r->rk_e = VectorNew(r->rk_stages);
    rk_a = r->rk_a->e;
    rk_b = r->rk_b->e;
    rk_c = r->rk_c->e;
    rk_e = r->rk_e->e;

    rk_a[1][ 0] = 0.5;

    rk_a[2][ 1] = 0.5;

    rk_a[3][ 2] = 1.0;

    rk_a[4][ 0] = 1.0 / 6.0;
    rk_a[4][ 1] = 1.0 / 3.0;
    rk_a[4][ 2] = 1.0 / 3.0;
    rk_a[4][ 3] = 1.0 / 6.0;

    rk_c[0] = 0.0;
    rk_c[1] = 0.5;
    rk_c[2] = 0.5;
    rk_c[3] = 1.0;
    rk_c[4] = 1.0;

    rk_b[0] = 1.0 / 6.0;
    rk_b[1] = 1.0 / 3.0;
    rk_b[2] = 1.0 / 3.0;
    rk_b[3] = 1.0 / 6.0;

    rk_e[3] = 1.0 / 6.0;
    rk_e[4] = -1.0 / 6.0;
}

void _getRKcoefficients2(RKOCP r) {
    // set rk_a, rk_b, rk_c, rk_e, r->rk_stages
    // Dormand-Prince 5(4)

    double **rk_a, *rk_b, *rk_c, *rk_e;
    r->rkMethodName = "Dormand-Prince-5(4)";
    r->rk_stages = 7;
    r->rk_order = 5;
    r->rk_a = MatrixNew(r->rk_stages, r->rk_stages);
    r->rk_b = VectorNew(r->rk_stages);
    r->rk_c = VectorNew(r->rk_stages);
    r->rk_e = VectorNew(r->rk_stages);
    rk_a = r->rk_a->e;
    rk_b = r->rk_b->e;
    rk_c = r->rk_c->e;
    rk_e = r->rk_e->e;

    rk_c[1] = 0.2;
    rk_c[2] = 0.3;
    rk_c[3] = 0.8;
    rk_c[4] = 8.0 / 9.0;
    rk_c[5] = 1.0;
    rk_c[6] = 1.0;

    rk_a[1][ 0] = 0.2;

    rk_a[2][ 0] = 3.0 / 40.0;
    rk_a[2][ 1] = 9.0 / 40.0;
    rk_a[3][ 0] = 44.0 / 45.0;
    rk_a[3][ 1] = -56.0 / 15.0;
    rk_a[3][ 2] = 32.0 / 9.0;
    rk_a[4][ 0] = 19372.0 / 6561.0;
    rk_a[4][ 1] = -25360.0 / 2187.0;
    rk_a[4][ 2] = 64448.0 / 6561.0;
    rk_a[4][ 3] = -212.0 / 729.0;
    rk_a[5][ 0] = 9017.0 / 3168.0;
    rk_a[5][ 1] = -355.0 / 33.0;
    rk_a[5][ 2] = 46732.0 / 5247.0;
    rk_a[5][ 3] = 49.0 / 176.0;
    rk_a[5][ 4] = -5103.0 / 18656.0;
    rk_a[6][ 0] = 35.0 / 384.0;
    rk_a[6][ 2] = 500.0 / 1113.0;
    rk_a[6][ 3] = 125.0 / 192.0;
    rk_a[6][ 4] = -2187.0 / 6784.0;
    rk_a[6][ 5] = 11.0 / 84.0;

    rk_b[0] = 35.0 / 384.0;
    rk_b[2] = 500.0 / 1113.0;
    rk_b[3] = 125.0 / 192.0;
    rk_b[4] = -2187.0 / 6784.0;
    rk_b[5] = 11.0 / 84.0;

    rk_e[0] = 71.0 / 57600.0;
    rk_e[2] = -71.0 / 16695.0;
    rk_e[3] = 71.0 / 1920.0;
    rk_e[4] = -17253.0 / 339200.0;
    rk_e[5] = 22.0 / 525.0;
    rk_e[6] = -1.0 / 40.0;
}

void _getRKcoefficients3(RKOCP r) {
    // set rk_a, rk_b, rk_c, rk_e, r->rk_stages
    // Verner robust 6(5)
    // http://www.math.sfu.ca/~jverner/
    int i;
    double **rk_a, *rk_b, *rk_c, *rk_e;
    r->rkMethodName = "Verner-robust-6(5)";
    r->rk_stages = 9;
    r->rk_order = 6;
    r->rk_a = MatrixNew(r->rk_stages, r->rk_stages);
    r->rk_b = VectorNew(r->rk_stages);
    r->rk_c = VectorNew(r->rk_stages);
    r->rk_e = VectorNew(r->rk_stages);
    rk_a = r->rk_a->e;
    rk_b = r->rk_b->e;
    rk_c = r->rk_c->e;
    rk_e = r->rk_e->e;

    rk_c[0] = 0;
    rk_c[1] = 9.0 / 50.0;
    rk_c[2] = 1.0 / 6.0;
    rk_c[3] = 1.0 / 4.0;
    rk_c[4] = 53.0 / 100.0;
    rk_c[5] = 3.0 / 5.0;
    rk_c[6] = 4.0 / 5.0;
    rk_c[7] = 1.0;
    rk_c[8] = 1.0;

    rk_a[1][0] = 9.0 / 50.0;

    rk_a[2][0] = 29.0 / 324.0;
    rk_a[2][1] = 25.0 / 324.0;

    rk_a[3][0] = 1.0 / 16.0;
    rk_a[3][1] = 0.0;
    rk_a[3][2] = 3.0 / 16.0;

    rk_a[4][0] = 79129.0 / 250000.0;
    rk_a[4][1] = 0.0;
    rk_a[4][2] = -261237.0 / 250000.0;
    rk_a[4][3] = 19663.0 / 15625.0;

    rk_a[5][0] = 1336883.0 / 4909125.0;
    rk_a[5][1] = 0.0;
    rk_a[5][2] = -25476.0 / 30875.0;
    rk_a[5][3] = 194159.0 / 185250.0;
    rk_a[5][4] = 8225.0 / 78546.0;

    rk_a[6][0] = -2459386.0 / 14727375.0;
    rk_a[6][1] = 0.0;
    rk_a[6][2] = 19504.0 / 30875.0;
    rk_a[6][3] = 2377474.0 / 13615875.0;
    rk_a[6][4] = -6157250.0 / 5773131.0;
    rk_a[6][5] = 902.0 / 735.0;

    rk_a[7][0] = 2699.0 / 7410.0;
    rk_a[7][1] = 0.0;
    rk_a[7][2] = -252.0 / 1235.0;
    rk_a[7][3] = -1393253.0 / 3993990.0;
    rk_a[7][4] = 236875.0 / 72618.0;
    rk_a[7][5] = -135.0 / 49.0;
    rk_a[7][6] = 15.0 / 22.0;

    rk_a[8][0] = 11.0 / 144.0;
    rk_a[8][1] = 0.0;
    rk_a[8][2] = 0.0;
    rk_a[8][3] = 256.0 / 693.0;
    rk_a[8][4] = 0.0;
    rk_a[8][5] = 125.0 / 504.0;
    rk_a[8][6] = 125.0 / 528.0;
    rk_a[8][7] = 5.0 / 72.0;

    rk_b[0] = 11.0 / 144.0;
    rk_b[1] = 0.0;
    rk_b[2] = 0.0;
    rk_b[3] = 256.0 / 693.0;
    rk_b[4] = 0.0;
    rk_b[5] = 125.0 / 504.0;
    rk_b[6] = 125.0 / 528.0;
    rk_b[7] = 5.0 / 72.0;
    rk_b[8] = 0.0;

    rk_e[0] = 28.0 / 477.0;
    rk_e[1] = 0.0;
    rk_e[2] = 0.0;
    rk_e[3] = 212.0 / 441.0;
    rk_e[4] = -312500.0 / 366177.0;
    rk_e[5] = 2125.0 / 1764.0;
    rk_e[6] = 0.0;
    rk_e[7] = -2105.0 / 35532.0;
    rk_e[8] = 2995.0 / 17766.0;
    for (i = 0; i < r->rk_stages; i++) {

        rk_e[i] -= rk_b[i];
    }
}

void _getRKcoefficients4(RKOCP r) {
    // set rk_a, rk_b, rk_c, rk_e, r->rk_stages
    // Verner 'robust' 7(6)
    // http://www.math.sfu.ca/~jverner/

    double **rk_a, *rk_b, *rk_c, *rk_e;
    r->rkMethodName = "Verner-robust-7(6)";
    r->rk_stages = 10;
    r->rk_order = 7;
    r->rk_a = MatrixNew(r->rk_stages, r->rk_stages);
    r->rk_b = VectorNew(r->rk_stages);
    r->rk_c = VectorNew(r->rk_stages);
    r->rk_e = VectorNew(r->rk_stages);
    rk_a = r->rk_a->e;
    rk_b = r->rk_b->e;
    rk_c = r->rk_c->e;
    rk_e = r->rk_e->e;

    rk_c[1] = 5.000000000000e-03;
    rk_c[2] = 1.088888888889e-01;
    rk_c[3] = 1.633333333333e-01;
    rk_c[4] = 4.550000000000e-01;
    rk_c[5] = 6.059617471463e-01;
    rk_c[6] = 8.350000000000e-01;
    rk_c[7] = 9.150000000000e-01;
    rk_c[8] = 1.000000000000e+00;
    rk_c[9] = 1.000000000000e+00;

    rk_b[0] = 4.742583783371e-02;
    rk_b[3] = 2.562236165937e-01;
    rk_b[4] = 2.695137683307e-01;
    rk_b[5] = 1.268662240909e-01;
    rk_b[6] = 2.488722594206e-01;
    rk_b[7] = 3.074483740820e-03;
    rk_b[8] = 4.802380998950e-02;
    rk_e[0] = -5.940986559288e-05;
    rk_e[3] = 2.294907067993e-04;
    rk_e[4] = -1.071012479935e-03;
    rk_e[5] = 1.810037246668e-03;
    rk_e[6] = -3.172427816838e-03;
    rk_e[7] = 3.074483740820e-03;
    rk_e[8] = 4.802380998950e-02;
    rk_e[9] = -4.883497152142e-02;
    rk_a[1][0] = 5.000000000000e-03;
    rk_a[2][0] = -1.076790123457e+00;
    rk_a[2][1] = 1.185679012346e+00;
    rk_a[3][0] = 4.083333333333e-02;
    rk_a[3][2] = 1.225000000000e-01;
    rk_a[4][0] = 6.360714285714e-01;
    rk_a[4][2] = -2.444464285714e+00;
    rk_a[4][3] = 2.263392857143e+00;
    rk_a[5][0] = -2.535121107935e+00;
    rk_a[5][2] = 1.029937465445e+01;
    rk_a[5][3] = -7.951303288599e+00;
    rk_a[5][4] = 7.930114892310e-01;
    rk_a[6][0] = 1.001876581252e+00;
    rk_a[6][2] = -4.166571282442e+00;
    rk_a[6][3] = 3.834343292913e+00;
    rk_a[6][4] = -5.023333356071e-01;
    rk_a[6][5] = 6.676847438842e-01;
    rk_a[7][0] = 2.725501835463e+01;
    rk_a[7][2] = -4.200461727841e+01;
    rk_a[7][3] = -1.053571312662e+01;
    rk_a[7][4] = 8.049553671141e+01;
    rk_a[7][5] = -6.734388227179e+01;
    rk_a[7][6] = 1.304865761078e+01;
    rk_a[8][0] = -3.039737805711e+00;
    rk_a[8][2] = 1.013816141033e+01;
    rk_a[8][3] = -6.429305674865e+00;
    rk_a[8][4] = -1.586437148341e+00;
    rk_a[8][5] = 1.892178184197e+00;
    rk_a[8][6] = 1.969933540761e-02;
    rk_a[8][7] = 5.441698982793e-03;
    rk_a[9][0] = -1.444951891678e+00;
    rk_a[9][2] = 8.031891385996e+00;
    rk_a[9][3] = -7.583174166340e+00;
    rk_a[9][4] = 3.581616935319e+00;
    rk_a[9][5] = -2.436972263220e+00;
    rk_a[9][6] = 8.515899999233e-01;
}

void _getRKcoefficients5(RKOCP r) {
    // set rk_a, rk_b, rk_c, rk_e, r->rk_stages
    // Verner 'robust' 13 stage 8(7)
    // http://www.math.sfu.ca/~jverner/

    double **rk_a, *rk_b, *rk_c, *rk_e;
    r->rkMethodName = "Verner-robust-8(7)";
    r->rk_stages = 13;
    r->rk_order = 8;
    r->rk_a = MatrixNew(r->rk_stages, r->rk_stages);
    r->rk_b = VectorNew(r->rk_stages);
    r->rk_c = VectorNew(r->rk_stages);
    r->rk_e = VectorNew(r->rk_stages);
    rk_a = r->rk_a->e;
    rk_b = r->rk_b->e;
    rk_c = r->rk_c->e;
    rk_e = r->rk_e->e;

    rk_c[1] = 2.500000000000e-01;
    rk_c[2] = 1.128884514436e-01;
    rk_c[3] = 1.693326771654e-01;
    rk_c[4] = 4.240000000000e-01;
    rk_c[5] = 5.090000000000e-01;
    rk_c[6] = 8.670000000000e-01;
    rk_c[7] = 1.500000000000e-01;
    rk_c[8] = 7.090680365139e-01;
    rk_c[9] = 3.200000000000e-01;
    rk_c[10] = 4.500000000000e-01;
    rk_c[11] = 1.000000000000e+00;
    rk_c[12] = 1.000000000000e+00;

    rk_b[0] = 4.472956466670e-02;
    rk_b[5] = 1.569103352771e-01;
    rk_b[6] = 1.846097340815e-01;
    rk_b[7] = 2.251638060209e-01;
    rk_b[8] = 1.479461565197e-01;
    rk_b[9] = 7.605554244496e-02;
    rk_b[10] = 1.227729023502e-01;
    rk_b[11] = 4.181195863899e-02;
    rk_e[0] = -1.117546733800e-03;
    rk_e[5] = -1.054085787644e-01;
    rk_e[6] = -7.083989297010e-03;
    rk_e[7] = 8.072082741844e-03;
    rk_e[8] = 2.056426027137e-02;
    rk_e[9] = -3.904976140870e-02;
    rk_e[10] = 1.227729023502e-01;
    rk_e[11] = 4.181195863899e-02;
    rk_e[12] = -4.056132779844e-02;
    rk_a[1][0] = 2.500000000000e-01;
    rk_a[2][0] = 8.740084650492e-02;
    rk_a[2][1] = 2.548760493865e-02;
    rk_a[3][0] = 4.233316929134e-02;
    rk_a[3][2] = 1.269995078740e-01;
    rk_a[4][0] = 4.260950588874e-01;
    rk_a[4][2] = -1.598795284659e+00;
    rk_a[4][3] = 1.596700225772e+00;
    rk_a[5][0] = 5.071933729671e-02;
    rk_a[5][3] = 2.543337726460e-01;
    rk_a[5][4] = 2.039468900573e-01;
    rk_a[6][0] = -2.900037471752e-01;
    rk_a[6][3] = 1.344187391026e+00;
    rk_a[6][4] = -2.864777943361e+00;
    rk_a[6][5] = 2.677594299511e+00;
    rk_a[7][0] = 9.853501133799e-02;
    rk_a[7][4] = 2.219268063075e-01;
    rk_a[7][5] = -1.814062291181e-01;
    rk_a[7][6] = 1.094441147256e-02;
    rk_a[8][0] = 3.871105254573e-01;
    rk_a[8][3] = -1.442445497486e+00;
    rk_a[8][4] = 2.905398189070e+00;
    rk_a[8][5] = -1.853771069630e+00;
    rk_a[8][6] = 1.400364809873e-01;
    rk_a[8][7] = 5.727394081150e-01;
    rk_a[9][0] = -1.612440344444e-01;
    rk_a[9][3] = -1.733960295736e-01;
    rk_a[9][4] = -1.301289281407e+00;
    rk_a[9][5] = 1.137950375174e+00;
    rk_a[9][6] = -3.174764966397e-02;
    rk_a[9][7] = 9.335129382493e-01;
    rk_a[9][8] = -8.378631833473e-02;
    rk_a[10][0] = -1.919944488159e-02;
    rk_a[10][3] = 2.733085726526e-01;
    rk_a[10][4] = -6.753497320694e-01;
    rk_a[10][5] = 3.415184981385e-01;
    rk_a[10][6] = -6.795006480338e-02;
    rk_a[10][7] = 9.659175224762e-02;
    rk_a[10][8] = 1.325308251118e-01;
    rk_a[10][9] = 3.685495936039e-01;
    rk_a[11][0] = 6.091877403645e-01;
    rk_a[11][3] = -2.272569085898e+00;
    rk_a[11][4] = 4.757898342694e+00;
    rk_a[11][5] = -5.516106706693e+00;
    rk_a[11][6] = 2.900596369680e-01;
    rk_a[11][7] = 5.691423963359e-01;
    rk_a[11][8] = 7.926795760332e-01;
    rk_a[11][9] = 1.547372045329e-01;
    rk_a[11][10] = 1.614970895662e+00;
    rk_a[12][0] = 8.873576220853e-01;
    rk_a[12][3] = -2.975459782109e+00;
    rk_a[12][4] = 5.600717009488e+00;
    rk_a[12][5] = -5.915607450537e+00;
    rk_a[12][6] = 2.202968915613e-01;
    rk_a[12][7] = 1.015509782446e-01;
    rk_a[12][8] = 1.151434564739e+00;
    rk_a[12][9] = 1.929710166527e+00;
}

void _getRKcoefficients6(RKOCP r) {
    // set rk_a, rk_b, rk_c, rk_e, rk_stages
    // Verner 'robust' 16 stage 9(8)
    // http://www.math.sfu.ca/~jverner/
    double **rk_a, *rk_b, *rk_c, *rk_e;
    r->rkMethodName = "Verner-robust-9(8)";
    r->rk_stages = 16;
    r->rk_order = 9;
    r->rk_a = MatrixNew(r->rk_stages, r->rk_stages);
    r->rk_b = VectorNew(r->rk_stages);
    r->rk_c = VectorNew(r->rk_stages);
    r->rk_e = VectorNew(r->rk_stages);
    rk_a = r->rk_a->e;
    rk_b = r->rk_b->e;
    rk_c = r->rk_c->e;
    rk_e = r->rk_e->e;

    rk_c[1] = 4.000000000000e-02;
    rk_c[2] = 9.648736013787e-02;
    rk_c[3] = 1.447310402068e-01;
    rk_c[4] = 5.760000000000e-01;
    rk_c[5] = 2.272326564619e-01;
    rk_c[6] = 5.407673435381e-01;
    rk_c[7] = 6.400000000000e-01;
    rk_c[8] = 4.800000000000e-01;
    rk_c[9] = 6.754000000000e-02;
    rk_c[10] = 2.500000000000e-01;
    rk_c[11] = 6.770920153543e-01;
    rk_c[12] = 8.115000000000e-01;
    rk_c[13] = 9.060000000000e-01;
    rk_c[14] = 1.000000000000e+00;
    rk_c[15] = 1.000000000000e+00;

    rk_b[0] = 1.458885278406e-02;
    rk_b[7] = 2.024197887889e-03;
    rk_b[8] = 2.178047084570e-01;
    rk_b[9] = 1.274895340854e-01;
    rk_b[10] = 2.244617745463e-01;
    rk_b[11] = 1.787254491260e-01;
    rk_b[12] = 7.594344758097e-02;
    rk_b[13] = 1.294845879198e-01;
    rk_b[14] = 2.947744761262e-02;
    rk_b[15] = 0.000000000000e+00;
    rk_e[0] = -5.757813768189e-03;
    rk_e[7] = -1.067593453095e+00;
    rk_e[8] = 1.409963613439e-01;
    rk_e[9] = 1.441171539691e-02;
    rk_e[10] = -3.079696125188e-02;
    rk_e[11] = 1.161315257818e+00;
    rk_e[12] = -3.222111348612e-01;
    rk_e[13] = 1.294845879198e-01;
    rk_e[14] = 2.947744761262e-02;
    rk_e[15] = -4.932600711507e-02;
    rk_a[1][0] = 4.000000000000e-02;
    rk_a[2][0] = -1.988527319182e-02;
    rk_a[2][1] = 1.163726333297e-01;
    rk_a[3][0] = 3.618276005170e-02;
    rk_a[3][2] = 1.085482801551e-01;
    rk_a[4][0] = 2.272114264290e+00;
    rk_a[4][2] = -8.526886447976e+00;
    rk_a[4][3] = 6.830772183686e+00;
    rk_a[5][0] = 5.094385535389e-02;
    rk_a[5][3] = 1.755865049809e-01;
    rk_a[5][4] = 7.022961270757e-04;
    rk_a[6][0] = 1.424783668683e-01;
    rk_a[6][3] = -3.541799434669e-01;
    rk_a[6][4] = 7.595315450295e-02;
    rk_a[6][5] = 6.765157656337e-01;
    rk_a[7][0] = 7.111111111111e-02;
    rk_a[7][5] = 3.279909287606e-01;
    rk_a[7][6] = 2.408979601283e-01;
    rk_a[8][0] = 7.125000000000e-02;
    rk_a[8][5] = 3.268842451575e-01;
    rk_a[8][6] = 1.156157548425e-01;
    rk_a[8][7] = -3.375000000000e-02;
    rk_a[9][0] = 4.822677322466e-02;
    rk_a[9][5] = 3.948559980495e-02;
    rk_a[9][6] = 1.058851161935e-01;
    rk_a[9][7] = -2.152006320474e-02;
    rk_a[9][8] = -1.045374260183e-01;
    rk_a[10][0] = -2.609113435755e-02;
    rk_a[10][5] = 3.333333333333e-02;
    rk_a[10][6] = -1.652504006638e-01;
    rk_a[10][7] = 3.434664118369e-02;
    rk_a[10][8] = 1.595758283215e-01;
    rk_a[10][9] = 2.140857321828e-01;
    rk_a[11][0] = -3.628423396256e-02;
    rk_a[11][5] = -1.096167597427e+00;
    rk_a[11][6] = 1.826035504321e-01;
    rk_a[11][7] = 7.082254444171e-02;
    rk_a[11][8] = -2.313647018482e-02;
    rk_a[11][9] = 2.711204726321e-01;
    rk_a[11][10] = 1.308133749423e+00;
    rk_a[12][0] = -5.074635056417e-01;
    rk_a[12][5] = -6.631342198657e+00;
    rk_a[12][6] = -2.527480100909e-01;
    rk_a[12][7] = -4.952612380036e-01;
    rk_a[12][8] = 2.932525545254e-01;
    rk_a[12][9] = 1.440108693768e+00;
    rk_a[12][10] = 6.237934498647e+00;
    rk_a[12][11] = 7.270192054527e-01;
    rk_a[13][0] = 6.130118256956e-01;
    rk_a[13][5] = 9.088803891640e+00;
    rk_a[13][6] = -4.073788156293e-01;
    rk_a[13][7] = 1.790733389490e+00;
    rk_a[13][8] = 7.149271667618e-01;
    rk_a[13][9] = -1.438580857842e+00;
    rk_a[13][10] = -8.263329312065e+00;
    rk_a[13][11] = -1.537570570809e+00;
    rk_a[13][12] = 3.453832827565e-01;
    rk_a[14][0] = -1.211697910344e+00;
    rk_a[14][5] = -1.905581871560e+01;
    rk_a[14][6] = 1.263060675390e+00;
    rk_a[14][7] = -6.913916969178e+00;
    rk_a[14][8] = -6.764622665095e-01;
    rk_a[14][9] = 3.367860445027e+00;
    rk_a[14][10] = 1.800675164313e+01;
    rk_a[14][11] = 6.838828926794e+00;
    rk_a[14][12] = -1.031516451922e+00;
    rk_a[14][13] = 4.129106232131e-01;
    rk_a[15][0] = 2.157389007494e+00;
    rk_a[15][5] = 2.380712219810e+01;
    rk_a[15][6] = 8.862779249217e-01;
    rk_a[15][7] = 1.313913039760e+01;
    rk_a[15][8] = -2.604415709288e+00;
    rk_a[15][9] = -5.193859949784e+00;
    rk_a[15][10] = -2.041234071154e+01;
    rk_a[15][11] = -1.230085625251e+01;
    rk_a[15][12] = 1.521553095009e+00;
}

void _getRKcoefficients7(RKOCP r) {
    // set rk_a, rk_b, rk_c, rk_e, r->rk_stages
    // Dormand-Prince 8(7)
    int i;
    double **rk_a, *rk_b, *rk_c, *rk_e;
    r->rkMethodName = "Dormand-Prince-8(7)";
    r->rk_stages = 13;
    r->rk_order = 8;
    r->rk_a = MatrixNew(r->rk_stages, r->rk_stages);
    r->rk_b = VectorNew(r->rk_stages);
    r->rk_c = VectorNew(r->rk_stages);
    r->rk_e = VectorNew(r->rk_stages);
    rk_a = r->rk_a->e;
    rk_b = r->rk_b->e;
    rk_c = r->rk_c->e;
    rk_e = r->rk_e->e;

    rk_a[1][0] = 5.55555555555555555555555555556e-2;

    rk_a[2][0] = 2.08333333333333333333333333333e-2;
    rk_a[2][1] = 6.25e-2;

    rk_a[3][0] = 3.125e-2;
    rk_a[3][2] = 9.375e-2;

    rk_a[4][0] = 3.125e-1;
    rk_a[4][2] = -1.171875e0;
    rk_a[4][3] = 1.171875e0;

    rk_a[5][0] = 3.75e-2;
    rk_a[5][3] = 1.875e-1;
    rk_a[5][4] = 1.5e-1;

    rk_a[6][0] = 4.79101371111111111111111111111e-2;
    rk_a[6][3] = 1.12248712777777777777777777778e-1;
    rk_a[6][4] = -2.55056737777777777777777777778e-2;
    rk_a[6][5] = 1.28468238888888888888888888889e-2;

    rk_a[7][0] = 1.6917989787292281181431107136e-2;
    rk_a[7][3] = 3.87848278486043169526545744159e-1;
    rk_a[7][4] = 3.59773698515003278967008896348e-2;
    rk_a[7][5] = 1.96970214215666060156715256072e-1;
    rk_a[7][6] = -1.72713852340501838761392997002e-1;

    rk_a[8][0] = 6.90957533591923006485645489846e-2;
    rk_a[8][3] = -6.34247976728854151882807874972e-1;
    rk_a[8][4] = -1.61197575224604080366876923982e-1;
    rk_a[8][5] = 1.38650309458825255419866950133e-1;
    rk_a[8][6] = 9.4092861403575626972423968413e-1;
    rk_a[8][7] = 2.11636326481943981855372117132e-1;

    rk_a[9][0] = 1.83556996839045385489806023537e-1;
    rk_a[9][3] = -2.46876808431559245274431575997e0;
    rk_a[9][4] = -2.91286887816300456388002572804e-1;
    rk_a[9][5] = -2.6473020233117375688439799466e-2;
    rk_a[9][6] = 2.84783876419280044916451825422e0;
    rk_a[9][7] = 2.81387331469849792539403641827e-1;
    rk_a[9][8] = 1.23744899863314657627030212664e-1;

    rk_a[10][0] = -1.21542481739588805916051052503e0;
    rk_a[10][3] = 1.66726086659457724322804132886e1;
    rk_a[10][4] = 9.15741828416817960595718650451e-1;
    rk_a[10][5] = -6.05660580435747094755450554309e0;
    rk_a[10][6] = -1.60035735941561781118417064101e1;
    rk_a[10][7] = 1.4849303086297662557545391898e1;
    rk_a[10][8] = -1.33715757352898493182930413962e1;
    rk_a[10][9] = 5.13418264817963793317325361166e0;

    rk_a[11][0] = 2.58860916438264283815730932232e-1;
    rk_a[11][3] = -4.77448578548920511231011750971e0;
    rk_a[11][4] = -4.3509301377703250944070041181e-1;
    rk_a[11][5] = -3.04948333207224150956051286631e0;
    rk_a[11][6] = 5.57792003993609911742367663447e0;
    rk_a[11][7] = 6.15583158986104009733868912669e0;
    rk_a[11][8] = -5.06210458673693837007740643391e0;
    rk_a[11][9] = 2.19392617318067906127491429047e0;
    rk_a[11][10] = 1.34627998659334941535726237887e-1;

    rk_a[12][0] = 8.22427599626507477963168204773e-1;
    rk_a[12][3] = -1.16586732572776642839765530355e1;
    rk_a[12][4] = -7.57622116690936195881116154088e-1;
    rk_a[12][5] = 7.13973588159581527978269282765e-1;
    rk_a[12][6] = 1.20757749868900567395661704486e1;
    rk_a[12][7] = -2.12765911392040265639082085897e0;
    rk_a[12][8] = 1.99016620704895541832807169835e0;
    rk_a[12][9] = -2.34286471544040292660294691857e-1;
    rk_a[12][10] = 1.7589857770794226507310510589e-1;

    rk_b[0] = 4.17474911415302462220859284685e-2;
    rk_b[5] = -5.54523286112393089615218946547e-2;
    rk_b[6] = 2.39312807201180097046747354249e-1;
    rk_b[7] = 7.0351066940344302305804641089e-1;
    rk_b[8] = -7.59759613814460929884487677085e-1;
    rk_b[9] = 6.60563030922286341461378594838e-1;
    rk_b[10] = 1.58187482510123335529614838601e-1;
    rk_b[11] = -2.38109538752862804471863555306e-1;
    rk_b[12] = 2.5e-1;

    rk_e[0] = 2.9553213676353496981964883112e-2;
    rk_e[5] = -8.28606276487797039766805612689e-1;
    rk_e[6] = 3.11240900051118327929913751627e-1;
    rk_e[7] = 2.46734519059988698196468570407e0;
    rk_e[8] = -2.54694165184190873912738007542e0;
    rk_e[9] = 1.44354858367677524030187495069e0;
    rk_e[10] = 7.94155958811272872713019541622e-2;
    rk_e[11] = 4.44444444444444444444444444445e-2;
    for (i = 0; i < r->rk_stages; i++) {
        rk_e[i] -= rk_b[i];
    }

    rk_c[1] = 5.55555555555555555555555555556e-2;
    rk_c[2] = 8.33333333333333333333333333334e-2;
    rk_c[3] = 1.25e-1;
    rk_c[4] = 3.125e-1;
    rk_c[5] = 3.75e-1;
    rk_c[6] = 1.475e-1;
    rk_c[7] = 4.65e-1;
    rk_c[8] = 5.64865451382259575398358501426e-1;
    rk_c[9] = 6.5e-1;
    rk_c[10] = 9.24656277640504446745013574318e-1;
    rk_c[11] = 1.0e0;
    rk_c[12] = 1.0e0;
}

void _getRKcoefficients8(RKOCP r) {
    // set rk_a, rk_b, rk_c, rk_e, r->rk_stages
    // Runge-Kutta-Verner 6(5) DVERK
    int i;
    double **rk_a, *rk_b, *rk_c, *rk_e;
    r->rkMethodName = "Runge-Kutta-Verner-6(5)-DVERK";
    r->rk_stages = 8;
    r->rk_order = 6;
    r->rk_a = MatrixNew(r->rk_stages, r->rk_stages);
    r->rk_b = VectorNew(r->rk_stages);
    r->rk_c = VectorNew(r->rk_stages);
    r->rk_e = VectorNew(r->rk_stages);
    rk_a = r->rk_a->e;
    rk_b = r->rk_b->e;
    rk_c = r->rk_c->e;
    rk_e = r->rk_e->e;

    rk_a[1][0] = 1.0 / 6.0;

    rk_a[2][0] = 4.0 / 75.0;
    rk_a[2][1] = 16.0 / 75.0;

    rk_a[3][0] = 5.0 / 6.0;
    rk_a[3][1] = -8.0 / 3.0;
    rk_a[3][2] = 5.0 / 2.0;

    rk_a[4][0] = -165.0 / 64.0;
    rk_a[4][1] = 55.0 / 6.0;
    rk_a[4][2] = -425.0 / 64.0;
    rk_a[4][3] = 85.0 / 96.0;

    rk_a[5][0] = 12.0 / 5.0;
    rk_a[5][1] = -8.0;
    rk_a[5][2] = 4015.0 / 612.0;
    rk_a[5][3] = -11.0 / 36.0;
    rk_a[5][4] = 88.0 / 255.0;

    rk_a[6][0] = -8263.0 / 15000.0;
    rk_a[6][1] = 124.0 / 75.0;
    rk_a[6][2] = -643.0 / 680.0;
    rk_a[6][3] = -81.0 / 250.0;
    rk_a[6][4] = 2484.0 / 10625.0;

    rk_a[7][0] = 3501.0 / 1720.0;
    rk_a[7][1] = -300.0 / 43.0;
    rk_a[7][2] = 297275.0 / 52632.0;
    rk_a[7][3] = -319.0 / 2322.0;
    rk_a[7][4] = 24068.0 / 84065.0;
    rk_a[7][5] = 0.0;
    rk_a[7][6] = 3850.0 / 26703.0;

    rk_c[0] = 0.0;
    rk_c[1] = 1.0 / 6.0;
    rk_c[2] = 4.0 / 15.0;
    rk_c[3] = 2.0 / 3.0;
    rk_c[4] = 5.0 / 6.0;
    rk_c[5] = 1.0;
    rk_c[6] = 1.0 / 15.0;
    rk_c[7] = 1.0;

    rk_b[0] = 3.0 / 40.0;
    rk_b[1] = 0.0;
    rk_b[2] = 875.0 / 2244.0;
    rk_b[3] = 23.0 / 72.0;
    rk_b[4] = 264.0 / 1955.0;
    rk_b[5] = 0.0;
    rk_b[6] = 125.0 / 11592.0;
    rk_b[7] = 43.0 / 616.0;

    rk_e[0] = 13.0 / 160.0;
    rk_e[1] = 0.0;
    rk_e[2] = 2375.0 / 5984.0;
    rk_e[3] = 5.0 / 16.0;
    rk_e[4] = 12.0 / 85.0;
    rk_e[5] = 3.0 / 44.0;
    rk_e[6] = 0.0;
    rk_e[7] = 0.0;
    for (i = 0; i < r->rk_stages; i++) {
        rk_e[i] -= rk_b[i];
    }

}

void _getRKcoefficients9(RKOCP r) {
    // set rk_a, rk_b, rk_c, rk_e, r->rk_stages
    // Runge-Kutta-Fehlberg 7(8)
    int i;
    double **rk_a, *rk_b, *rk_c, *rk_e;
    r->rkMethodName = "Runge-Kutta-Fehlberg-7(8)";
    r->rk_stages = 13;
    r->rk_order = 8;
    r->rk_a = MatrixNew(r->rk_stages, r->rk_stages);
    r->rk_b = VectorNew(r->rk_stages);
    r->rk_c = VectorNew(r->rk_stages);
    r->rk_e = VectorNew(r->rk_stages);
    rk_a = r->rk_a->e;
    rk_b = r->rk_b->e;
    rk_c = r->rk_c->e;
    rk_e = r->rk_e->e;

    rk_a[1][0] = 2.0 / 27.0;

    rk_a[2][0] = 1.0 / 36.0;
    rk_a[2][1] = 1.0 / 12.0;

    rk_a[3][0] = 1.0 / 24.0;
    rk_a[3][2] = 1.0 / 8.0;

    rk_a[4][0] = 5.0 / 12.0;
    rk_a[4][2] = -25.0 / 16.0;
    rk_a[4][3] = 25.0 / 16.0;

    rk_a[5][0] = 1.0 / 20.0;
    rk_a[5][3] = 1.0 / 4.0;
    rk_a[5][4] = 1.0 / 5.0;

    rk_a[6][0] = -25.0 / 108.0;
    rk_a[6][3] = 125.0 / 108.0;
    rk_a[6][4] = -65.0 / 27.0;
    rk_a[6][5] = 125.0 / 54.0;


    rk_a[7][0] = 31.0 / 300.0;
    rk_a[7][4] = 61.0 / 225.0;
    rk_a[7][5] = -2.0 / 9.0;
    rk_a[7][6] = 13.0 / 900.0;

    rk_a[8][0] = 2.0;
    rk_a[8][3] = -53.0 / 6.0;
    rk_a[8][4] = 704.0 / 45.0;
    rk_a[8][5] = -107.0 / 9.0;
    rk_a[8][6] = 67.0 / 90.0;
    rk_a[8][7] = 3.0;

    rk_a[9][0] = -91.0 / 108.0;
    rk_a[9][3] = 23.0 / 108.0;
    rk_a[9][4] = -976.0 / 135.0;
    rk_a[9][5] = 311.0 / 54.0;
    rk_a[9][6] = -19.0 / 60.0;
    rk_a[9][7] = 17.0 / 6.0;
    rk_a[9][8] = -1.0 / 12.0;

    rk_a[10][0] = 2383.0 / 4100.0;
    rk_a[10][3] = -341.0 / 164.0;
    rk_a[10][4] = 4496.0 / 1025.0;
    rk_a[10][5] = -301.0 / 82.0;
    rk_a[10][6] = 2133.0 / 4100.0;
    rk_a[10][7] = 45.0 / 82.0;
    rk_a[10][8] = 45.0 / 164.0;
    rk_a[10][9] = 18.0 / 41.0;

    rk_a[11][0] = 3.0 / 205;
    rk_a[11][5] = -6.0 / 41.0;
    rk_a[11][6] = -3.0 / 205.0;
    rk_a[11][7] = -3.0 / 41.0;
    rk_a[11][8] = 3.0 / 41.0;
    rk_a[11][9] = 6.0 / 41.0;

    rk_a[12][0] = -1777.0 / 4100.0;
    rk_a[12][3] = -341.0 / 164.0;
    rk_a[12][4] = 4496.0 / 1025.0;
    rk_a[12][5] = -289.0 / 82.0;
    rk_a[12][6] = 2193.0 / 4100.0;
    rk_a[12][7] = 51.0 / 82.0;
    rk_a[12][8] = 33.0 / 164.0;
    rk_a[12][9] = 19.0 / 41.0;
    rk_a[12][11] = 1.0;

    rk_c[1] = 2.0 / 27.0;
    rk_c[2] = 1.0 / 9.0;
    rk_c[3] = 1.0 / 6.0;
    rk_c[4] = 5.0 / 12.0;
    rk_c[5] = 0.5;
    rk_c[6] = 5.0 / 6.0;
    rk_c[7] = 1.0 / 6.0;
    rk_c[8] = 2.0 / 3.0;
    rk_c[9] = 1.0 / 3.0;
    rk_c[10] = 1.0;
    rk_c[12] = 1.0;


    rk_b[0] = 41.0 / 840.0;
    rk_b[5] = 34.0 / 105.0;
    rk_b[6] = 9.0 / 35.0;
    rk_b[7] = 9.0 / 35.0;
    rk_b[8] = 9.0 / 280.0;
    rk_b[9] = 9.0 / 280.0;
    rk_b[10] = 41.0 / 840.0;

    rk_e[5] = 34.0 / 105.0;
    rk_e[6] = 9.0 / 35.0;
    rk_e[7] = 9.0 / 35.0;
    rk_e[8] = 9.0 / 280.0;
    rk_e[9] = 9.0 / 280.0;
    rk_e[11] = 41.0 / 840.0;
    rk_e[12] = 41.0 / 840.0;
    for (i = 0; i < r->rk_stages; i++) {
        rk_e[i] -= rk_b[i];
    }
}

void _getRKcoefficients10(RKOCP r) {
    // set rk_a, rk_b, rk_c, rk_e, rk_stages
    // Merson 4("5")
    // (3rd/5th-oder error estimator)
    int i;
    double **rk_a, *rk_b, *rk_c, *rk_e;
    r->rkMethodName = "Runge-Kutta-Merson-4(\"5\")";
    r->rk_stages = 5;
    r->rk_order = 4;
    r->rk_a = MatrixNew(r->rk_stages, r->rk_stages);
    r->rk_b = VectorNew(r->rk_stages);
    r->rk_c = VectorNew(r->rk_stages);
    r->rk_e = VectorNew(r->rk_stages);
    rk_a = r->rk_a->e;
    rk_b = r->rk_b->e;
    rk_c = r->rk_c->e;
    rk_e = r->rk_e->e;

    rk_a[1][0] = 1.0 / 3.0;

    rk_a[2][0] = 1.0 / 6.0;
    rk_a[2][1] = 1.0 / 6.0;

    rk_a[3][0] = 1.0 / 8.0;
    rk_a[3][2] = 3.0 / 8.0;

    rk_a[4][0] = 1.0 / 2.0;
    rk_a[4][1] = 0.0;
    rk_a[4][2] = -3.0 / 2.0;
    rk_a[4][3] = 2.0;

    rk_c[0] = 0.0;
    rk_c[1] = 1.0 / 3.0;
    rk_c[2] = 1.0 / 3.0;
    rk_c[3] = 1.0 / 2.0;
    rk_c[4] = 1.0;

    rk_b[0] = 1.0 / 6.0;
    rk_b[1] = 0.0;
    rk_b[2] = 0.0;
    rk_b[3] = 2.0 / 3.0;
    rk_b[4] = 1.0 / 6.0;

    rk_e[0] = 1.0 / 10.0;
    rk_e[1] = 0.0;
    rk_e[2] = 3.0 / 10.0;
    rk_e[3] = 2.0 / 5.0;
    rk_e[4] = 1.0 / 5.0;
    for (i = 0; i < r->rk_stages; i++) {
        rk_e[i] -= rk_b[i];
    }
}

double _compute_LTE(RKOCP r, double *eta, double *y0, double *y) {
    // Compute the local truncation error
    int i;
    double w, xe, lte_sum = 0.0;
    for (i = 0; i < r->n_states + 1; i++) {
        w = r->tolerance + r->tolerance
                * fmax(fabs(y0[i]), fabs(y[i]));
        xe = eta[i] / w;
        lte_sum += xe * xe;
    }
    return sqrt(lte_sum / ((double) r->n_states + 1));
}

void _setUdata(RKOCP r, Vector x) {
    switch (r->c_type) {
        case CONSTANT:
            _setUdata0(r, x);
            break;
        case LINEAR:
            _setUdata2(r, x);
            break;
        case CUBIC:
            _setUdata8(r, x);
            break;
    }
}

void _setUdata0(RKOCP r, Vector x) {
    int i, j, k;
    int k0 = r->n_parameters + r->n_states;
    for (i = 0; i < r->n_nodes - 1; i++) {
        for (k = 0; k < r->n_controls; k++) r->U->e[i][k] = x->e[k0 + k];
        for (j = 0; j < r->rk_stages; j++) {
            for (k = 0; k < r->n_controls; k++) r->Udata[i]->e[j][k] = r->U->e[i][k];
        }
        k0 += r->n_controls;
    }
}

void _setUdata2(RKOCP r, Vector x) {
    double *U0, *U1, tau, s0;//, s1;
    int n_parameters = r->n_parameters;
    int n_controls = r->n_controls;
    int n_states = r->n_states;
    int rk_stages = r->rk_stages;
    int n_nodes = r->n_nodes;
    double *rk_c = r->rk_c->e;
    int i_node, i, j, k;
    int k0 = n_parameters + n_states;
    for (k = 0; k < n_controls; k++) r->U->e[0][k] = x->e[k0 + k];
    k0 += n_controls;
    for (i_node = 0; i_node < (n_nodes - 1); i_node++) {
        U0 = r->U->e[i_node];
        U1 = r->U->e[i_node + 1];
        for (k = 0; k < n_controls; k++) r->U->e[i_node + 1][k] = x->e[k0 + k];
        k0 += n_controls;
        for (j = 0; j < rk_stages; j++) {
            tau = rk_c[j];
            s0 = 1.0 - tau;
            for (i = 0; i < n_controls; i++) {
                r->Udata[i_node]->e[j][i] = s0 * U0[i] + tau * U1[i];
            }
        }
    }
}

void _setUdata8(RKOCP r, Vector x) {
    int n_parameters = r->n_parameters;
    int n_controls = r->n_controls;
    int n_states = r->n_states;
    int rk_stages = r->rk_stages;
    int n_nodes = r->n_nodes;
    int k0 = n_parameters + n_states;
    int i, j, k, i_node;
    for (i_node = 0; i_node < n_nodes; i_node++) {
        for (k = 0; k < n_controls; k++) r->U->e[i_node][k] = x->e[k0 + k];
        k0 += n_controls;
    }
    double *Us[6];
    for (i_node = 0; i_node < (n_nodes - 1); i_node++) {
        if ((i_node == 0) || (i_node == 1)) {
            for (i = 0; i < 6; i++) {
                Us[i] = r->U->e[i];
            }
        } else if ((i_node == (n_nodes - 3)) || (i_node == (n_nodes - 2))) {
            for (i = 0; i < 6; i++) {
                Us[i] = r->U->e[n_nodes - 6 + i];
            }
        } else {
            for (i = 0; i < 6; i++) {
                Us[i] = r->U->e[i_node - 2 + i];
            }
        }
        MatrixSetAllTo(r->Udata[i_node], 0.0);
        for (j = 0; j < rk_stages; j++) {
            for (k = 0; k < 6; k++) {
                _axpy(r->Udata[i_node]->e[j], r->Ucubic[i_node]->e[j][k], Us[k], n_controls);
            }
        }
    }
}

// y = a * x + y
void _axpy(double *y, double a, double *x, int n) {
    int i;
    for (i = 0; i < n; i++) {
        y[i] += a * x[i];
    }
}

void _addDphi(RKOCP r, Matrix Yxi, Vector nlp_df) {
    int i, j;
    double s;
    int n_states = r->n_states;
    int n_parameters = r->n_parameters;
    int n = r->nlp_n;
    for (i = 0; i < n; i++) {
        s = 0;
        for (j = 0; j < n_states; j++) {
            s += r->phiy->e[j] * Yxi->e[j][i];
        }
        nlp_df->e[i] += s;
    }

    for (i = 0; i < n_parameters; i++) {
        nlp_df->e[i] += r->phip->e[i];
    }
}

void _setDhGamma(RKOCP r, Vector ys, Vector ps, Matrix nlp_dh) {
    // get Gamma_y and Gamma_p
    int i, j;
    int n_initial = r->n_initial;
    int n_parameters = r->n_parameters;
    int n_states = r->n_states;

    if (r->ocp->Dinitial_constraints == NULL) {
        OCPFDGamma(r->ocp, ys, ps, r->Gammay, r->Gammap);
    } else {
        r->ocp->Dinitial_constraints(ys, ps, r->Gammay, r->Gammap);
    }
    for (i = 0; i < n_initial; i++) {
        for (j = 0; j < n_parameters; j++) {
            nlp_dh->e[i][j] = r->Gammap->e[i][j];
        }
        for (j = 0; j < n_states; j++) {
            nlp_dh->e[i][j + n_parameters] = r->Gammay->e[i][j];
        }
    }
}

void _setDhPsi(RKOCP r, Matrix Yxi, Matrix nlp_dh) {
    int n_initial = r->n_initial;
    int n_terminal = r->n_terminal;
    int n_parameters = r->n_parameters;
    int n_states = r->n_states;
    int n = r->nlp_n;
    int i, j, k;
    double s;
    for (i = 0; i < n_terminal; i++) {
        for (j = 0; j < n; j++) {
            s = 0.0;
            for (k = 0; k < n_states; k++) {
                s += r->Psiy->e[i][k] * Yxi->e[k][j];
            }
            nlp_dh->e[i + n_initial][j] += s;
        }
    }
    if (n_parameters > 0) {
        for (i = 0; i < n_terminal; i++) {
            for (j = 0; j < n_parameters; j++) {
                nlp_dh->e[i + n_initial][j] += r->Psip->e[i][j];
            }
        }
    }
}

void _setDg(RKOCP r, Matrix Yxi, Vector ys, Vector us,
        Vector ps, double t, int i_node, Matrix nlp_dg) {
    // Dg = g_y * Yx + g_u * Ux + g_p * Px
    int i, j, k, row, col0, col1;
    double s;
    int n_inequality = r->n_inequality;
    int n_controls = r->n_controls;
    int n_parameters = r->n_parameters;
    int n_states = r->n_states;
    int n_nodes = r->n_nodes;

    // get g_y, g_u, g_p
    if (r->ocp->Dinequality_constraints == NULL) {
        OCPFDd(r->ocp, ys, us, ps, t, r->Gy, r->Gu, r->Gp);
    } else {
        r->ocp->Dinequality_constraints(ys, us, ps, t, r->Gy, r->Gu, r->Gp);
    }
    row = i_node * n_inequality;

    // g_y * Yx
    col1 = 0;
    switch (r->c_type) {
        case CONSTANT:
            col1 = n_parameters + n_states + i_node * n_controls;
            break;
        case LINEAR:
            col1 = n_parameters + n_states + (i_node + 1) * n_controls;
            break;
        case CUBIC:
            col1 = n_parameters + n_states + (i_node + 1) * n_controls;
            break;
    }
    for (i = 0; i < n_inequality; i++) {
        for (j = 0; j < col1; j++) {
            s = 0;
            for (k = 0; k < n_states; k++) {
                s += r->Gy->e[i][k] * Yxi->e[k][j];
            }
            nlp_dg->e[i + row][j] += s;
        }
    }
    // g_u * Ux
    if (n_controls > 0) {
        col0 = 0;
        int ind;
        switch (r->c_type) {
            case CONSTANT:
                ind = i_node > (n_nodes - 2) ? (n_nodes - 2) : i_node;
                col0 = n_parameters + n_states + ind * n_controls;
                break;
            case LINEAR:
                col0 = n_parameters + n_states + i_node * n_controls;
                break;
            case CUBIC:
                col0 = n_parameters + n_states + i_node * n_controls;
                break;
        }
        for (i = 0; i < n_inequality; i++) {
            for (j = 0; j < n_controls; j++) {
                nlp_dg->e[i + row][j + col0] += r->Gu->e[i][j];
            }
        }
    }
    // g_p * Px
    if (n_parameters > 0) {
        for (i = 0; i < n_inequality; i++) {
            for (j = 0; j < n_parameters; j++) {
                nlp_dg->e[i + row][j] += r->Gp->e[i][j];
            }
        }
    }
}

void _setYx_stage(RKOCP r, Matrix Yx_stage_l, Matrix Yxi,
        int i_node, int l, double delta, Matrix* Kx_stage) {
    int rr;
    // Yx_stage[l] = Yx[i_node] + delta_{i_node} *
    //      Sum_{r = 0; to l-1} a(i,r) * Kx_stage[r]
    MatrixSetEqual(Yx_stage_l, Yxi);
    for (rr = 0; rr < l; rr++) {
        _Mapsx(r, Yx_stage_l, delta * r->rk_a->e[l][rr], Kx_stage[rr], i_node);
    }
}

void _setYx(RKOCP r, Matrix Yxi, Matrix *Kstagei, int i_node,
        double delta) {
    int rr;
    // Y += delta * b[r] * Kx_stage
    for (rr = 0; rr < r->rk_stages; rr++) {
        _Mapsx(r, Yxi, delta * r->rk_b->e[rr], Kstagei[rr], i_node);
    }
}

void _Mapsx(RKOCP r, Matrix Yxi, double s, Matrix B,
        int i_node) {
    int n_controls = r->n_controls;
    int n_parameters = r->n_parameters;
    int n_nodes = r->n_nodes;
    int n_states = r->n_states;
    int c0, c1, i, j;
    int col = 0;
    double *_Yxii, *_Bi;
    switch (r->c_type) {
        case CONSTANT:
            col = n_parameters + n_states + (i_node + 1) * n_controls;
            break;
        case LINEAR:
            col = n_parameters + n_states + (i_node + 2) * n_controls;
            break;
        case CUBIC:
            c0 = (i_node + 6) * n_controls;
            c1 = n_nodes * n_controls;
            c0 = c0 < c1 ? c0 : c1;
            col = n_parameters + n_states + c0;
            break;
    }
    for (i = 0; i < Yxi->r; i++) {
        _Yxii = Yxi->e[i];
        _Bi = B->e[i];
        for (j = 0; j < col; j++) {
            _Yxii[j] += s * _Bi[j];
        }
    }
}

void _setKx_stage(RKOCP r, Matrix Kx_stage_l, Matrix Yx_stage_l, int i_node, int l, Vector ys, Vector us, Vector ps, double t) {
    int *kn8 = NULL;
    int kn80[] = {0, 1, 2, 3, 4};
    int kn8j[] = {0, 1, 2, 3, 4, 5};
    int kn8n[] = {1, 2, 3, 4, 5};
    int knlength; // 80: 5, 8j : 6, 8n : 5
    int n_parameters = r->n_parameters;
    int n_states = r->n_states;
    int n_controls = r->n_controls;
    int n_nodes = r->n_nodes;
    // K = [Fy * Yx + Fu * Ux + Fp * Px; Ly * Yx + Lu * Ux + Lp * Px]
    int i, j, k, col, c0, c1, ir;
    double s;
    double *_Fyi, *_Kxi, *_Fui, *_Fpi;

    // get Fy, Fu, Fp
    //ocp.Dstates(ys, us, ps, t, Fy, Fu, Fp, Ly, Lu, Lp);
    if (r->ocp->Ddifferential_equations == NULL) {
        OCPFDfL(r->ocp, ys, us, ps, t, r->Fy, r->Fu, r->Fp, r->Ly, r->Lu, r->Lp);
    } else {
        r->ocp->Ddifferential_equations(ys, us, ps, t, r->Fy, r->Fu, r->Fp, r->Ly, r->Lu, r->Lp);
    }

    MatrixSetAllTo(Kx_stage_l, 0.0);

    //  [Fy * Yx; Ly * Yx]
    col = 0;
    switch (r->c_type) {
        case CONSTANT:
            col = n_parameters + n_states + (i_node + 1) * n_controls;
            break;
        case LINEAR:
            col = n_parameters + n_states + (i_node + 2) * n_controls;
            break;
        case CUBIC:
            c0 = (i_node + 6) * n_controls;
            c1 = n_nodes * n_controls;
            c0 = c0 < c1 ? c0 : c1;
            col = n_parameters + n_states + c0;
            break;
    }

    for (j = 0; j < col; j++) {
        for (i = 0; i < n_states; i++) {
            s = 0;
            _Fyi = r->Fy->e[i];
            for (k = 0; k < n_states; k++) {
                s += _Fyi[k] * Yx_stage_l->e[k][j];
            }
            Kx_stage_l->e[i][j] += s;
        }
        s = 0;
        for (k = 0; k < n_states; k++) {
            s += r->Ly->e[k] * Yx_stage_l->e[k][j];
        }
        Kx_stage_l->e[n_states][j] += s;
    }
    
    // [Fu; Lu] * Ux
    if (n_controls > 0) {
        switch (r->c_type) {
            case CONSTANT:
                col = n_parameters + n_states + i_node * n_controls;
                for (i = 0; i < n_states; i++) {
                    _Kxi = Kx_stage_l->e[i];
                    _Fui = r->Fu->e[i];
                    for (j = 0; j < n_controls; j++) {
                        _Kxi[j + col] += _Fui[j];
                    }
                }
                for (j = 0; j < n_controls; j++) {
                    Kx_stage_l->e[n_states][j + col] += r->Lu->e[j];
                }
                break;
            case LINEAR:
                col = n_parameters + n_states + i_node * n_controls;
                s = 1.0 - r->rk_c->e[l];
                for (i = 0; i < n_states; i++) {
                    _Fui = r->Fu->e[i];
                    _Kxi = Kx_stage_l->e[i];
                    for (j = 0; j < n_controls; j++) {
                        _Kxi[j + col] += s * _Fui[j];
                    }
                }
                for (j = 0; j < n_controls; j++) {
                    Kx_stage_l->e[n_states][j + col] += s * r->Lu->e[j];
                }
                col += n_controls;
                s = r->rk_c->e[l];
                for (i = 0; i < n_states; i++) {
                    _Fui = r->Fu->e[i];
                    _Kxi = Kx_stage_l->e[i];
                    for (j = 0; j < n_controls; j++) {
                        _Kxi[j + col] += s * _Fui[j];
                    }
                }
                for (j = 0; j < n_controls; j++) {
                    Kx_stage_l->e[n_states][j + col] += s * r->Lu->e[j];
                }
                break;
            case CUBIC:
            {
                if ((i_node == 0) || (i_node == 1)) {
                    // t0 t1 t2 t3 t4
                    kn8 = kn80;
                    knlength = 5; // 80: 5, 8j : 6, 8n : 5
                    col = n_parameters + n_states;
                } else if (i_node > (n_nodes - 4)) {
                    // i_node == n_nodes - 2
                    // i_node == n_nodes - 3
                    // tN-5 tN-4 tN-3 tN-2 tN-1
                    kn8 = kn8n;
                    knlength = 5;
                    col = n_parameters + n_states + (n_nodes - 5) * n_controls;
                } else {
                    // tj-2 tj-1 tj tj+1 tj+2 tj+3
                    kn8 = kn8j;
                    knlength = 6;
                    col = n_parameters + n_states + (i_node - 2) * n_controls;
                }
                for (ir = 0; ir < knlength; ir++) {
                    s = r->Ucubic[i_node]->e[l][kn8[ir]];
                    for (i = 0; i < n_states; i++) {
                        _Fui = r->Fu->e[i];
                        _Kxi = Kx_stage_l->e[i];
                        for (j = 0; j < n_controls; j++) {
                            _Kxi[j + col] += s * _Fui[j];
                        }
                    }
                    for (j = 0; j < n_controls; j++) {
                        Kx_stage_l->e[n_states][j + col] += s * r->Lu->e[j];
                    }
                    col += n_controls;
                }
            }
            break;
        }
    }

    // Fp * Px
    if (n_parameters > 0) {
        for (i = 0; i < n_states; i++) {
            _Fpi = r->Fp->e[i];
            _Kxi = Kx_stage_l->e[i];
            for (j = 0; j < n_parameters; j++) {
                _Kxi[j] += _Fpi[j];
            }
        }
        for (j = 0; j < n_parameters; j++) {
            Kx_stage_l->e[n_states][j] += r->Lp->e[j];
        }
    }
}

