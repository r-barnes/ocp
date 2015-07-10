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
#include "equations.h"
#include "nlp.h"
#include "ocp.h"

static OCP thisOCP;

void _OCPallocateStorage(OCP o);
void _OCPFGamma(Vector x, Vector F);
void _OCPFphiPsi(Vector x, Vector F);
void _OCPFfL(Vector x, Vector F);
void _OCPFd(Vector x, Vector F);

// create a new OCP data structure
//
// Public members:
//
//    int n_states;       // Number of states > 0  
//    int n_controls;     // Number of controls >= 0
//    int n_parameters;   // Number of parameters >= 0 
//    int n_initial;      // Number of initial conditions >= 0 
//    int n_terminal;     // Number of terminal conditions >= 0 
//    int n_inequality;   // Number of inequality constraints >= 0 
//    void (*initial_constraints)(Vector y, Vector p, Vector Gamma);  // Computer Gamma(y,p)
//    double (*terminal_constraints)(Vector y, Vector p, Vector Psi); // Compute phi, Psi
//    void (*inequality_constraints)(Vector y, Vector u, Vector p,
//        double t, Vector d);                                        // Compute d
//    double (*differential_equations)(Vector y, Vector u, Vector p,
//        double t, Vector F);                                        // Compute L, F
//    void (*Dinitial_constraints)(Vector y, Vector p, Matrix Gammay,
//        Matrix Gammap);                                             // derivatives of Gamma
//    void (*Dterminal_constraints)(Vector y, Vector p,
//            Vector phiy, Vector phip, Matrix Psiy, Matrix Psip);    // derivatives of phi and Psi
//    void (*Dinequality_constraints)(Vector y, Vector u, Vector p,
//        double t, Matrix dy, Matrix du, Matrix dp);                 // derivatives of d
//    void (*Ddifferential_equations)(Vector y, Vector u, Vector p,
//        double t, Matrix Fy, Matrix Fu, Matrix Fp, Vector Ly,
//        Vector Lu, Vector Lp);                                      // derivatives of L and F
OCP OCPNew() {
    OCP o;
    o = (OCP) malloc(sizeof (struct ocp_));
    o->n_states = 0;
    o->n_controls = 0;
    o->n_parameters = 0;
    o->n_initial = 0;
    o->n_terminal = 0;
    o->n_inequality = 0;
    o->initial_constraints = NULL;
    o->differential_equations = NULL;
    o->inequality_constraints = NULL;
    o->terminal_constraints = NULL;
    o->Ddifferential_equations = NULL;
    o->Dinequality_constraints = NULL;
    o->Dinitial_constraints = NULL;
    o->Dterminal_constraints = NULL;

    // private
    o->_t = 0.0;

    o->_y = NULL;
    o->_u = NULL;
    o->_p = NULL;
    o->_yp = NULL;
    o->_yup = NULL;

    o->_G = NULL; //Equations
    o->_dgamma = NULL; //Matrix

    o->_phiPsi = NULL; //Equations
    o->_dphiPsi = NULL; //Matrix

    o->_fL = NULL; //Equations
    o->_dfL = NULL; //Matrix

    o->_d = NULL; //Equations
    o->_dd = NULL; //Matrix

    o->_memory_initilized = NO;

    return o;
}

// delete the OCP data structure
void OCPDelete(OCP o) {
    if (o->_memory_initilized == YES) {
        VectorDelete(o->_y);
        VectorDelete(o->_yp);
        VectorDelete(o->_yup);
        if (o->n_controls > 0) {
            VectorDelete(o->_u);
        }
        if (o->n_parameters > 0) {
            VectorDelete(o->_p);
        }
        EquationsDelete(o->_fL);
        MatrixDelete(o->_dfL);
        if (o->n_initial > 0) {
            EquationsDelete(o->_G);
            MatrixDelete(o->_dgamma);
        }
        EquationsDelete(o->_phiPsi);
        MatrixDelete(o->_dphiPsi);

        if (o->n_inequality > 0) {
            EquationsDelete(o->_d);
            MatrixDelete(o->_dd);
        }
    }
    free(o);
}

void _OCPallocateStorage(OCP o) {
    int n_states = o->n_states;
    int n_parameters = o->n_parameters;
    int n_initial = o->n_initial;
    int n_terminal = o->n_terminal;
    int n_inequality = o->n_inequality;
    int n_controls = o->n_controls;
    // _y, _u, _p, _yp, _yup
    if (o->n_states < 1) {
        RuntimeError("OCPallocateStorage: n_states < 1");
    }
    o->_y = VectorNew(n_states);
    o->_yp = VectorNew(n_states + n_parameters);
    o->_yup = VectorNew(n_states + n_controls + n_parameters);
    if (n_controls > 0) {
        o->_u = VectorNew(n_controls);
    }
    if (n_parameters > 0) {
        o->_p = VectorNew(n_parameters);
    }
    if (n_initial > 0) {
        o->_G = EquationsNew();
        o->_G->F = _OCPFGamma;
        o->_dgamma = MatrixNew(n_initial, n_states + n_parameters);
    }
    o->_fL = EquationsNew();
    o->_fL->F = _OCPFfL;
    o->_dfL = MatrixNew(1 + n_states, n_states + n_controls + n_parameters);
    if (n_inequality > 0) {
        o->_d = EquationsNew();
        o->_d->F = _OCPFd;
        o->_dd = MatrixNew(n_inequality, n_states + n_controls + n_parameters);
    }
    o->_phiPsi = EquationsNew();
    o->_phiPsi->F = _OCPFphiPsi;
    o->_dphiPsi = MatrixNew(1 + n_terminal, n_states + n_parameters);
    o->_memory_initilized = YES;
}

// Finite difference Gamma_y, Gamma_p
void _OCPFGamma(Vector x, Vector F) {
    OCP o = thisOCP;
    int n_states = o->n_states;
    int n_parameters = o->n_parameters;
    int i;
    for (i = 0; i < n_states; i++) {
        o->_y->e[i] = x->e[i];
    }
    for (i = 0; i < n_parameters; i++) {
        o->_p->e[i] = x->e[n_states + i];
    }
    o->initial_constraints(o->_y, o->_p, F);
}

void OCPFDGamma(OCP o, Vector y, Vector p, Matrix Gy, Matrix Gp) {
    int n_states = o->n_states;
    int n_parameters = o->n_parameters;
    int n_initial = o->n_initial;
    int i, j;
    thisOCP = o;
    if (o->_memory_initilized == NO) {
        _OCPallocateStorage(o);
    }
    for (i = 0; i < n_states; i++) {
        o->_yp->e[i] = y->e[i];
    }
    for (i = 0; i < n_parameters; i++) {
        o->_yp->e[n_states + i] = p->e[i];
    }

    EquationsDF(o->_G, o->_yp, o->_dgamma);

    for (i = 0; i < n_initial; i++) {
        for (j = 0; j < n_states; j++) {
            Gy->e[i][j] = o->_dgamma->e[i][j];
        }
        for (j = 0; j < n_parameters; j++) {
            Gp->e[i][j] = o->_dgamma->e[i][n_states + j];
        }
    }
}

// Finite difference phi_y, phi_p, Psi_y, Psi_p
void _OCPFphiPsi(Vector x, Vector F) {
    OCP o = thisOCP;
    int n_states = o->n_states;
    int n_parameters = o->n_parameters;
    int n_terminal = o->n_terminal;
    int i;
    for (i = 0; i < n_states; i++) {
        o->_y->e[i] = x->e[i];
    }
    for (i = 0; i < n_parameters; i++) {
        o->_p->e[i] = x->e[n_states + i];
    }
    F->e[n_terminal] = o->terminal_constraints(o->_y, o->_p, F);
}

void OCPFDphiPsi(OCP o, Vector y, Vector p, Vector phiy, Vector phip,
        Matrix Psiy, Matrix Psip) {
    int n_states = o->n_states;
    int n_parameters = o->n_parameters;
    int n_terminal = o->n_terminal;
    int i, j;
    thisOCP = o;
    if (o->_memory_initilized == NO) {
        _OCPallocateStorage(o);
    }
    for (i = 0; i < n_states; i++) {
        o->_yp->e[i] = y->e[i];
    }
    for (i = 0; i < n_parameters; i++) {
        o->_yp->e[n_states + i] = p->e[i];
    }

    EquationsDF(o->_phiPsi, o->_yp, o->_dphiPsi);

    for (i = 0; i < n_terminal; i++) {
        for (j = 0; j < n_states; j++) {
            Psiy->e[i][j] = o->_dphiPsi->e[i][j];
        }
        for (j = 0; j < n_parameters; j++) {
            Psip->e[i][j] = o->_dphiPsi->e[i][n_states + j];
        }
    }
    for (i = 0; i < n_states; i++) {
        phiy->e[i] = o->_dphiPsi->e[n_terminal][i];
    }
    for (i = 0; i < n_parameters; i++) {
        phip->e[i] = o->_dphiPsi->e[n_terminal][n_states + i];
    }
}

// Finite difference L_y, L_u, L_p, F_y, F_u, F_p
void _OCPFfL(Vector x, Vector F) {
    OCP o = thisOCP;
    int n_states = o->n_states;
    int n_parameters = o->n_parameters;
    int n_controls = o->n_controls;
    int i;
    for (i = 0; i < n_states; i++) {
        o->_y->e[i] = x->e[i];
    }
    for (i = 0; i < n_controls; i++) {
        o->_u->e[i] = x->e[n_states + i];
    }
    for (i = 0; i < n_parameters; i++) {
        o->_p->e[i] = x->e[n_states + n_controls + i];
    }
    F->e[n_states] = o->differential_equations(o->_y, o->_u, o->_p, o->_t, F);
}

void OCPFDfL(OCP o, Vector y, Vector u, Vector p, double t,
        Matrix fy, Matrix fu, Matrix fp,
        Vector Ly, Vector Lu, Vector Lp) {
    int n_states = o->n_states;
    int n_parameters = o->n_parameters;
    int n_controls = o->n_controls;
    int i, j;

    thisOCP = o;

    if (o->_memory_initilized == NO) {
        _OCPallocateStorage(o);
    }

    o->_t = t;

    for (i = 0; i < n_states; i++) {
        o->_yup->e[i] = y->e[i];
    }
    for (i = 0; i < n_controls; i++) {
        o->_yup->e[n_states + i] = u->e[i];
    }
    for (i = 0; i < n_parameters; i++) {
        o->_yup->e[n_states + n_controls + i] = p->e[i];
    }

    EquationsDF(o->_fL, o->_yup, o->_dfL);

    for (i = 0; i < n_states; i++) {
        for (j = 0; j < n_states; j++) {
            fy->e[i][j] = o->_dfL->e[i][j];
        }

        for (j = 0; j < n_controls; j++) {
            fu->e[i][j] = o->_dfL->e[i][n_states + j];
        }
        for (j = 0; j < n_parameters; j++) {
            fp->e[i][j] = o->_dfL->e[i][n_states + n_controls + j];
        }
    }
    for (i = 0; i < n_states; i++) {
        Ly->e[i] = o->_dfL->e[n_states][i];
    }

    for (i = 0; i < n_controls; i++) {
        Lu->e[i] = o->_dfL->e[n_states][n_states + i];
    }
    for (i = 0; i < n_parameters; i++) {
        Lp->e[i] = o->_dfL->e[n_states][n_states + n_controls + i];
    }
}

// Finite difference d_y, d_u, d_p
void _OCPFd(Vector x, Vector F) {
    OCP o = thisOCP;
    int n_states = o->n_states;
    int n_parameters = o->n_parameters;
    int n_controls = o->n_controls;
    int i;
    for (i = 0; i < n_states; i++) {
        o->_y->e[i] = x->e[i];
    }
    for (i = 0; i < n_controls; i++) {
        o->_u->e[i] = x->e[n_states + i];
    }
    for (i = 0; i < n_parameters; i++) {
        o->_p->e[i] = x->e[n_states + n_controls + i];
    }
    o->inequality_constraints(o->_y, o->_u, o->_p, o->_t, F);
}

void OCPFDd(OCP o, Vector y, Vector u, Vector p, double t,
        Matrix dy, Matrix du, Matrix dp) {
    int n_states = o->n_states;
    int n_parameters = o->n_parameters;
    int n_controls = o->n_controls;
    int n_inequality = o->n_inequality;
    int i, j;

    thisOCP = o;

    if (o->_memory_initilized == NO) {
        _OCPallocateStorage(o);
    }

    o->_t = t;

    for (i = 0; i < n_states; i++) {
        o->_yup->e[i] = y->e[i];
    }
    for (i = 0; i < n_controls; i++) {
        o->_yup->e[n_states + i] = u->e[i];
    }
    for (i = 0; i < n_parameters; i++) {
        o->_yup->e[n_states + n_controls + i] = p->e[i];
    }

    EquationsDF(o->_d, o->_yup, o->_dd);

    for (i = 0; i < n_inequality; i++) {
        for (j = 0; j < n_states; j++) {
            dy->e[i][j] = o->_dd->e[i][j];
        }
        for (j = 0; j < n_controls; j++) {
            du->e[i][j] = o->_dd->e[i][n_states + j];
        }
        for (j = 0; j < n_parameters; j++) {
            dp->e[i][j] = o->_dd->e[i][n_states + n_controls + j];
        }
    }
}

// Check the validity of the data structure
// return 0 if OK, -1 otherwise
int OCPCheck(OCP o) {
    if (o == NULL) {
        RuntimeWarning("OCPcheck: ocp == NULL");
        return -1;
    }
    if (o->n_states < 1) {
        RuntimeWarning("OCPcheck: ocp->n_states < 1");
        return -1;
    }
    if (o->n_controls < 0) {
        RuntimeWarning("OCPcheck: ocp->n_controls < 0");
        return -1;
    }
    if (o->n_parameters < 0) {
        RuntimeWarning("OCPcheck: ocp->n_controls < 0");
        return -1;
    }
    if (o->n_initial < 0) {
        RuntimeWarning("OCPcheck: ocp->n_initial < 0");
        return -1;
    }
    if (o->n_terminal < 0) {
        RuntimeWarning("OCPcheck: ocp->n_terminal < 0");
        return -1;
    }
    if (o->n_inequality < 0) {
        RuntimeWarning("OCPcheck: ocp->n_inequality < 0");
        return -1;
    }
    if (o->differential_equations == NULL) {
        RuntimeWarning("OCPcheck: ocp->differential_equations == NULL");
        return -1;
    }
    if ((o->n_initial > 0) && (o->initial_constraints == NULL)) {
        RuntimeWarning("OCPcheck: ocp->initial_constraints == NULL");
        return -1;
    }
    if ((o->n_terminal > 0) && (o->terminal_constraints == NULL)) {
        RuntimeWarning("OCPcheck: ocp->terminal_constraints == NULL");
        return -1;
    }
    if ((o->n_inequality > 0) && (o->inequality_constraints == NULL)) {
        RuntimeWarning("OCPcheck: ocp->inequality_constraints == NULL");
        return -1;
    }
    return 0;
}

int OCPWriteData(OCP o, FILE *file_id,
        Vector T, Matrix Y, Matrix U, Vector P,
        int control_type, double tol) {
    int i, j;
    int n_nodes;
    if (o == NULL) {
        RuntimeWarning("OCPWriteData: o == NULL\n");
        return -1;
    }
    if (T == NULL) {
        RuntimeWarning("OCPWriteData: T == NULL\n");
        return -1;
    }
    if (Y == NULL) {
        RuntimeWarning("OCPWriteData: Y == NULL\n");
        return -1;
    }
    if (file_id == NULL) {
        RuntimeWarning("OCPWriteData: file_id == NULL\n");
        return -1;
    }
    n_nodes = T->r;
    fprintf(file_id, "n_states\n%d\n", o->n_states);
    fprintf(file_id, "n_controls\n%d\n", o->n_controls);
    fprintf(file_id, "n_parameters\n%d\n", o->n_parameters);
    fprintf(file_id, "n_initial\n%d\n", o->n_initial);
    fprintf(file_id, "n_terminal\n%d\n", o->n_terminal);
    fprintf(file_id, "n_inequality\n%d\n", o->n_inequality);
    fprintf(file_id, "n_nodes\n%d\n", n_nodes);
    fprintf(file_id, "tolerance\n%e\n", tol);
    fprintf(file_id, "control_type\n%d\n", control_type);
    fprintf(file_id, "T\n");
    for (i = 0; i < n_nodes; i++)
        fprintf(file_id, "%g ", T->e[i]);
    fprintf(file_id, "\n");
    fprintf(file_id, "Y\n");
    for (i = 0; i < o->n_states; i++) {
        for (j = 0; j < n_nodes; j++)
            fprintf(file_id, "%e ", Y->e[j][i]);
        fprintf(file_id, "\n");
    }
    //fprintf(file_id, "\n");
    fprintf(file_id, "U\n");
    for (i = 0; i < o->n_controls; i++) {
        for (j = 0; j < n_nodes; j++)
            fprintf(file_id, "%e ", U->e[j][i]);
        fprintf(file_id, "\n");
    }
    //fprintf(file_id, "\n");
    fprintf(file_id, "P\n");
    for (i = 0; i < o->n_parameters; i++)
        fprintf(file_id, "%e ", P->e[i]);
    fprintf(file_id, "\n");

    fclose(file_id);

    return 0;
}

int OCPReadData(OCP o, FILE *file_id, Vector T, Matrix Y, Matrix U, Vector P) {
    int i, j, err;
    int n_states, n_nodes, n_controls, n_parameters;
    //int ignore;
    double x;
    char s[256];

    if (file_id == NULL) {
    	RuntimeWarning("OCPReadData(): file_id == NULL\n");
        return -1;
    }
    err = fscanf(file_id, "%s", s);
    err = fscanf(file_id, "%d", &i); //printf("%s\n%d\n", s,i);
    n_states = i;
    err = fscanf(file_id, "%s", s);
    err = fscanf(file_id, "%d", &i); //printf("%s\n%d\n", s,i);
    n_controls = i;
    err = fscanf(file_id, "%s", s);
    err = fscanf(file_id, "%d", &i); //printf("%s\n%d\n", s,i);
    n_parameters = i;
    err = fscanf(file_id, "%s", s);
    err = fscanf(file_id, "%d", &i); //printf("%s\n%d\n", s,i);
    //o->n_initial = i;
    err = fscanf(file_id, "%s", s);
    err = fscanf(file_id, "%d", &i); //printf("%s\n%d\n", s,i);
    //o->n_terminal = i;
    err = fscanf(file_id, "%s", s);
    err = fscanf(file_id, "%d", &i); //printf("%s\n%d\n", s,i);
    //o->n_inequality = i;
    err = fscanf(file_id, "%s", s);
    err = fscanf(file_id, "%d", &i); //printf("%s\n%d\n", s,i);
    n_nodes = i;
    err = fscanf(file_id, "%s", s);
    err = fscanf(file_id, "%le", &x); //printf("%s\n%g\n", s,x);
    //p->data->tolerance = x;
    err = fscanf(file_id, "%s", s);
    err = fscanf(file_id, "%d", &i); //printf("%s\n%d\n", s,i); pause();
    //c_type = i;

    //check T, Y, U, P
    if ((n_nodes != T->r) || (n_nodes != Y->r)) {
        RuntimeWarning("OCPReadData: the number of nodes in the data file does not match the dimensions of T and Y");
        return -1;
    }
    if (n_states != Y->c) {
        RuntimeWarning("OCPReadData: the number of states in the data file does not match Y");
        return -1;
    }
    if (U != NULL) {
        if ((n_controls != U->c) || (n_nodes != U->r)) {
            RuntimeWarning("OCPReadData: the number of nodes and controls in the data file does not match U");
            return -1;
        }
    }
    if ((U == NULL) && (n_controls > 0)) {
        RuntimeWarning("OCPReadData: the number of controls in the data file does not match U");
        return -1;
    }
    if (P != NULL) {
        if (n_parameters != P->r) {
            RuntimeWarning("OCPReadData: the number of parameters in the data file does not match P");
            n_parameters = 0;
        }
    }
    err = fscanf(file_id, "%s", s); /* T */
    for (i = 0; i < n_nodes; i++) {
        err = fscanf(file_id, "%le", &x);
        T->e[i] = x;
    }
    err = fscanf(file_id, "%s", s); /* Y */
    for (i = 0; i < n_states; i++)
        for (j = 0; j < n_nodes; j++) {
            err = fscanf(file_id, "%le", &x);
            Y->e[j][i] = x;
        }
    err = fscanf(file_id, "%s", s); /* U */
    for (i = 0; i < n_controls; i++)
        for (j = 0; j < n_nodes; j++) {
            err = fscanf(file_id, "%le", &x);
            U->e[j][i] = x;
        }
    err = fscanf(file_id, "%s", s); /* P */
    for (i = 0; i < n_parameters; i++) {
        err = fscanf(file_id, "%le", &x);
        P->e[i] = x;
    }
    fclose(file_id);

    return 0*err; // surpress warning that err is not used
}
