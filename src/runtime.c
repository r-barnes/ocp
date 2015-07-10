// Copyright (c) 2009,2010,2011,2012 Brian C. Fabien
// All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the license.txt file.

#include <stdio.h>
#include <stdlib.h>
#include <time.h>

#include "runtime.h"

static int first_call = 0;
/*
double fmax(double a, double b) {
    return a > b ? a : b;
}

double fmin(double a, double b) {
    return a < b ? a : b;
}
*/
double sign(double x) {
    if (x < 0.0) {
        return -1.0;
    }
    
    return 1.0;
}

double drand() {
    double d, rmax_2;
    //time_t tp;
    //struct tm *lt;
    unsigned int seed;

    if (first_call == 0) {
        //tp = time((time_t *) NULL);
        //lt = localtime(&tp);
        //seed = abs(lt->tm_sec + lt->tm_min + lt->tm_hour + lt->tm_mday);
        //seed = 123456;
        seed = 1;
        srand(seed);
        first_call = 1;
    }
    rmax_2 = ((double) RAND_MAX / 2);
    d = ((double) rand() - rmax_2) / rmax_2;
    return d;
}

int irand(int M) {
    //double d = (drand() + ((double) RAND_MAX / 2))*M;
    //double d = drand()*((double)M);
    double d = (((double) rand()) / ((double)RAND_MAX + 1.0))*((double)M); // must call drand() first
    return abs((int)d + 1);
}

void displayTime() {
    time_t tp;
    struct tm *lt;
    tp = time((time_t *) NULL);
    lt = localtime(&tp);
    printf("%d:%d:%d\n", lt->tm_hour, lt->tm_min, lt->tm_sec);
}

void pause() {
    char c;
    printf("\nPause ...");
    c = getchar();
    c = c;
}

void RuntimeError(const char *string) {
	printf("Fatal error: %s\n", string);
	exit(1);
}

void RuntimeWarning(const char *string) {
	printf("Warning: %s\n", string);
}
