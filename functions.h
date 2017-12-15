#ifndef FUNCTIONS_H
#define FUNCTIONS_H

#include "matrix.h"
#include "constants.h"
#include "lu.h"

void input(double &del_x, double &del_t);

void dufort(double dx, double dt, int n);

void richardson(double dx, double dt, int n);

void laasonen(double delx, double delt, double n);

void crank(double delx, double delt, double n);

void crank1(double delx, double delt, double n);

void analytical(double n);

void write_file(Vector V, const char *filename);

#endif
