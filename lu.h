#ifndef LU_H
#define LU_H

#include "matrix.h"

void lu_fact(const Matrix& a, Matrix& l, Matrix& u, int n);
void lu_solve(const Matrix& l, const Matrix& u, const Vector& b, int n, Vector& x);
void reorder(const Matrix& a, int n, Matrix& p);

#endif
