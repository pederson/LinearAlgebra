#ifndef _LINALGSOLVERS_H
#define _LINALGSOLVERS_H

// collection of solvers for linear algebra problems

// upper triangular matrix Ux = b
Vector & upper_triangular_solve(const Matrix & U, const Vector & b);

// lower triangular matrix Lx = b
Vector & lower_triangular_solve(const Matrix & L, const Vector & b);

// diagonal matrix Dx = b
Vector & diagonal_solve(const Matrix & D, const Vector & b);

// unitary matrix Ux = b
Vector & unitary_solve(const Matrix & U, const Vector & b);

// tridiagonal matrix Tx = b
Vector & tridiagonal_solve(const Matrix & T, const Vector & b);





// collection of helper routines

// qr factorization using Gram-Schmidt algorithm A=QR
void qr_gram_schmidt(const Matrix & A, Matrix & Q, Matrix & R);

// qr factorization using modified G-S algorithm
void qr_gram_schmidt_mod(const Matrix & A, Matrix & Q, Matrix & R);

// qr factorization using Householder reflections (stable)
void qr_householder(const Matrix & A, Matrix & Q, Matrix & R);

// qr factorization using Givens rotations
void qr_givens(const Matrix & A, Matrix & Q, Matrix & R);

// svd A = USV*
void svd(const Matrix & A, Matrix & U, Matrix & S, Matrix & V);

// LU decomposition
void lu(const Matrix & A, Matrix & L, Matrix & U);
#endif