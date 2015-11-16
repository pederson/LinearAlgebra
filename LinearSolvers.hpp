#ifndef _LINALGSOLVERS_H
#define _LINALGSOLVERS_H

#include <stdio.h>
#include <stdlib.h>
#include <random>

#include "Matrix.hpp"

// collection of solvers for linear algebra problems

// upper triangular matrix Ux = b
Vector upper_triangular_solve(const Matrix & U, const Vector & b)
{
	if (U.rows() != U.cols())
	{
		throw "Matrix is not square!";
	}
	if (b.length() != U.rows())
	{
		throw "Matrix-Vector dimensions do not match!";
	}

	Vector out(b);
	double sum;
	// last element
	out(U.rows()-1) = b(U.rows()-1)/U(U.rows()-1, U.rows()-1);
	// backward substitution
	for (int i=U.rows()-2; i>=0; i--)
	{
		sum = Vector_Proxy::dot(U.subrow(i, i+1, U.rows()-1), out(i+1, U.rows()-1));
		out(i) = (b(i) - sum)/U(i,i);
	}
	return out;
}

// lower triangular matrix Lx = b
Vector lower_triangular_solve(const Matrix & L, const Vector & b)
{
	if (L.rows() != L.cols())
	{
		throw "Matrix is not square!";
	}

	Vector out(b);
	double sum;
	// first element
	out(0) = b(0)/L(0,0);
	// forward substitution
	for (auto i=1; i<L.rows(); i++)
	{
		sum = Vector_Proxy::dot(L.subrow(i, 0, i-1), out(0, i-1));
		out(i) = (b(i) - sum)/L(i,i);
	}
	return out;
}


// diagonal matrix Dx = b
Vector diagonal_solve(const Matrix & D, const Vector & b)
{

	if (D.rows() != D.cols())
	{
		throw "Matrix is not square!";
	}

	Vector out(b);
	for (auto i=0; i<D.rows(); i++)
	{
		out(i) /= D(i,i);
	}
	return out;
}

// unitary matrix Ux = b
Vector unitary_solve(const Matrix & U, const Vector & b)
{
	if (U.rows() != b.length())
	{
		throw "Matrix-Vector dimensions do not match!";
	}

	Vector out(U.cols());
	for (auto i=0; i<U.cols(); i++)
	{
		out(i) = Vector_Proxy::dot(U.col(i), b.col(0));
	}

	return out;
}

// tridiagonal matrix Tx = b (Thomas algorithm)
Vector tridiagonal_solve(const Matrix & T, const Vector & b);




// collection of helper routines

// qr factorization using Gram-Schmidt algorithm A=QR
void qr_gram_schmidt(const Matrix & A, Matrix & Qout, Matrix & Rout)
{
	Matrix Q(A.rows(), A.cols());

	Vector q(A.rows());
	Vector z(A.rows());

	// first one
	q = A.col(0);
	q /= q.norm();
	Q.col(0) = q;

	// loop over the rest
	for (auto j=1; j<A.cols(); j++)
	{

		z.fill(0);
		// project out the already found columns
		for (auto k=0; k<j; k++)
		{
			z += Q.col(k)*Vector_Proxy::dot(A.col(j), Q.col(k));
		}

		// subtract
		q = A.col(j);
		q -= z;

		// normalize
		q /= q.norm();

		// insert into output Q
		Q.col(j) = q;
	}

	// calculate R
	Matrix R = (~Q)*A;

	swap(Qout, Q);
	swap(Rout, R);

	return;
}

// qr factorization using modified G-S algorithm
void qr_gram_schmidt_mod(const Matrix & A, Matrix & Q, Matrix & R);

// qr factorization using Householder reflections (stable)
void qr_householder(const Matrix & A, Matrix & Q, Matrix & R);

// qr factorization using Householder with pivoting
void qr_householder(const Matrix & A, Matrix & Q, Matrix & R, Matrix & P);

// qr factorization using Givens rotations
void qr_givens(const Matrix & A, Matrix & Q, Matrix & R);

// svd A = USV*
void svd(const Matrix & A, Matrix & U, Matrix & S, Matrix & V);

// LU decomposition
void lu(const Matrix & A, Matrix & L, Matrix & U);





// collection of matrix/vector generators

// identity matrix
Matrix eye(unsigned int size)
{
	Matrix out(size, size);
	out.fill(0);
	for (auto i=0; i<size; i++)
	{
		out(i,i) = 1;
	}
	return out;
}

// random matrix uniformly distributed [0,1]
Matrix randmat(unsigned int rows, unsigned int cols)
{
	// seed
	std::default_random_engine generator;
	std::uniform_real_distribution<double> distrib(0.0,1.0);

	Matrix out(rows, cols);
	for (auto i=0; i<rows; i++){
		for (auto j=0; j<cols; j++){
			out(i,j) = distrib(generator);
		}
	}

	return out;
}

// random matrix normally distributed
Matrix randmatn(unsigned int rows, unsigned int cols)
{
	std::default_random_engine generator;
	std::normal_distribution<double> distrib(0.0,1.0);

	Matrix out(rows, cols);
	for (auto i=0; i<rows; i++){
		for (auto j=0; j<cols; j++){
			out(i,j) = distrib(generator);
		}
	}

	return out;
}

// random vector uniformly distributed [0.1]
Vector randvec(unsigned int length)
{
	Matrix m = randmat(length,1);
	Vector out = m.col(0);
	return out;
}

// random vector normally distributed
Vector randvecn(unsigned int length)
{
	Matrix m=randmatn(length,1);
	Vector out = m.col(0);
	return out;
}


#endif