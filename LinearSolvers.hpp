#ifndef _LINALGSOLVERS_H
#define _LINALGSOLVERS_H

// collection of solvers for linear algebra problems

// upper triangular matrix Ux = b
Vector upper_triangular_solve(const Matrix & U, const Vector & b);
{
	if (L.rows() != U.cols())
	{
		throw "Matrix is not square!";
	}

	Vector out(b);
	double sum;
	// last element
	out(U.rows()-1) = b(U.rows()-1)/U(U.rows()-1, U.rows()-1);
	// backward substitution
	for (auto i=U.rows()-2; i>=0; i--)
	{
		sum = Vector_Proxy::dot(U(i, i, i+1, U.rows()-1), out(i+1, U.rows()-1));
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
		sum = Vector_Proxy::dot(L(i, i, 0, i-1), out(0, i-1));
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
	if (U.rows() != U.cols())
	{
		throw "Matrix is not square (and therefore not unitary)!";
	}

	Vector out(U.rows());
	for (auto i=0; i<U.rows(); i++)
	{
		out(i) = Vector_Proxy::dot(U.col(i)*b.col(0));
	}

	return out;
}

// tridiagonal matrix Tx = b (Thomas algorithm)
Vector tridiagonal_solve(const Matrix & T, const Vector & b);




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





// collection of matrix/vector generators
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


#endif