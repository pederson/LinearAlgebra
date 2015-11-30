#ifndef _LINALGSOLVERS_H
#define _LINALGSOLVERS_H

#include <stdio.h>
#include <stdlib.h>
#include <random>

#include "Matrix.hpp"

// miscellaneous
template <typename T> int sgn(T val) {
    return (T(0) < val) - (val < T(0));
}


// collection of matrix/vector generators

// identity matrix
Matrix eye(unsigned int rows, unsigned int cols)
{
	Matrix out(rows, cols);
	out.fill(0);
	for (auto i=0; i<cols; i++)
	{
		out(i,i) = 1;
	}
	return out;
}

// hilbert matrix
// entries are H(i,j) = 1/(i+j+1) where i, j start at 0
Matrix hilb(unsigned int rows)
{
	Matrix out(rows, rows);
	for (auto i=0; i<rows; i++)
	{
		for (auto j=0; j<rows; j++)
		{
			out(i,j) = 1.0/(i+j+1);
		}
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

// random vector uniformly distributed [0,1]
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
	if (b.length() != L.rows())
	{
		throw "Matrix-Vector dimensions do not match!";
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
	if (b.length() != D.rows())
	{
		throw "Matrix-Vector dimensions do not match!";
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
void qr_gram_schmidt(const Matrix & A, Matrix & Qout)
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

	swap(Qout, Q);

	return;
}


// qr factorization using Gram-Schmidt algorithm A=QR
void qr_gram_schmidt(const Matrix & A, Matrix & Qout, Matrix & Rout)
{

	qr_gram_schmidt(A, Qout);

	// calculate R
	Matrix R = (~Qout)*A;
	swap(Rout, R);

	return;
}

// qr factorization using modified G-S algorithm
void qr_gram_schmidt_mod(const Matrix & A, Matrix & Q, Matrix & R);

// householder utility function
void householder_reflect_to_e1(const Vector & w, Vector & uout, double & b)
{
	Vector u = w;
	double sigma = sgn(w(0))*w.norm();
	u(0) += sigma;
	b = 2.0/(sigma*u(0));
	u /= u.norm();
	
	swap(u,uout);
	return;
}

// qr factorization using Householder reflections (stable)
void qr_householder(const Matrix & A, Matrix & Uout, Matrix & Rout)
{

// non-recursive way
	// std::size_t m,n;
	// A.size(m,n);

	// Matrix U(A);
	// U.fill(0);
	// Q = U;
	// R(A);
	// Vector w, u;

	// double sigma;
	// for (auto k=0; k<n; k++)
	// {
	// 	w = R.subcol(k, k, m-1);
	// 	sigma = sgn(w(1))*w.norm();

	// 	// calcuate the reflector
	// 	u = w;
	// 	u(1) += sigma;
	// 	u /= u.norm();

	// 	// apply the transformation
	// 	R(k, m-1, k, n-1) = R(k, m-1, k, n-1) - u*(~u*R(k, m-1, k, n-1))*2.0;
	// 	U.subcol(k, k, m-1) = u;
	// }

	// // construct Q



	// // reflect the first column to e_1
	// Vector u = A.col(0);
	// double sigma = 
	// u(0) += u.norm();
	// u /= u.norm();

// recursive way
	
	// reflect to -e_1
	Vector u;
	double b;
	householder_reflect_to_e1(A.col(0), u, b);

	// multiply through by H = I - 2u*u'
	Matrix B = A - 2.0*u*(~u*A);

	// recompose the full matrices
	// R
	Matrix R=A;
	R.col(0) = B.col(0);

	// U
	Matrix U(A.rows(), A.cols());
	U.fill(0);
	U.col(0) = u;

	if (A.cols() > 1 && A.rows() > 1)
	{
		// recurse on submatrix
		Matrix Rsub;
		Matrix Usub;
		Matrix Bsub = B(1, B.rows()-1, 1, B.cols()-1);

		qr_householder(Bsub, Usub, Rsub);

		R(1, R.rows()-1, 1, R.cols()-1) = Rsub;
		R.subrow(0, 1, B.cols()-1) = B.subrow(0, 1, B.cols()-1);
	
		U(1,U.rows()-1, 1,U.cols()-1) = Usub;
	}

	// swap
	swap(R, Rout);
	swap(U, Uout);

	return;
}

// qr factorization using Householder reflections (stable)
void qr_householder(const Matrix & A, Matrix & Uout, Matrix & Rout, Matrix & Qout)
{
	Matrix Uh, Rh;
	qr_householder(A, Uh, Rh);
	Vector u=Uh.col(0);
	Matrix I = eye(A.rows(), A.rows());
	Matrix Qh = I - 2.0*u*~u;

	// extract Q from Uh
	for (auto j=1; j<Uh.cols(); j++)
	{
		u = Uh.col(j);
		Matrix Qj = I - 2.0*u*~u;
		Qh *= Qj;
	}

	// extract the square R
	Matrix R = Rh(0, A.cols()-1, 0, A.cols()-1);

	// extract the appropriate size Q
	Matrix Q = Qh(0, Qh.rows()-1, 0, A.cols()-1);

	swap(Rout, R);
	swap(Qout, Q);
	swap(Uout, Uh);
	return;
}

// qr factorization using Householder with pivoting
void qr_householder(const Matrix & A, Matrix & U, Matrix & R, Matrix & Q, Matrix & P);

// qr factorization using Givens rotations
void qr_givens(const Matrix & A, Matrix & Q, Matrix & R);

// svd A = USV*
void svd(const Matrix & A, Matrix & U, Matrix & S, Matrix & V);

// LU decomposition
void lu(const Matrix & A, Matrix & Lout, Matrix & Uout)
{

	if (A.cols() != A.rows())
	{
		std::cout << "LU decomposition will fail for non-square matrices!" << std::endl;
		throw -1;
	}
	Matrix L(A);
	Matrix U(L(0, A.cols()-1, 0, A.cols()-1));

	//L(0,0) = 1.0;
	L.subcol(0, 0, L.rows()-1) /= A(0,0);
	

	if (A.cols() > 1 && A.rows()>1)
	{
		// zero out some stuff
		L.subrow(0, 1, L.cols()-1) *= 0.0;
		U.subcol(0, 1, L.cols()-1) *= 0.0;

		// form submatrix A
		Vector v;
		Vector l_21 = L.subcol(0, 1, L.rows()-1);
		Vector a_12t = A.subrow(0, 1, A.cols()-1);
		Matrix Asub(A.rows()-1, A.cols()-1);

		for (auto j=1; j<A.cols(); j++)
		{
			v = A.subcol(j, 1, A.cols()-1);
			Asub.col(j-1) = v;
		}
		Asub -= l_21*a_12t;

		// recurse on submatrix
		Matrix Lsub, Usub;
		lu(Asub, Lsub, Usub);

		// std::cout << "Asub: " << Asub << std::endl;
		// std::cout << "Lsub: " << Lsub << std::endl;
		// std::cout << "Usub: " << Usub << std::endl;
		// std::cout << "Lsub*Usub: " << Lsub*Usub << std::endl;

		// form the sub L and U
		L(1, L.rows()-1, 1, L.cols()-1) = Lsub;
		U(1, U.rows()-1, 1, U.cols()-1) = Usub;

	}

	swap(L, Lout);
	swap(U, Uout);

	return;
}


// pivoted lu decomposition
void lu(const Matrix & A, Matrix & Lout, Matrix & Uout, Matrix & P);


// randomized method for basis
void rand_basis(const Matrix & A, Matrix & Qout, unsigned int rank)
{

	unsigned int q=3;

	// get random matrix
	Matrix G = randmatn(A.cols(), rank+10);
	Matrix Y = A*G;
	Matrix U, Q, R;

	qr_householder(Y, U, R, Q);

	// power iteration correction
	for (auto j=1; j<q; j++)
	{
		Matrix Y = A*~A*Q;
		qr_householder(Y, U, R, Q);
	}

	// truncate basis for output
	Matrix Qr = Q(0, Q.rows()-1, 0, rank-1);
	swap(Qout, Qr);


	return;
}

// reduction to hessenberg form
// uses orthogonal similarity transformation
// this tridiagonalizes for hermitian matrices
void hessenberg(const Matrix & A, Matrix & Tout)
{

	std::size_t m,n;
	A.size(m,n);
	
	if (m != n)
	{
		std::cout << "Matrix must be square to reduce to Hessenberg form!" << std::endl;
		throw -1;
	}

	Matrix T(A);

	// loop and apply householder transformations successively
	for (auto k=0; k<m-2; k++)
	{
		Vector x = T.subcol(k, k+1, m-1);
		Vector v(x.length()); v.fill(0);
		v(0) = sgn(x(0))*x.norm();
		v += x;
		v /= v.norm(); // normalize
		T(k+1, m-1, k, n-1) -= 2.0*v*(~v*T(k+1,m-1,k,n-1));
		T(0, m-1, k+1, m-1) -= 2.0*(T(0,m-1, k+1, m-1)*v)*~v;
	}
	swap(Tout, T);
	return;
}

// qr algorithm (unshifted)
// reveals the smallest eigenvalue
// requires a tridiagonal matrix input
void qr_alg_tridiag_unshifted(const Matrix & A, Matrix & Tout)
{
	std::size_t m,n;
	A.size(m,n);

	Matrix T(A);
	double crit = abs(T(m-1, m-2));

	Matrix Q, R, U;
	unsigned int ctr=0;
	while(crit >= 1.0e-12)
	{
		qr_householder(T,U,R,Q);
		Matrix Tnew = R*Q;

		swap(Tnew, T);
		crit = abs(T(m-1, m-2));
	}

	swap(T,Tout);
	return;
}

// qr algorithm (Wilkinson Shift)
// reveals the smallest eigenvalue
// requires a tridiagonal matrix input
void qr_alg_tridiag_shifted(const Matrix & A, Matrix & Tout)
{
	std::size_t m,n;
	A.size(m,n);

	Matrix T(A);
	double crit = abs(T(m-1, m-2));
	double delta, mu;

	Matrix Q, R, U;
	unsigned int ctr=0;
	while(crit >= 1.0e-12)
	{

		delta = (T(m-2, m-2) - T(m-1,m-1))/2.0;
		mu = T(m-1, m-1) 
		   - sgn(delta)*T(m-1, m-2)*T(m-1,m-2)
		   /(abs(delta) + sqrt(delta*delta + T(m-1, m-2)*T(m-1, m-2)));

		qr_householder(T-mu*eye(n,n),U,R,Q);
		Matrix Tnew = R*Q + mu*eye(n,n);

		swap(Tnew, T);
		crit = abs(T(m-1, m-2));
	}

	swap(T,Tout);
	return;
}

// real, symmetric matrix eigenvalue decomp
void eig_symm(const Matrix & A, Matrix & Lout)
{

	std::size_t m,n;
	A.size(m,n);

	Matrix L = eye(n,n);

	// tridiagonalize
	Matrix T;
	hessenberg(A, T);

	for (auto i=0; i<m-1; i++)
	{
		// get eigenvalue
		Matrix Tnew;
		qr_alg_tridiag_shifted(T, Tnew);

		// capture eigenvalue
		L(m-i-1, m-i-1) = Tnew(m-i-1, m-i-1);

		// deflate
		T = Tnew(0,m-i-2, 0, m-i-2);
	}
	L(0,0) = T(0,0);

	swap(L,Lout);
	return;
}


// generalized complex eigenvalue decomposition

// schur decomposition

// compute the determinant

// estimate the condition number

// cholesky factorization

// direct inverse (slow in general)




#endif