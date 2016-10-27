#ifndef _LINALGSOLVERS_H
#define _LINALGSOLVERS_H

#include <stdio.h>
#include <stdlib.h>
#include <random>
#include <complex>

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

// return the diagonals of a matrix as a vector
Vector diag(const Matrix & mtx)
{
	unsigned int dim = std::min(mtx.rows(), mtx.cols());
	Vector out(dim);

	for (auto i=0; i<dim; i++)
	{
		out(i) = mtx(i,i);
	}
	return out;
}

// convert vector into square diagonal matrix
Matrix diag(const Vector & vct)
{
	Matrix out = eye(vct.length(), vct.length());

	for (auto i=0; i<vct.length(); i++)
	{
		out(i,i) = vct(i);
	}
	return out;
}


double trace(const Matrix & A){
	double s=0;
	for (auto j=0; j<A.rows(); j++) s+= A(j,j);
	return s;
}


double det2x2(const Matrix & A){
	if (A.rows() != 2 || A.cols() != 2){
		throw "Matrix is not 2x2!";
	}

	return A(0,0)*A(1,1) - A(0,1)*A(1,0);
}


void eig2x2(const Matrix & A, std::complex<double> & l1, std::complex<double> & l2){
	if (A.rows() != 2 || A.cols() != 2){
		throw "Matrix is not 2x2!";
	}

	l1 = 0.5*(A(0,0)+A(1,1)) + 0.5*sqrt(std::complex<double>((A(0,0)+A(1,1))*(A(0,0)+A(1,1)) - 4*(A(0,0)*A(1,1)-A(0,1)*A(1,0))));
	l2 = 0.5*(A(0,0)+A(1,1)) - 0.5*sqrt(std::complex<double>((A(0,0)+A(1,1))*(A(0,0)+A(1,1)) - 4*(A(0,0)*A(1,1)-A(0,1)*A(1,0))));
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




// collection of helper routines and matrix decompositions

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

// recursive golub-kahan bidiagonalization
// returns matrices C and D that contain the householder reflectors
// used to generate bidiagonal B
void golub_kahan(const Matrix & A, Matrix & Bout, Matrix & Cout, Matrix & Dout){

	Matrix B=A;
	Matrix G=A;
	Matrix C(A.rows(), A.cols());
	Matrix D(A.rows(), A.cols());

	// left multiply
	if (A.rows() > 1){
		// reflect to -e_1
		Vector u;
		double b;
		householder_reflect_to_e1(A.col(0), u, b);

		// multiply through by H = I - 2u*u'
		G -= 2.0*u*(~u*G);

		// C
		C.fill(0);
		C.col(0) = u;
	}


	// right multiply
	if (A.cols() > 1){
		Vector v;
		double b2;
		householder_reflect_to_e1(A.subrow(0, 1, A.cols()-1), v, b2);

		// multiply through by H = I - 2v*v'
		G -= 2.0*(G*v)*~v;

		// D
		D.fill(0);
		D.col(0) = v;
	}

	B.col(0) = G.col(0);
	B.row(0) = G.row(0);

	if (A.cols() > 1 && A.rows() > 1)
	{
		// recurse on submatrix
		Matrix Bsub;
		Matrix Csub;
		Matrix Dsub;
		Matrix Gsub = G(1, G.rows()-1, 1, G.cols()-1);

		golub_kahan(Gsub, Bsub, Csub, Dsub);

		B(1, B.rows()-1, 1, B.cols()-1) = Bsub;
		B.subrow(0, 1, G.cols()-1) = G.subrow(0, 1, G.cols()-1);
		B.subcol(0, 1, G.rows()-1) = G.subcol(0, 1, G.rows()-1);
	
		C(1,C.rows()-1, 1,C.cols()-1) = Csub;
		D(1,D.rows()-1, 1,D.cols()-1) = Dsub;

	}

	// swap
	swap(B, Bout);
	swap(C, Cout);
	swap(D, Dout);

	return;
}

// bidiagonalization Golub-Kahan
void bidiagonalize_golub_kahan(const Matrix & A, Matrix & B, Matrix & U, Matrix & V){

	// call golub_kahan

	// reform the U and V matrices

	// swap
}

// svd A = USV*
void svd(const Matrix & A, Matrix & U, Matrix & S, Matrix & V){

	// bidiagonalize
	Matrix U1, V1, B;
	bidiagonalize_golub_kahan(A, B, U1, V1);

	// eigenvalue decomposition on bidiagonal matrix
	
}

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


// cholesky factorization
// only works for symmetric, positive definite matrices
void cholesky(const Matrix & A, Matrix & Rout)
{

	std::size_t m,n;
	A.size(m,n);

	Matrix R = A;
	for (auto k=0; k<m; k++)
	{
		for (auto j=k+1; j<m; j++)
		{
			R.subrow(j, j, m-1) -= R.subrow(k, j, m-1)*R(k,j)/R(k,k); 
		}
		R.subrow(k,k,m-1) /= sqrt(R(k,k));
	}

	swap(R, Rout);
}


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

// qr algorithm (Double Shift) = Francis Algorithm
// reveals the smallest eigenvalue, possibly complex
// but stays in real arithmetic
// requires a hessenberg matrix input
void qr_alg_double_shifted(const Matrix & A, Matrix & Tout)
{
	std::size_t m,n;
	A.size(m,n);

	Matrix T(A);
	// double crit = abs(T(m-1, m-2));
	double delta, mu;

	double s, t, x, y, z;
	Vector u, v(3), w(2), uw;
	double b;
	unsigned int r;

	Matrix Q, R, U;
	unsigned int ctr=0;
	double eps = 1.0e-16;
	std::size_t p = n, q;
	while(p > 2)
	{

		// std::cout << "H(" << ctr << "): " << T << '\n' << std::endl;
		
		q = p-1;

		// apply double shift
		s = trace(T(q-1, p-1, q-1, p-1)); 	// trace of submatrix
		t = det2x2(T(q-1, p-1, q-1, p-1));	// det of submatrix
		x = T(0,0)*T(0,0) + T(0,1)*T(1,0) - s*T(0,0) + t;
		y = T(1,0)*(T(0,0)+T(1,1)-s);
		z = T(1,0)*T(2,1);

		// chase the bulge ... heh
		for (std::size_t k=0; k<=p-3; k++){

			// std::cout << "p: " << p << std::endl;


			v(0) = x;
			v(1) = y;
			v(2) = z;
			householder_reflect_to_e1(v, u, b);

			// apply to left
			r = std::max(std::size_t(1),k);
			T(k, k+2, r-1,n-1) -= 2.0*u*(~u*T(k, k+2, r-1, n-1));
			// apply to right
			r = std::min(k+4, p);
			T(0, r-1, k, k+2) -= 2.0*(T(0, r-1, k, k+2)*u)*~u;

			// std::cout << "b: " << b << std::endl;

			// calc new x, y, z
			x = T(k+1, k);
			y = T(k+2, k);

			if (k<p-3) z=T(k+3, k);
		}

		// // givens rotation P for final submatrix (or householder)
		// w(0) = x;
		// w(1) = y;
		// householder_reflect_to_e1(w, uw, b);
		// // apply to left
		// T(q-1, p-1, p-3,n-1) -= 2.0*uw*((~uw)*T(q-1, p-1, p-3, n-1));
		// // apply to right
		// T(0,p-1, p-2, p-1) -= 2.0*(T(0, p-1, p-2, p-1)*uw)*~uw;


		// determine if solution meets criteria
		if (abs(T(p-1, q-2)) < eps*(abs(T(q-1, q-1)) + abs(T(p-1, p-1)))){
			
			if (abs(T(p-1, q-1) < eps)){
				T(p-1, q-1) = 0;
				p--;
			}
			else{
				T(p-2, q-2) = 0;
				p -=2;
			}
			// q = p-1;
		}
		else if (abs(T(p-2, q-2)) < eps*(abs(T(q-2, q-2)) + abs(T(q-1, q-1)))){
			T(p-2, q-2) = 0;
			p -= 2;
			// q = p-1;
		}
		ctr++;
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


// real, nonhermitian matrix eigenvalue decomp
void eig(const Matrix & A, Matrix & Lout, std::vector<std::complex<double>> & eigout)
{

	std::size_t m,n;
	A.size(m,n);

	Matrix L = eye(n,n);

	// hessenberg form
	Matrix T, H;
	hessenberg(A, H);

	qr_alg_double_shifted(H, T);

	std::vector<std::complex<double>> eigs(n);
	std::complex<double> l1, l2;
	unsigned int i=0;
	unsigned int block;
	while (i<T.rows()){
		if (i==T.rows()-1){
			eigs[i] = T(i,i);
			break;
		}


		if (abs(T(i+1, i)) < 1.0e-16){
			block = 1;
			eigs[i] = T(i,i);
			i++;
		} 
		else{
			block = 2;
			eig2x2(T(i, i+1, i, i+1), l1, l2);
			eigs[i] = l1;
			eigs[i+1] = l2;
			i+=2;
		}
	}
	
	swap(T, Lout);
	eigout = eigs;
	return;
}

// *** iterative methods *** //

// steepest descent
// for square, SPD matrices only
// uses the x argument as x0
void steepest_descent(const Matrix & mtx, const Vector & b, Vector & x)
{

	Vector rv = b-mtx*x;
	double resid = rv.norm();
	double threshold = std::max(resid*1.0e-10, 1.0e-6);
	double alpha;
	Vector matvec;

	unsigned int ctr=0;
	while (resid > threshold)
	{
		// one matrix-vector product
		matvec = mtx*rv;

		// direction is residual
		alpha = Vector::dot(rv,rv)/Vector::dot(rv,matvec);

		// update x
		x += alpha*rv;

		if (ctr%50 == 0)
		{
			// recalculate the exact resid every 50 iters
			rv = b - mtx*x;
		}
		else
		{
			// update residual
			rv -= alpha*matvec;
		}

		resid = rv.norm();
		//std::cout << "residual: " << resid << std::endl;

		ctr++;
	}

	std::cout << "iterated: " << ctr << " times" << std::endl;


	return;
}


// conjugate gradient
// for square, SPD matrices only
// uses the x argument as x0
void conjugate_gradient(const Matrix & mtx, const Vector & b, Vector & x)
{
	Vector rv = b-mtx*x;

	double resid = rv.norm();
	double threshold = std::max(resid*1.0e-16, 1.0e-15);

	// first iteration here
	Vector d = rv;
	Vector matvec = mtx*d;
	double alpha = Vector::dot(d,rv)/Vector::dot(d,matvec);
	x = alpha*b;

	// continue iterating
	unsigned int ctr=0;
	while (resid > threshold)
	{
		// calculate residual
		rv = b - mtx*x;
		resid = rv.norm();

		// pick new direction
		d = rv - Vector::dot(rv, matvec)/Vector::dot(d, matvec)*d;

		// one matrix-vector product
		matvec = mtx*d;

		// calculate new step size
		alpha = Vector::dot(d,rv)/Vector::dot(d,matvec);

		// update x
		x += alpha*d;
		
		ctr++;
	}

	std::cout << "iterated: " << ctr << " times" << std::endl;

	return;
}


// generalized complex eigenvalue decomposition

// schur decomposition

// compute the determinant

// estimate the condition number

// direct inverse (slow in general)




#endif