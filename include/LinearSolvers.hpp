#ifndef _LINALGSOLVERS_H
#define _LINALGSOLVERS_H

#include <stdio.h>
#include <stdlib.h>
#include <complex>
#include <algorithm>
#include <set>
#include <iterator>

#include "Matrix.hpp"
#include "Vector.hpp"
#include "SparseVector.hpp"
#include "SparseMatrix.hpp"
#include "Preconditioner.hpp"

#include "VectorTools.hpp"
#include "Traits.hpp"

// using namespace libra;

namespace libra{

template <typename MatrixType>
decltype(auto) diag(const MatrixType & mtx)
{
	typedef decltype(mtx(0,0)) 	ScalarType;
	std::size_t dim = std::min(mtx.rows(), mtx.cols());
	Vector<ScalarType, libra::dynamic_size> out(dim);

	for (auto i=0; i<dim; i++)
	{
		out(i) = mtx(i,i);
	}
	return out;
}




// sparse BiConjugate Gradient Stabilized method
template <typename MatrixType, typename VectorTypeIn, typename VectorTypeOut>
unsigned int bicgstab(const MatrixType & A, const VectorTypeIn & b, VectorTypeOut & x, unsigned int max_iters, double res_thresh=1.0e-15){
	typedef typename std::remove_const<typename std::remove_reference<typename type_traits::contained_type<VectorTypeIn>::type>::type>::type 	ScalarType;

	static_assert(!std::is_const<ScalarType>::value, "Type Must not be const!");

	// initialize stuff
	libra::Vector<ScalarType, libra::dynamic_size> tmp = x;
	A.vmult(x, tmp);
	libra::Vector<ScalarType, libra::dynamic_size> r = b - tmp;
	double resid = 1.0;
	unsigned int it = 0;
	libra::Vector<ScalarType, libra::dynamic_size> rhat = r;
	ScalarType rho=1, alpha=1, omega=1;
	libra::Vector<ScalarType, libra::dynamic_size> v = r; vector::fill(v, 0);
	libra::Vector<ScalarType, libra::dynamic_size> p = r; vector::fill(p, 0);
	libra::Vector<ScalarType, libra::dynamic_size> s = r;
	libra::Vector<ScalarType, libra::dynamic_size> t = r;
	ScalarType beta, rho_old;

	// std::cout << "starting iterations" << std::endl;

	while (resid > res_thresh && it < max_iters){

		// enforce biorthogonality and biconjugacy
		rho_old = rho;
		rho = vector::inner_product(rhat, r);
		beta = (rho/rho_old)*(alpha/omega);
		p = r + beta*(p-omega*v);

		A.vmult(p, v);
		// v = A*p;

		alpha = rho/vector::inner_product(rhat,v);
		s = r - alpha*v;

		A.vmult(s, t);
		// t = A*s;
		
		omega = vector::inner_product(t,s)/vector::inner_product(t,t);

		x += alpha*p + omega*s;

		// recalculate residual
		r = s - omega*t;

		// update counters
		it++;
		resid = std::abs(vector::norm_2(r));

		std::cout << "it: " << it << " resid: " << resid << std::endl;
	}

	std::cout << "iterated: " << it << " times" << std::endl;
	return it;
}




// BiConjugate Gradient Stabilized L method
template <typename MatrixType, typename VectorTypeIn, typename VectorTypeOut>
unsigned int bicgstab_l(unsigned int l, const MatrixType & A, const VectorTypeIn & b, VectorTypeOut & x, unsigned int max_iters, double res_thresh=1.0e-15){
	typedef typename std::remove_const<typename std::remove_reference<typename type_traits::contained_type<VectorTypeIn>::type>::type>::type 	ScalarType;

	static_assert(!std::is_const<ScalarType>::value, "Type Must not be const!");
	// std::cout << "pre matmult" << std::endl;


	// initialize stuff
	libra::Vector<ScalarType, libra::dynamic_size> tmp(vector::length(x));// = x;
	// std::cout << "pre matmult" << std::endl;

	A.vmult(x, tmp);
	// std::cout << " matmult" << std::endl;
	// tmp = A*x;
	libra::Vector<ScalarType, libra::dynamic_size> r = b - tmp; // residual vector
	libra::Vector<ScalarType, libra::dynamic_size> rtilde = r;
	double resid = 1.0;
	unsigned int it = 0;

	// std::cout << "resid" << std::endl;

	ScalarType beta, rho0 = 1, rho1, alpha = 0, omega = 1;
	
	unsigned int N = vector::length(x);
	libra::matrix::Matrix<ScalarType, libra::dynamic_size, libra::dynamic_size> Tau(l+1,l+1); // fill(Tau, 0);
	libra::matrix::Matrix<ScalarType, libra::dynamic_size, libra::dynamic_size> R(l+1, N);// fill(R, 0);
	libra::matrix::Matrix<ScalarType, libra::dynamic_size, libra::dynamic_size> U(l+1, N);// fill(U, 0);
	libra::Vector<ScalarType, libra::dynamic_size> gamma(l+1);
	libra::Vector<ScalarType, libra::dynamic_size> gammap(l+1);
	libra::Vector<ScalarType, libra::dynamic_size> gammapp(l+1);
	libra::Vector<ScalarType, libra::dynamic_size> sigma(l+1);

	// std::cout << "allocations" << std::endl;

	R.row(0) = r;
	vector::fill(U.row(0), 0);

	// std::cout << "set row 0 to resid" << std::endl;

	// std::cout << "starting iterations" << std::endl;

	while (resid > res_thresh && it < max_iters){
		
		rho0 = -omega*rho0;

		// ***  BiCG part *** //
		for (auto j=0; j<l; j++){
			// std::cout << " ********* j = " << j << std::endl;
			rho1 = vector::inner_product(R.row(j), rtilde);
			beta = alpha*rho1/rho0;
			rho0 = rho1;
			

			// std::cout << "did inner prod" << std::endl;

			for (auto i=0; i<=j; i++){
				U.row(i) = R.row(i) - beta*U.row(i);
			}

			// std::cout << "did R - beta*U" << std::endl;

			// U.row(j+1) = A*U.row(j);
			A.vmult(U.row(j), U.row(j+1));
			// libra::vector::write<true>(U.row(j+1));
			// libra::vector::write<true>(rtilde);
			// std::cout << "did vmult " << vector::inner_product(U.row(j+1), rtilde) << std::endl;
			alpha = rho0/vector::inner_product(U.row(j+1), rtilde);

			// std::cout << "did inner prod 2" << std::endl;

			for (auto i=0; i<=j; i++){
				R.row(i) = R.row(i) - alpha*U.row(i+1);
			}

			// R.row(j+1) = A*R.row(j);
			A.vmult(R.row(j), R.row(j+1));
			x += alpha*U.row(0);
		}



		// ***  Modified Gram-Schmidt part *** //
		for (auto j=1; j<=l; j++){
			for (auto i=1; i<j; i++){
				Tau(i,j) = 1.0/sigma(i)*vector::inner_product(R.row(j), R.row(i));
				R.row(j) = R.row(j) - Tau(i,j)*R.row(i);
			}

			sigma(j) = vector::inner_product(R.row(j), R.row(j));
			gammap(j) = 1.0/sigma(j)*vector::inner_product(R.row(0), R.row(j));
		}
		gamma(l) = gammap(l);
		omega = gamma(l);


		// ***  solve Tau*gamma = gammap *** //
		for (auto j=l-1; j>=1; j--){
			ScalarType sum_Taugamma = Tau(j, j+1)*gamma(j+1);
			for (auto i=j+2; i<=l; i++) sum_Taugamma += Tau(j,i)*gamma(i);
			gamma(j) = gammap(j) - sum_Taugamma;
		}

		// ***  solve gammapp = Tau*S*gamma
		for (auto j=1; j<l; j++){
			ScalarType sum_Taugamma = Tau(j,j+1)*gamma(j+2);
			for (auto i=j+2; i<l; i++) sum_Taugamma += Tau(j,i)*gamma(i+1);
			gammapp(j) = gamma(j+1) + sum_Taugamma;
		}


		// ***  update  *** //
		x = x + gamma(1)*R.row(0);
		R.row(0) = R.row(0) - gammap(l)*R.row(l);
		U.row(0) = U.row(0) - gamma(l)*U.row(l);
		for (auto j=1; j<l; j++){
			U.row(0) = U.row(0) - gamma(j)*U.row(j);
			x += gammapp(j)*R.row(j);
			R.row(0) = R.row(0) - gammap(j)*R.row(j);
		}

		// update counters
		it+=l;
		resid = std::abs(norm_2(R.row(0)));

		std::cout << "it: " << it << " resid: " << resid << std::endl;
	}

	std::cout << "iterated: " << it << " times" << std::endl;
	return it;
}




} // end namespace libra



























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

// return the diagonals of a sparse matrix as a vector
Vector diag(const SparseMatrix & mtx){
	unsigned int dim = std::min(mtx.rows(), mtx.cols());
	Vector out(dim);

	for (auto i=0; i<dim; i++)
	{
		out(i) = mtx.get(i,i);
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

// convert vector into square diagonal sparse matrix
// off is the offset of the diagonal (+ is up, - is down)
SparseMatrix spdiag(const Vector & vct, int off)
{
	SparseMatrix out(vct.length()+abs(off), vct.length()+abs(off));

	unsigned int ioff=0, joff=0;
	if (off < 0) ioff=abs(off);
	if (off > 0) joff=abs(off);
	for (auto i=0; i<vct.length(); i++)
	{
		out.set(i+ioff,i+joff,vct(i));
	}
	return out;
}

// return strictly upper triangular part of a matrix
Matrix strictly_upper(const Matrix & mtx){
	Matrix out(mtx.rows(), mtx.cols());
	out.fill(0);
	for (auto i=0; i<mtx.rows(); i++){
		for (auto j=i+1; j<mtx.cols(); j++){
			out(i,j) = mtx(i,j);
		}
	}
	return out;
}

// return strictly upper triangular part of a sparse matrix
SparseMatrix strictly_upper(const SparseMatrix & mtx){
	SparseMatrix out(mtx.rows(), mtx.cols());
	auto data = mtx.data();
	auto rowptr = mtx.row_ptr();

	for (auto i=0; i<mtx.rows(); i++){
		auto it = rowptr[i];
		while (it != rowptr[i+1] && it !=data.end()){
			if (it->first > i) out.set(i,it->first,it->second);
			it++;
		}
	}

	return out;
}

// return strictly lower triangular part of a matrix
Matrix strictly_lower(const Matrix & mtx){
	Matrix out(mtx.rows(), mtx.cols());
	out.fill(0);
	for (auto i=0; i<mtx.rows(); i++){
		for (auto j=0; j<i; j++){
			out(i,j) = mtx(i,j);
		}
	}
	return out;
}

// return strictly lower triangular part of a sparse matrix
SparseMatrix strictly_lower(const SparseMatrix & mtx){
	SparseMatrix out(mtx.rows(), mtx.cols());
	auto data = mtx.data();
	auto rowptr = mtx.row_ptr();

	for (auto i=0; i<mtx.rows(); i++){
		auto it = rowptr[i];
		while (it != rowptr[i+1] && it !=data.end()){
			if (it->first < i) out.set(i,it->first,it->second);
			it++;
		}
	}

	return out;
}



double trace(const Matrix & A){
	double s=0;
	for (auto j=0; j<A.rows(); j++) s+= A(j,j);
	return s;
}

double trace(const SparseMatrix & A){
	double s=0;
	auto data = A.data();
	auto rowptr = A.row_ptr();

	for (auto i=0; i<A.rows(); i++){
		auto it = rowptr[i];
		while (it != rowptr[i+1] && it !=data.end()){
			if (it->first == i) s+=it->second;
			it++;
		}
	}
	return s;
}


double det2x2(const Matrix & A){
	if (A.rows() != 2 || A.cols() != 2){
		std::cout << "Matrix is not 2x2" << std::endl;
		throw -1;
	}

	return A(0,0)*A(1,1) - A(0,1)*A(1,0);
}


void eig2x2(const Matrix & A, std::complex<double> & l1, std::complex<double> & l2){
	if (A.rows() != 2 || A.cols() != 2){
		std::cout << "Matrix is not 2x2" << std::endl;
		throw -1;
	}

	l1 = 0.5*(A(0,0)+A(1,1)) + 0.5*sqrt(std::complex<double>((A(0,0)+A(1,1))*(A(0,0)+A(1,1)) - 4*(A(0,0)*A(1,1)-A(0,1)*A(1,0))));
	l2 = 0.5*(A(0,0)+A(1,1)) - 0.5*sqrt(std::complex<double>((A(0,0)+A(1,1))*(A(0,0)+A(1,1)) - 4*(A(0,0)*A(1,1)-A(0,1)*A(1,0))));
}

void invert2x2(const Matrix & A, Matrix & A_inv){
	double d = det2x2(A);
	Matrix inv(2,2);
	inv(0,0) = A(1,1);
	inv(1,1) = A(0,0);
	inv(0,1) = -A(0,1);
	inv(1,0) = -A(1,0);
	inv = inv/d;
	swap(inv, A_inv);
}



// collection of solvers for linear algebra problems

// upper triangular matrix Ux = b
Vector upper_triangular_solve(const Matrix & U, const Vector & b)
{
	if (U.rows() != U.cols())
	{
		std::cout << "Matrix must be square to upper triangular solve!" << std::endl;
		throw -1;
	}
	if (b.length() != U.rows())
	{
		std::cout << "Matrix-Vector dimensions do not match in upper triangular solve!" << std::endl;
		throw -1;
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


// upper triangular sparse matrix Ux = b
Vector upper_triangular_solve(const SparseMatrix & U, const Vector & b)
{
	if (U.rows() != U.cols())
	{
		std::cout << "Matrix must be square to upper triangular solve!" << std::endl;
		throw -1;
	}
	if (b.length() != U.rows())
	{
		std::cout << "Matrix-Vector dimensions do not match in upper triangular solve!" << std::endl;
		throw -1;
	}

	Vector out(b);
	SparseVector row;
	double sum;
	// last element
	out(U.rows()-1) = b(U.rows()-1)/U.get(U.rows()-1, U.rows()-1);
	// backward substitution
	for (int i=U.rows()-2; i>=0; i--)
	{
		row = U.row(i);
		sum = 0;
		auto data = row.data();
		for (auto it=data.begin(); it != data.end(); it++){
			if (it->first >= i+1) sum += it->second*out(it->first);
		}
		out(i) = (b(i) - sum)/U.get(i,i);
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

// lower triangular sparse matrix Lx = b
Vector lower_triangular_solve(const SparseMatrix & L, const Vector & b)
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
	SparseVector row;
	double sum;
	// first element
	out(0) = b(0)/L.get(0,0);
	// forward substitution
	for (auto i=1; i<L.rows(); i++)
	{
		row = L.row(i);
		sum = 0;
		auto data = row.data();
		for (auto it=data.begin(); it != data.end(); it++){
			if (it->first < i) sum += it->second*out(it->first);
		}
		//sum = Vector_Proxy::dot(L.subrow(i, 0, i-1), out(0, i-1));
		out(i) = (b(i) - sum)/L.get(i,i);
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

// diagonal sparse matrix Dx = b
Vector diagonal_solve(const SparseMatrix & D, const Vector & b)
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
		out(i) /= D.get(i,i);
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

// // tridiagonal matrix Tx = b (Thomas algorithm)
// Vector tridiagonal_solve(const Matrix & T, const Vector & b){
// 	if (T.rows() != T.cols())
// 	{
// 		throw "Matrix is not square!";
// 	}
// 	if (b.length() != T.rows())
// 	{
// 		throw "Matrix-Vector dimensions do not match!";
// 	}
// }




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


// givens utility function
void givens_rotate_to_e1(double a, double b, double & c, double & s)
{
	double r;
	if (b==0){
		c=1; 
		s=0;
	}
	else{
		if (abs(b) > abs(a)){
			r = -a/b;
			s = 1.0/sqrt(1+r*r);
			c = s*r;
		}
		else{
			r = -b/a;
			c = 1.0/sqrt(1+r*r);
			s = c*r;
		}
	}
	return;
}



// givens rotation matrix 
// premultiplied on a matrix A -> G*A at position (i,k)
// where i and k are indices STARTING AT 0
void givens_premultiply(Matrix & A, double c, double s, std::size_t i, std::size_t k)
{
	double t1, t2;
	for (auto j=0; j<A.cols(); j++){
		t1 = A(i, j);
		t2 = A(k, j);
		A(i,j) = c*t1 - s*t2;
		A(k,j) = s*t1 + c*t2;
	}
}


// givens rotation matrix 
// postmultiplied on a matrix A -> A*G at position (i,k)
// where i and k are indices STARTING AT 0
void givens_postmultiply(Matrix & A, double c, double s, std::size_t i, std::size_t k)
{
	double t1, t2;
	for (auto j=0; j<A.rows(); j++){
		t1 = A(j,i);
		t2 = A(j,k);
		A(j,i) = c*t1 - s*t2;
		A(j,k) = s*t1 + c*t2;
	}
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


// incomplete LU factorization
// only works for square
void ilu(const SparseMatrix & A, SparseMatrix & Lout, SparseMatrix & Uout)
{

	std::size_t m,n;
	A.size(m,n);

	SparseMatrix L = A;
	auto data = L.data();
	auto rowptr = L.row_ptr();
	Vector dg = diag(L);
	double e;

	
	// loop over rows
	for (auto i=0; i<m; i++){
		auto it = rowptr[i];
		// loop over k columns
		while (it != rowptr[i+1] && it !=data.end()){
			unsigned int k = it->first;
			if(k > i) break;

			double Lii = dg(k);
			double Lik = it->second;
			L.set(i,k,Lik/Lii);

			// loop over columns j
			auto iti = rowptr[i];
			auto itk = rowptr[k];
			while (iti != rowptr[i+1] && iti != data.end()){
				unsigned int j = iti->first;
				if(j <= k){
					iti++;
					continue;
				}
				// find item Aik and Akj and subtract from item Aij
				L.add(i,j, -Lik/Lii*L.get(k,j));
				iti++;
			}

			it++;
		}


	}
	//*/

	SparseMatrix U = strictly_upper(A) + spdiag(diag(A));
	L = strictly_lower(L) + spdiag(diag(L));
	swap(L, Lout);
	swap(U, Uout);
}


// cholesky factorization
// only works for symmetric, positive definite matrices
// returns a lower triangular matrix L
void cholesky(const Matrix & A, Matrix & Lout)
{

	std::size_t m,n;
	A.size(m,n);

	Matrix L(A);
	for (auto k=1; k<=m; k++)
	{
		for (auto j=k+1; j<=m; j++)
		{
			L.subrow(j-1, j-1, m-1) -= L.subrow(k-1, j-1, m-1)*L(k-1,j-1)/L(k-1,k-1); 
		}
		L.subrow(k-1,k-1,m-1) = L.subrow(k-1,k-1,m-1)/sqrt(L(k-1,k-1));
	}

	for (auto k=1; k<n; k++){
		L.subcol(k-1, k, m-1).fill(0);
	}
	L = ~L;
	swap(L, Lout);
}


// incomplete cholesky factorization
// only works for symmetric, positive definite matrices
// returns LOWER TRIANGULAR MATRIX (L)
void icholesky(const SparseMatrix & A, SparseMatrix & Lout)
{
	std::size_t m,n;
	A.size(m,n);

	SparseMatrix L = A;
	auto data = L.data();
	auto rowptr = L.row_ptr();
	Vector dg = diag(L);
	double e;

	
	// loop over rows
	for (auto i=0; i<m; i++){
		auto it = rowptr[i];
		// loop over k columns
		while (it != rowptr[i+1] && it !=data.end()){
			unsigned int k = it->first;
			if(k > i) break;

			double Lii = dg(k);
			double Lik = it->second;
			L.set(i,k,Lik/Lii);

			// loop over columns j
			auto iti = rowptr[i];
			auto itk = rowptr[k];
			while (iti != rowptr[i+1] && iti != data.end()){
				unsigned int j = iti->first;
				if(j <= k){
					iti++;
					continue;
				}
				// find item Aik and Akj and subtract from item Aij
				L.add(i,j, -Lik/Lii*L.get(k,j));
				iti++;
			}

			it++;
		}


	}
	//*/

	swap(L, Lout);
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

		// double c,s;
		// givens_rotate_to_e1(x, y, c, s);
		// givens_premultiply(T, c, s, p-1, p-3);
		// givens_postmultiply(T, c, s, p-1, p-3);



		// determine if solution meets criteria
		if (abs(T(p-1, q-2)) < eps*(abs(T(q-1, q-1)) + abs(T(p-1, p-1)))){
			
			// T(p-1, q-1) = 0;
			// p--;

			if (abs(T(p-1, q-1) < eps)){
				T(p-1, q-1) = 0;
				p--;
			}
			else{
				T(p-2, q-2) = 0;
				p -=2;
			}

		}
		else if (abs(T(p-2, q-2)) < eps*(abs(T(q-2, q-2)) + abs(T(q-1, q-1)))){
			T(p-2, q-2) = 0;
			p -= 2;
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


void singular_values_sq(const Matrix & A, Matrix & Sout){
	Matrix C = (~A)*A;

	Matrix S;
	eig_symm(C, S);
	swap(S, Sout);
}


// use a QR decomposition to solve A*A_inv = I ----> Q*R*A_inv = I
void invert(const Matrix & A, Matrix & A_inv){
	// check dimensions
	if (A.rows() != A.cols()){
		std::cout << "Matrix must be square to invert!" << std::endl;
		throw -1;
	}

	// first, QR factorize A
	Matrix Q, R, U;
	qr_householder(A, U, R, Q);
	Matrix inv(A);
	Vector rhs;

	for (auto i=0; i<A.rows(); i++){
		rhs = Q.subrow(i, 0, A.cols()-1);
		inv.subcol(i, 0, A.rows()-1) = upper_triangular_solve(R, ~rhs);
	}

	swap(inv, A_inv);
}

Matrix inv(const Matrix & A){
	Matrix out(A);
	invert(A, out);
	return out;
}

// *** iterative methods *** //

// jacobi iteration
// for square diagonally dominant matrices
unsigned int jacobi(const Matrix & A, const Vector & b, Vector & x, unsigned int max_iters, double res_thresh=1.0e-15){

	Vector r = b-A*x;
	Vector d;
	double resid = norm_2(r);
	unsigned int it = 0;
	
	while (resid > res_thresh && it < max_iters){
		// d = b - R*xk;
		x += diagonal_solve(A, r);

		it++;
		r = b - A*x;
		resid = norm_2(r);
	}

	return it;
}

// sparse jacobi iteration
// for square diagonally dominant matrices
unsigned int jacobi(const SparseMatrix & A, const Vector & b, Vector & x, unsigned int max_iters, double res_thresh=1.0e-15){

	Vector r = b-A*x;
	Vector d;
	double resid = norm_2(r);
	unsigned int it = 0;
	
	while (resid > res_thresh && it < max_iters){
		// d = b - R*xk;
		x += diagonal_solve(A, r);

		it++;
		r = b - A*x;
		resid = norm_2(r);
	}

	return it;
}


// gauss-seidel iteration
// for square diagonally dominant matrices
unsigned int gauss_seidel(const Matrix & A, const Vector & b, Vector & x, unsigned int max_iters, double res_thresh=1.0e-15){

	// decompose the matrix A
	Matrix U = strictly_upper(A);
	Matrix L = A-U;
	// Vector xk(b); xk.fill(0);

	double resid = 1.0;
	unsigned int it = 0;
	Vector r;
	Vector d;
	while (resid > res_thresh && it < max_iters){
		d = b - U*x;
		x = lower_triangular_solve(L, d);

		it++;
		r = b - A*x;
		resid = norm_2(r);
	}

	return it;
}


// sparse gauss-seidel iteration
// for square diagonally dominant matrices
unsigned int gauss_seidel(const SparseMatrix & A, const Vector & b, Vector & x, unsigned int max_iters, double res_thresh=1.0e-15){

	// decompose the matrix A
	SparseMatrix U = strictly_upper(A);
	SparseMatrix L = strictly_lower(A) + spdiag(diag(A));

	double resid = 1.0;
	unsigned int it = 0;
	Vector r;
	Vector d;
	while (resid > res_thresh && it < max_iters){
		d = b - U*x;
		x = lower_triangular_solve(L, d);

		it++;
		r = b - A*x;
		resid = norm_2(r);
	}

	return it;
}


// successive over-relaxation iteration
// for square diagonally dominant matrices
unsigned int sor(const Matrix & A, const Vector & b, Vector & x, double w, unsigned int max_iters, double res_thresh=1.0e-15){

	// decompose the matrix A
	Matrix U = strictly_upper(A);
	Vector dg = diag(A);
	Matrix D = diag(dg);
	Matrix L = A-U-D;

	double resid = 1.0;
	unsigned int it = 0;
	Vector r;
	Vector d;
	Matrix N;
	while (resid > res_thresh && it < max_iters){
		d = w*b - (w*U+(w-1.0)*D)*x;
		N = D + w*L;
		x = lower_triangular_solve(N, d);

		it++;
		r = b - A*x;
		resid = norm_2(r);
	}

	return it;
}

// sparse successive over-relaxation iteration
// for square diagonally dominant matrices
unsigned int sor(const SparseMatrix & A, const Vector & b, Vector & x, double w, unsigned int max_iters, double res_thresh=1.0e-15){

	// decompose the matrix A
	SparseMatrix U = strictly_upper(A);
	SparseMatrix D = spdiag(diag(A));
	SparseMatrix L = strictly_lower(A);

	double resid = 1.0;
	unsigned int it = 0;
	Vector r;
	Vector d;
	SparseMatrix N = D + w*L;
	while (resid > res_thresh && it < max_iters){
		d = w*b - w*(U*x)-(w-1.0)*(D*x);
		//N = D + w*L;
		x = lower_triangular_solve(N, d);

		it++;
		r = b - A*x;
		resid = norm_2(r);
	}

	return it;
}

// steepest descent
// for square, SPD matrices only
// uses the x argument as x0
unsigned int steepest_descent(const Matrix & mtx, const Vector & b, Vector & x, unsigned int max_iters, double res_thresh=1.0e-15)
{

	Vector rv = b-mtx*x;
	double resid = rv.norm();
	double alpha;
	Vector matvec;

	unsigned int ctr=0;
	while (resid > res_thresh && ctr < max_iters)
	{
		// one matrix-vector product
		matvec = mtx*rv;

		// direction is residual
		alpha = Vector::dot(rv,rv)/Vector::dot(rv,matvec);

		// update x
		x += alpha*rv;

		// recalculate residual
		rv = b - mtx*x;

		// if (ctr%50 == 0)
		// {
		// 	// recalculate the exact resid every 50 iters
		// 	rv = b - mtx*x;
		// }
		// else
		// {
		// 	// update residual
		// 	rv -= alpha*matvec;
		// }

		resid = rv.norm();
		//std::cout << "residual: " << resid << std::endl;

		ctr++;
	}

	return ctr;
}


// conjugate gradient
// for square, SPD matrices only
// uses the x argument as x0
unsigned int conjugate_gradient(const Matrix & mtx, const Vector & b, Vector & x, unsigned int max_iters, double res_thresh=1.0e-15)
{
	Vector rv = b-mtx*x;

	double resid = rv.norm();

	// first iteration here
	Vector d = rv;
	Vector matvec = mtx*d;
	double alpha = Vector::dot(d,rv)/Vector::dot(d,matvec);
	x = alpha*b;

	// continue iterating
	unsigned int ctr=0;
	while (resid > res_thresh && ctr < max_iters)
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

	return ctr;
}


// sparse conjugate gradient
// for square, SPD matrices only
// uses the x argument as x0
unsigned int conjugate_gradient(const SparseMatrix & mtx, const Vector & b, Vector & x, unsigned int max_iters, double res_thresh=1.0e-15)
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
	while (resid > res_thresh && ctr < max_iters)
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

	// std::cout << "iterated: " << ctr << " times" << std::endl;

	return ctr;
}


// preconditioned sparse conjugate gradient
// for square, SPD matrices only
// uses the x argument as x0
unsigned int conjugate_gradient(const Preconditioner * pc, const SparseMatrix & A, const Vector & b, Vector & x, unsigned int max_iters, double res_thresh=1.0e-15)
{
	Vector r = b-A*x;
	double resid = r.norm();

	// setup
	Vector z = pc->solve(r);
	Vector d = 1*z;
	Vector matvec;
	double alpha, beta, zdr;

	// continue iterating
	unsigned int it=0;
	while (resid > res_thresh && it < max_iters)
	{

		// one matrix-vector product
		matvec = A*d;

		zdr = Vector::dot(z,r);
		alpha = zdr/Vector::dot(d,matvec);

		// update x
		x += alpha*d;

		// update residual
		r -= alpha*matvec;
		resid = r.norm();

		// update z
		z = pc->solve(r);

		// update beta
		beta = Vector::dot(z,r)/zdr;

		// new direction
		d = z + beta*d;

		it++;
	}

	// std::cout << "iterated: " << it << " times" << std::endl;

	return it;
}



// conjugate residual
// for square, symmetric matrices only
// uses the x argument as x0
unsigned int conjugate_residual(const Matrix & A, const Vector & b, Vector & x, unsigned int max_iters, double res_thresh=1.0e-15)
{
	Vector r = b-A*x;
	Vector rold;
	Vector Arold;
	Vector Ar = A*r;
	Vector Ap = Ar;
	double resid = r.norm();

	// first iteration here
	Vector d = r;
	double alpha, beta;

	// continue iterating
	unsigned int it=0;
	while (resid > res_thresh && it < max_iters)
	{
		// calc alpha
		Ap = A*d;
		alpha = Vector::dot(Ar,r)/Vector::dot(Ap,Ap);

		// update x
		x += alpha*d;

		// calculate residual
		rold = 1*r;
		r -= alpha*Ap;
		resid = r.norm();

		Arold = 1*Ar;
		Ar = A*r;
		beta = Vector::dot(r,Ar)/Vector::dot(rold,Arold);

		// update direction
		d = r+beta*d;
		Ap = Ar+beta*Ap;
		
		it++;
	}

	// std::cout << "iterated: " << it << " times" << std::endl;

	return it;
}


// sparse conjugate residual
// for square, symmetric matrices only
// uses the x argument as x0
unsigned int conjugate_residual(const SparseMatrix & A, const Vector & b, Vector & x, unsigned int max_iters, double res_thresh=1.0e-15)
{
	Vector r = b-A*x;
	Vector rold;
	Vector Arold;
	Vector Ar = A*r;
	Vector Ap = Ar;
	double resid = r.norm();

	// first iteration here
	Vector d = r;
	double alpha, beta;

	// continue iterating
	unsigned int it=0;
	while (resid > res_thresh && it < max_iters)
	{
		// calc alpha
		Ap = A*d;
		alpha = Vector::dot(Ar,r)/Vector::dot(Ap,Ap);

		// update x
		x += alpha*d;

		// calculate residual
		rold = 1*r;
		r -= alpha*Ap;
		resid = r.norm();

		Arold = 1*Ar;
		Ar = A*r;
		beta = Vector::dot(r,Ar)/Vector::dot(rold,Arold);

		// update direction
		d = r+beta*d;
		Ap = Ar+beta*Ap;
		
		it++;
	}

	// std::cout << "iterated: " << it << " times" << std::endl;

	return it;
}

/*
// preconditioned sparse conjugate residual
// for square, symmetric matrices only
// uses the x argument as x0
void conjugate_residual(const Preconditioner * pc, const SparseMatrix & A, const Vector & b, Vector & x, unsigned int max_iters, double res_thresh=1.0e-15)
{
	Vector r = b-A*x;
	Vector rhat = pc->solve(r);
	Vector rold;
	Vector Arold;
	Vector Ar = A*r;
	Vector Ap = 1*Ar;
	Vector dhat;
	double resid = r.norm();

	// first iteration here
	Vector d = r;
	double alpha, beta;

	// continue iterating
	unsigned int it=0;
	while (resid > res_thresh && it < max_iters)
	{
		// calc alpha
		dhat = pc->solve(d);
		Ap = A*dhat;
		alpha = Vector::dot(Ar,r)/Vector::dot(Ap,Ap);

		// update x
		x += alpha*dhat;

		// calculate residual
		rold = 1*r;
		r -= alpha*Ap;
		resid = r.norm();

		Arold = 1*Ar;
		rhat = pc->solve(r);
		Ar = A*rhat;
		beta = Vector::dot(r,Ar)/Vector::dot(rold,Arold);

		// update direction
		d = r+beta*d;
		Ap = Ar+beta*Ap;
		
		it++;
	}

	std::cout << "iterated: " << it << " times" << std::endl;

	return;
}
*/

// BiConjugate Gradient method
unsigned int bicg(const Matrix & A, const Vector & b, Vector & x, unsigned int max_iters, double res_thresh=1.0e-15){
	// initialize stuff
	Matrix At = ~A;
	Vector r = b - A*x;
	Vector bc = b; 
	Vector rhat = ~bc - (~x)*At; rhat.transpose();
	Vector xhat = ~x;
	double resid = 1.0;
	unsigned int it = 0;

	Vector p = r;
	Vector phat = rhat; 
	Vector Ap;
	double alpha, beta, rho;


	while (resid > res_thresh && it < max_iters){

		Ap = A*p;
		rho = Vector::dot(rhat, r);
		alpha = rho/Vector::dot(phat, Ap);
		x += alpha*p;
		xhat += alpha*(phat);
		r -= alpha*Ap;
		// std::cout << At.rows() << ", " << At.cols() << std::endl;
		// std::cout << rhat.rows() << ", " << rhat.cols() << std::endl;
		// std::cout << phat.rows() << ", " << phat.cols() << std::endl;
		rhat = rhat - alpha*phat*At; rhat.transpose();

		beta = Vector::dot(rhat, r)/rho;
		p = r+beta*p;
		phat = rhat + beta*phat; phat.transpose();

		// update counters
		it++;
		resid = norm_2(r);
	}

	// std::cout << "iterated: " << it << " times" << std::endl;
	return it;
}


// BiConjugate Residual method
unsigned int bicr(const Matrix & A, const Vector & b, Vector & x, unsigned int max_iters, double res_thresh=1.0e-15){
	// initialize stuff
	Matrix At = ~A;
	Vector r = b - A*x;
	Vector bc = b; 
	Vector rhat = 1*r;
	Vector p =  0*r;
	Vector phat = 0*r;

	double resid = 1.0;
	unsigned int it = 0;

	Vector Ap=A*p, Atphat;
	Vector Ar = A*r;
	double alpha, beta=0;


	while (resid > res_thresh && it < max_iters){

		p = r+beta*p;
		phat = rhat+beta*phat;
		Ap = Ar + beta*Ap;
		Atphat = At*phat;

		alpha = Vector::dot(r,Ar)/Vector::dot(Atphat, Ap);

		x += alpha*p;
		r -= alpha*Ap;
		rhat -= alpha*Atphat;

		beta = Vector::dot(rhat,A*r)/Vector::dot(rhat,Ar);
		Ar = A*r;

		// update counters
		it++;
		resid = norm_2(r);
	}

	// std::cout << "iterated: " << it << " times" << std::endl;
	return it;
}


// BiConjugate Gradient Stabilized method
unsigned int bicgstab(const Matrix & A, const Vector & b, Vector & x, unsigned int max_iters, double res_thresh=1.0e-15){
	// initialize stuff
	Vector r = b - A*x;
	double resid = 1.0;
	unsigned int it = 0;
	Vector rhat = r;
	double rho=1, alpha=1, omega=1;
	Vector v(r); v.fill(0);
	Vector p(r); p.fill(0);
	Vector s,t;
	double beta, rho_old;

	while (resid > res_thresh && it < max_iters){

		// enforce biorthogonality and biconjugacy
		rho_old = rho;
		rho = Vector::dot(rhat, r);
		beta = rho/rho_old*(alpha/omega);
		p = r + beta*(p-omega*v);
		v = A*p;
		alpha = rho/Vector::dot(rhat,v);
		s = r - alpha*v;
		t = A*s;
		omega = Vector::dot(s,t)/Vector::dot(t,t);
		x += alpha*p + omega*s;

		// recalculate residual
		r = s - omega*t;

		// update counters
		it++;
		resid = norm_2(r);
	}

	// std::cout << "iterated: " << it << " times" << std::endl;
	return it;
}


// sparse BiConjugate Gradient Stabilized method
unsigned int bicgstab(const SparseMatrix & A, const Vector & b, Vector & x, unsigned int max_iters, double res_thresh=1.0e-15){
	// initialize stuff
	Vector r = b - A*x;
	double resid = 1.0;
	unsigned int it = 0;
	Vector rhat = r;
	double rho=1, alpha=1, omega=1;
	Vector v(r); v.fill(0);
	Vector p(r); p.fill(0);
	Vector s,t;
	double beta, rho_old;

	while (resid > res_thresh && it < max_iters){

		// enforce biorthogonality and biconjugacy
		rho_old = rho;
		rho = Vector::dot(rhat, r);
		beta = rho/rho_old*(alpha/omega);
		p = r + beta*(p-omega*v);
		v = A*p;
		alpha = rho/Vector::dot(rhat,v);
		s = r - alpha*v;
		t = A*s;
		omega = Vector::dot(s,t)/Vector::dot(t,t);
		x += alpha*p + omega*s;

		// recalculate residual
		r = s - omega*t;

		// update counters
		it++;
		resid = norm_2(r);
	}

	// std::cout << "iterated: " << it << " times" << std::endl;
	return it;
}


// preconditioned sparse BiConjugate Gradient Stabilized method
unsigned int bicgstab(const Preconditioner * pc, const SparseMatrix & A, const Vector & b, Vector & x, unsigned int max_iters, double res_thresh=1.0e-15){
	// initialize stuff
	Vector r = b - A*x;
	double resid = 1.0;
	unsigned int it = 0;
	Vector rhat = 1*r;
	double rho=1, alpha=1, omega=1;
	Vector v(r); v.fill(0);
	Vector p(r); p.fill(0);
	Vector s,t, shat, phat;
	double beta, rho_old=1;
	// Vector r0hat = 1*r;

	while (resid > res_thresh && it < max_iters){


		// enforce biorthogonality and biconjugacy
		
		rho = Vector::dot(rhat, r);
		beta = rho/rho_old*(alpha/omega);
		p = r + beta*(p-omega*v);
		phat = pc->solve(p);
		v = A*phat;
		alpha = rho/Vector::dot(rhat,v);
		s = r - alpha*v;
		shat = pc->solve(s);
		t = A*shat;
		omega = Vector::dot(s,t)/Vector::dot(t,t);
		x += alpha*phat + omega*shat;

		// recalculate residual
		r = s - omega*t;

		rho_old = rho;

		// update counters
		it++;
		resid = norm_2(r);
	}

	// std::cout << "iterated: " << it << " times" << std::endl;
	return it;
}


// Restarted Generalized Minimal Residual Method (GMRES)
// restarts every k iterations
unsigned int gmres_k(const Matrix & A, const Vector & b, Vector & x, unsigned int k, unsigned int max_iters, double res_thresh=1.0e-15){
	// initialize stuff
	k = std::min(std::size_t(k),A.cols());
	Vector r, p, pr, w, c(k), s(k), x0=1*x, yr, y(k); y.fill(0);
	Vector e1(k+1); e1.fill(0); e1(0)=1;
	Matrix Q(b.length(), k+1); Q.fill(0);
	Matrix Hr;
	Matrix H(k+1,k); 
	unsigned int red_k = k-1;
	double beta, gamma;

	double resid = 1.0;
	unsigned int it = 0;

	while (resid > res_thresh && it < max_iters ){
		r = b - A*x0;
		beta = norm_2(r);
		resid = beta;
		Q.fill(0);
		Q.col(0) = r/beta;
		H.fill(0);
		
		// initialize the rhs and the solution vector
		p = beta*e1;
		y.fill(0);

		// initialize givens rotators
		s.fill(0);
		c.fill(1);

		// inner loop
		for (auto j=0; j<k; j++){
			// if (it >= max_iters || it >= x.length() || resid <= res_thresh) break;
			if (it >= max_iters || resid <= res_thresh) break;

			// std::cout << "inner loop iter: " << j << std::endl;

			// Part of the Arnoldi iteration!
			// matvec
			w = A*Q.col(j);

			// gram-schmidt process
			for (auto i=0; i<=j; i++){
				H(i,j) = Vector::dot(w,Q.col(i));
				w -= H(i,j)*Q.col(i);
			}
			H(j+1,j) = norm_2(w);
			Q.col(j+1) = w/H(j+1,j);

			// solve reduced system least squares problem
			// apply previous givens rotation acting on H
			double tmp;
			for (auto i=0; i<j; i++){
				tmp = c(i)*H(i,j)+s(i)*H(i+1,j);
				H(i+1,j) = -s(i)*H(i,j)+c(i)*H(i+1,j);
				H(i,j) = tmp;
			}
			// now compute the jth rotation to annihilate the (j+1,j) entry
			gamma = sqrt(H(j,j)*H(j,j)+H(j+1,j)*H(j+1,j));
			c(j) = H(j,j)/gamma;
			s(j) = H(j+1,j)/gamma;

			// apply current givens rotation on H
			H(j,j) = gamma;
			H(j+1,j) = 0;

			// apply current givens rotation on p
			p(j+1) = -s(j)*p(j);
			p(j) = c(j)*p(j);

			it++;
			
			resid = fabs(p(j+1));
			if (resid < res_thresh){
				red_k = j;
				break;
			}
		}

		// std::cout << "red_k: " << red_k << " k: " << k << std::endl;


		// reduced system matrix
		Hr = H(0,red_k,0,red_k);

		// do back substitution to get y
		pr = 1*p.subcol(0,0,red_k);
		yr = upper_triangular_solve(Hr, pr);
		y.subcol(0,0,red_k) = 1*yr;

		// form approximate solution
		x = 1*x0;
		for (auto i=0; i<=red_k; i++) x += y(i)*Q.col(i);

		// set initial value as current x value
		x0 = 1*x;
	}

	// std::cout << "iterated: " << it << " times" << std::endl;
	return it;
}


// Sparse Restarted Generalized Minimal Residual Method (GMRES)
// restarts every k iterations
unsigned int gmres_k(const SparseMatrix & A, const Vector & b, Vector & x, unsigned int k, unsigned int max_iters, double res_thresh=1.0e-15){
	// initialize stuff
	k = std::min(std::size_t(k),A.cols());
	Vector r, p, pr, w, c(k), s(k), x0=1*x, yr, y(k); y.fill(0);
	Vector e1(k+1); e1.fill(0); e1(0)=1;
	Matrix Q(b.length(), k+1); Q.fill(0);
	Matrix Hr;
	Matrix H(k+1,k); 
	unsigned int red_k = k-1;
	double beta, gamma;

	double resid = 1.0;
	unsigned int it = 0;

	while (resid > res_thresh && it < max_iters ){
		r = b - A*x0;
		beta = norm_2(r);
		resid = beta;
		Q.fill(0);
		Q.col(0) = r/beta;
		H.fill(0);
		
		// initialize the rhs and the solution vector
		p = beta*e1;
		y.fill(0);

		// initialize givens rotators
		s.fill(0);
		c.fill(1);

		// inner loop
		for (auto j=0; j<k; j++){
			// if (it >= max_iters || it >= x.length() || resid <= res_thresh) break;
			if (it >= max_iters || resid <= res_thresh) break;

			// std::cout << "inner loop iter: " << j << std::endl;

			// Part of the Arnoldi iteration!
			// matvec
			w = A*Q.col(j);

			// gram-schmidt process
			for (auto i=0; i<=j; i++){
				H(i,j) = Vector::dot(w,Q.col(i));
				w -= H(i,j)*Q.col(i);
			}
			H(j+1,j) = norm_2(w);
			Q.col(j+1) = w/H(j+1,j);

			// solve reduced system least squares problem
			// apply previous givens rotation acting on H
			double tmp;
			for (auto i=0; i<j; i++){
				tmp = c(i)*H(i,j)+s(i)*H(i+1,j);
				H(i+1,j) = -s(i)*H(i,j)+c(i)*H(i+1,j);
				H(i,j) = tmp;
			}
			// now compute the jth rotation to annihilate the (j+1,j) entry
			gamma = sqrt(H(j,j)*H(j,j)+H(j+1,j)*H(j+1,j));
			c(j) = H(j,j)/gamma;
			s(j) = H(j+1,j)/gamma;

			// apply current givens rotation on H
			H(j,j) = gamma;
			H(j+1,j) = 0;

			// apply current givens rotation on p
			p(j+1) = -s(j)*p(j);
			p(j) = c(j)*p(j);

			it++;
			
			resid = fabs(p(j+1));
			if (resid < res_thresh){
				red_k = j;
				break;
			}
		}

		// std::cout << "red_k: " << red_k << " k: " << k << std::endl;


		// reduced system matrix
		Hr = H(0,red_k,0,red_k);

		// do back substitution to get y
		pr = 1*p.subcol(0,0,red_k);
		yr = upper_triangular_solve(Hr, pr);
		y.subcol(0,0,red_k) = 1*yr;

		// form approximate solution
		x = 1*x0;
		for (auto i=0; i<=red_k; i++) x += y(i)*Q.col(i);

		// set initial value as current x value
		x0 = 1*x;
	}

	// std::cout << "iterated: " << it << " times" << std::endl;
	return it;
}


// Preconditioned Sparse Restarted Generalized Minimal Residual Method (GMRES)
// restarts every k iterations
unsigned int gmres_k(const Preconditioner * pc, const SparseMatrix & A, const Vector & b, Vector & x, unsigned int k, unsigned int max_iters, double res_thresh=1.0e-15){
	// initialize stuff
	k = std::min(std::size_t(k),A.cols());
	Vector r, p, pr, w, c(k), s(k), x0=1*x, yr, y(k); y.fill(0);
	Vector e1(k+1); e1.fill(0); e1(0)=1;
	Matrix Q(b.length(), k+1); Q.fill(0);
	Matrix Hr;
	Matrix H(k+1,k); 
	unsigned int red_k = k-1;
	double beta, gamma;

	double resid = 1.0;
	unsigned int it = 0;

	while (resid > res_thresh && it < max_iters ){
		r = pc->solve(b - A*x0);
		beta = norm_2(r);
		resid = beta;
		Q.fill(0);
		Q.col(0) = r/beta;
		H.fill(0);
		
		// initialize the rhs and the solution vector
		p = beta*e1;
		y.fill(0);

		// initialize givens rotators
		s.fill(0);
		c.fill(1);

		// inner loop
		for (auto j=0; j<k; j++){
			// if (it >= max_iters || it >= x.length() || resid <= res_thresh) break;
			if (it >= max_iters || resid <= res_thresh) break;

			// std::cout << "inner loop iter: " << j << std::endl;

			// Part of the Arnoldi iteration!
			// matvec
			w = pc->solve(A*Q.col(j));

			// gram-schmidt process
			for (auto i=0; i<=j; i++){
				H(i,j) = Vector::dot(w,Q.col(i));
				w -= H(i,j)*Q.col(i);
			}
			H(j+1,j) = norm_2(w);
			Q.col(j+1) = w/H(j+1,j);

			// solve reduced system least squares problem
			// apply previous givens rotation acting on H
			double tmp;
			for (auto i=0; i<j; i++){
				tmp = c(i)*H(i,j)+s(i)*H(i+1,j);
				H(i+1,j) = -s(i)*H(i,j)+c(i)*H(i+1,j);
				H(i,j) = tmp;
			}
			// now compute the jth rotation to annihilate the (j+1,j) entry
			gamma = sqrt(H(j,j)*H(j,j)+H(j+1,j)*H(j+1,j));
			c(j) = H(j,j)/gamma;
			s(j) = H(j+1,j)/gamma;

			// apply current givens rotation on H
			H(j,j) = gamma;
			H(j+1,j) = 0;

			// apply current givens rotation on p
			p(j+1) = -s(j)*p(j);
			p(j) = c(j)*p(j);

			it++;
			
			resid = fabs(p(j+1));
			if (resid < res_thresh){
				red_k = j;
				break;
			}
		}

		// std::cout << "red_k: " << red_k << " k: " << k << std::endl;


		// reduced system matrix
		Hr = H(0,red_k,0,red_k);

		// do back substitution to get y
		pr = 1*p.subcol(0,0,red_k);
		yr = upper_triangular_solve(Hr, pr);
		y.subcol(0,0,red_k) = 1*yr;

		// form approximate solution
		x = 1*x0;
		for (auto i=0; i<=red_k; i++) x += y(i)*Q.col(i);

		// set initial value as current x value
		x0 = 1*x;
	}

	// std::cout << "iterated: " << it << " times" << std::endl;
	return it;
}



inline void get_strong_neighbors(const SparseMatrix & A, unsigned int row, std::set<unsigned int> & S, double theta){
	auto rowptr = A.row_ptr();
	auto data = A.data();
	std::set<unsigned int> Spts;
	auto it = rowptr[row];

	// for this row, find the max negative value
	double rowmax = 0.0;
	while (it != rowptr[row+1]){
		unsigned int j = it->first;
		if (j==row){
			it++;
			continue;
		}

		rowmax = std::max(rowmax, -it->second);
		it++;
	}

	// A point i is strongly connected to j if -Aij >= theta*max(-Aik)
	it = rowptr[row];
	while (it != rowptr[row+1]){
		unsigned int j = it->first;
		if (j==row){
			it++;
			continue;
		}

		// these points are strongly connected
		if (-it->second >= theta*rowmax){
			Spts.insert(it->first);
		}
		it++;
	}

	std::swap(Spts, S);
}


// get coarse point neighbor indices for point "row"
// support function for AMG
inline void get_C_neighbors(const SparseMatrix & A, const std::set<unsigned int> & C, 
				  unsigned int row, std::set<unsigned int> & rowC){
	auto rowptr = A.row_ptr();
	auto data = A.data();
	std::set<unsigned int> Cpts;
	auto it = rowptr[row];

	// for this row, find the max negative value
	double rowmax = 0.0;
	while (it != rowptr[row+1]){
		unsigned int j = it->first;

		if (C.find(j) != C.end()) Cpts.insert(j);
		it++;
	}

	std::swap(Cpts, rowC);
}

// get fine point neighbor indices for point "row"
// support function for AMG
inline void get_F_neighbors(const SparseMatrix & A, const std::set<unsigned int> & C, 
				  unsigned int row, std::set<unsigned int> & rowF){
	auto rowptr = A.row_ptr();
	auto data = A.data();
	std::set<unsigned int> Fpts;
	auto it = rowptr[row];

	// for this row, find the max negative value
	double rowmax = 0.0;
	while (it != rowptr[row+1]){
		unsigned int j = it->first;

		if (C.find(j) == C.end()) Fpts.insert(j);
		it++;
	}

	std::swap(Fpts, rowF);
}

// standard 2-pass coarsening
// support function for AMG
void amg_standard_coarsening(const SparseMatrix & A, std::set<unsigned int> & C, std::set<unsigned int> & F, double theta)
{


	// some preliminaries
	auto rowptr = A.row_ptr();
	auto data = A.data();
	
	// get vector of row maxima
	// std::cout << "...Constructing row maxima..." << std::endl;
	std::vector<double> rowmax(A.rows(), 0);
	for (auto i=0; i<A.rows(); i++){
		auto it = rowptr[i];
		// for this row, find the max negative valu
		while (it != rowptr[i+1]){
			if (it->first==i){
				it++;
				continue;
			}

			rowmax[i] = std::max(rowmax[i], -it->second);
			it++;
		}
		rowmax[i] *= -theta;
	}

	// get vector boolean of which points are coarse
	std::vector<bool> iscoarse(A.rows(), false);

	// identify point in C and F where
	// C is the set of coarse grid points, and F is 
	// the set of fine grid points
	std::vector<unsigned int> lambda(A.cols());
	std::set<unsigned int> Cpts, Fpts, Upts, S;
	// std::cout << "...Constructing lambda measure..." << std::endl;
	for (auto i=0; i<A.rows(); i++){
		//std::cout << i << std::endl;

		auto it = rowptr[i];
		// for this row, find strongly connected points
		while (it != rowptr[i+1]){
			if (it->first==i){
				it++;
				continue;
			}

			if (it->second <= rowmax[i]) lambda[it->first]++;
			it++;
		}
	}
	// std::cout << "DONE" << std::endl;
	// std::cout << "Constructed Lambda..." << std::endl;
	// for (auto i=lambda.begin(); i!=lambda.end(); i++){
	// 	std::cout << "lambda: " << *i << std::endl;
	// }

	// first pass separates C and F points
	// std::cout << "Starting first pass..." << std::endl ;
	for (auto i=0; i<A.cols(); i++) Upts.insert(i);
	unsigned int imax = std::distance(lambda.begin(), std::max_element(lambda.begin(), lambda.end()));//.begin(), lambda.end()));
	while (lambda[imax] > 0 && Upts.size() > 0){
		// std::cout << "imax: " << imax << " value: " << lambda[imax] << std::endl;

		// set lambda to zero
		lambda[imax] = 0;
		Cpts.insert(imax);
		Upts.erase(imax);
		iscoarse[imax] = true;

		// std::cout << "inserted coarse point..." << imax << std::endl;

		// get the strong connections to this point... these become fine points
		auto it = rowptr[imax];
		// A point i is strongly connected to j if -Aij >= theta*max(-Aik)
		while (it != rowptr[imax+1]){
			unsigned int j = it->first;
			if (j==imax){
				it++;
				continue;
			}

			// these points are strongly connected
			if (it->second <= rowmax[imax] && lambda[j] > 0){
				lambda[j] = 0;
				Fpts.insert(j);
				Upts.erase(j);

				// std::cout << "inserted fine point..." << j << std::endl;

				// get points strongly connected to these points
				// and increment lambda for them
				auto itj = rowptr[j];
				while (itj != rowptr[j+1]){
					unsigned int k = itj->first;
					if (k==j){
						itj++;
						continue;
					}

					// these points are strongly connected
					if (itj->second <= rowmax[j] && lambda[k] > 0){
						lambda[k] += 1;
					}
					itj++;
				}
			}
			it++;
		}

		// get new max lambda point to become coarse point
		// imax = std::distance(std::max_element(lambda.begin(), lambda.end()));
		imax = std::distance(lambda.begin(), std::max_element(lambda.begin(), lambda.end()));

		// std::cout << "imax: " << imax << " value: " << lambda[imax] << std::endl;
		// std::cout << "Csize: " << Cpts.size() << std::endl;
		// std::cout << "Fsize: " << Fpts.size() << std::endl;
		// std::cout << "Usize: " << Upts.size() << std::endl;

	}
	// std::cout << "DONE" << std::endl;
	// std::cout << "Csize: " << Cpts.size() << std::endl;
	// std::cout << "Fsize: " << Fpts.size() << std::endl;
	// std::cout << "Usize: " << Upts.size() << std::endl;



	// second pass checks that strong F-F dependencies share a common C point
	auto itF=Fpts.begin(); 
	std::set<unsigned int> iC;
	// std::cout << "Starting second pass..." << std::endl ;
	std::vector<bool> Cneighb(A.rows(), false);
	while (itF!=Fpts.end()){
		unsigned int i=*itF;
		//std::cout << "Fine point: " << i << std::endl;

		// mark C neighbors for this row
		auto it = rowptr[i];
		while (it != rowptr[i+1]){
			// unsigned int j = it->first;
			if (it->first==i){
				it++;
				continue;
			}
			if (iscoarse[it->first]) Cneighb[it->first] = true;
			it++;
		}
		
		// for strongly dependent fine neighbors, check for a common C point
		// iterate over strongly dependent fine neighbors
		it = rowptr[i];
		while (it != rowptr[i+1]){
			unsigned int j = it->first;
			if (j==i){
				it++;
				continue;
			}

			// these neighbors are strong fine points
			if (it->second <= rowmax[i] && !iscoarse[j]){
				unsigned int cc=0;

				// get C neighbors for row j and
				// find how many neighbors are matching
				auto itj = rowptr[j];
				while (itj != rowptr[j+1]){
					if (itj->first == j){
						itj++;
						continue;
					}
					if (iscoarse[itj->first] && Cneighb[itj->first]) cc++;
					itj++;
				}

				// if there are 0 matching C neighbors, change from F to C point
				if (cc==0){
					Cpts.insert(*itF);
					iscoarse[*itF] = true;
					Fpts.erase(*itF);
					break;
				}
			}
			it++;
		}

		// unmark C neighbors for this row
		it = rowptr[i];
		while (it != rowptr[i+1]){
			if (it->first==i){
				it++;
				continue;
			}
			if (iscoarse[it->first]) Cneighb[it->first] = false;
			it++;
		}

		itF++;
	}
	// std::cout << "DONE" << std::endl;
	// std::cout << "Csize: " << Cpts.size() << std::endl;
	// std::cout << "Fsize: " << Fpts.size() << std::endl;
	// std::cout << "Usize: " << Upts.size() << std::endl;

	std::swap(Cpts,C);
	std::swap(Fpts,F);
}


// direct interpolation (simplest)
// support function for AMG
void amg_direct_interpolation(const SparseMatrix & A, 
							  const std::set<unsigned int> & C, const std::set<unsigned int> & F, 
							  SparseMatrix & W)
{

	// setup vector that identifies coarse points
	std::vector<bool> iscoarse(A.cols(), false);
	for (auto it=C.begin(); it!=C.end(); it++) iscoarse[*it] = true;

	// setup a vector that identifies the index of each coarse
	// point in the restricted vector
	std::vector<unsigned int> coarse_ordering(A.cols(), 0);
	unsigned int ct=0;
	for (auto it=C.begin(); it!=C.end(); it++) coarse_ordering[*it] = ct++;

	// determine interpolation matrix W
	// and its transpose (the restriction matrix) Wt
	SparseMatrix Wi(A.rows(), C.size());
	auto rowptr = A.row_ptr();
	auto data = A.data();
	for (auto i=0; i<A.rows(); i++){
		// if the ith point is coarse, set the interpolation to return unity
		if (iscoarse[i]){
			Wi.set(i,coarse_ordering[i], 1.0);
			continue;
		}
		auto it = rowptr[i];

		// determine weighting for fine points
		// get rowsum, coarsesum, and diag
		double rowsum=0;
		double corsum=0;
		double rowdiag=0;
		while (it != rowptr[i+1]){
			unsigned int j = it->first;
			if (j==i){
				rowdiag=it->second;	// the diagonal in the row
				it++;
				continue;
			}

			// add to sums
			rowsum += it->second;		// sum over whole row
			if (iscoarse[j]) corsum += it->second;	// sum over coarse points in the row

			it++;
		}

		// scale row sum
		double rowscale=0;
		it = rowptr[i];
		while (it != rowptr[i+1]){
			unsigned int j = it->first;

			// add to sums
			if (iscoarse[j]) rowscale += -rowsum*it->second/(corsum*rowdiag);

			it++;
		}

		// set the values in this row
		it = rowptr[i];
		while (it != rowptr[i+1]){
			unsigned int j = it->first;

			// add to sums
			if (iscoarse[j]) Wi.set(i,coarse_ordering[j],-rowsum*it->second/(corsum*rowdiag*rowscale));

			it++;
		}
	}

	swap(Wi,W);
}


// standard interpolation (better than direct)
// support function for AMG
void amg_standard_interpolation(const SparseMatrix & A, 
							  const std::set<unsigned int> & C, const std::set<unsigned int> & F, 
							  SparseMatrix & W)
{
	double theta=0.25; // THIS SHOULD BE A PARAMETER!!

	// setup vector that identifies coarse points
	std::vector<bool> iscoarse(A.cols(), false);
	for (auto it=C.begin(); it!=C.end(); it++) iscoarse[*it] = true;

	// setup a vector that identifies the index of each coarse
	// point in the restricted vector
	std::vector<unsigned int> coarse_ordering(A.cols(), 0);
	unsigned int ct=0;
	for (auto it=C.begin(); it!=C.end(); it++) coarse_ordering[*it] = ct++;

	// determine interpolation matrix W
	// and its transpose (the restriction matrix) Wt
	SparseMatrix Wi(A.rows(), C.size());
	auto rowptr = A.row_ptr();
	auto data = A.data();
	std::set<unsigned int> Ci, Fi, Si, FSi; 
	for (auto i=0; i<A.rows(); i++){
		// if the ith point is coarse, set the interpolation to return unity
		if (iscoarse[i]){
			Wi.set(i,coarse_ordering[i], 1.0);
			continue;
		}
		
		// set of fine neighbors
		get_F_neighbors(A, C, i, Fi);

		// set of coarse neighbors
		get_C_neighbors(A, C, i, Ci);

		// set of strongly influencing neighbors
		get_strong_neighbors(A, i, Si, theta);

		// get set of strongly influencing fine neighbors
		FSi.clear();
		std::set_intersection(Fi.begin(), Fi.end(), Si.begin(), Si.end(), std::inserter(FSi, FSi.begin()));

		// get sum over weakly influencing fine neighbors
		double denom=A.get(i,i);
		auto it = rowptr[i];
		while (it != rowptr[i+1]){
			if (!iscoarse[it->first] && Si.find(it->first)==Si.end()) denom += it->second;
			it++;
		}

		// can precalculate sum(A_mk) for k in Ci for each m in FSi
		std::map<unsigned int, double> sigmaAmk;
		for (auto itm=FSi.begin(); itm!=FSi.end(); itm++){
			double val=0;
			unsigned int m=*itm;

			// get sum of intersecting coarse neighbors for m
			auto im=rowptr[m];
			while (im!=rowptr[m+1]){
				unsigned int k=im->first;
				if (iscoarse[k]){
					if (Ci.find(k)!=Ci.end()) val+=im->second;
				}
				im++;
			}

			// set the value of the sum
			sigmaAmk[m] = val;
		}

		// iterate over existing columns j in row i
		it = rowptr[i];
		while (it != rowptr[i+1]){
			unsigned int j = it->first;
			double Aij = it->second;

			// sum over the set m in FSi
			double sigmam=0;
			for (auto itm=FSi.begin(); itm!=FSi.end(); itm++){
				unsigned int m=*itm;

				sigmam += A.get(i,m)*A.get(m,j)/sigmaAmk[m];
			}

			// set the values in the interpolation matrix
			Wi.set(i, coarse_ordering[j], -(Aij+sigmam)/denom);

			it++;
		}
	}
	// std::cout << "Csize: " << C.size() << std::endl;
	//std::cout << "FINISHED: " << Wi.nnz() << "/" << Wi.rows()*Wi.cols() << std::endl;


	swap(Wi,W);
}



void amg_setup(const SparseMatrix & A, std::vector<SparseMatrix *> & Ws, std::vector<SparseMatrix *> & As){

	std::cout << "Level: " << Ws.size() << std::endl;

	// if A is "small enough", then quit setup
	if (A.cols() < 40) return;

	// perform coarsening
	std::set<unsigned int> Cpts, Fpts;
	amg_standard_coarsening(A, Cpts, Fpts, 0.25);

	std::cout << "coarsened..." << std::endl;

	// get interpolation matrix from C-F splitting
	SparseMatrix * W = new SparseMatrix();
	amg_direct_interpolation(A, Cpts, Fpts, *W);
	// amg_standard_interpolation(A, Cpts, Fpts, *W);

	std::cout << "interpolated..." << std::endl;

	// generate restricted matrix
	SparseMatrix * Ar = new SparseMatrix();
	(*Ar) = W->Tmult(A*(*W));
	As.push_back(Ar);

	std::cout << "restricted --> " << Ar->rows() << " x " << Ar->cols() << std::endl;


	// W->mmwrite("SparseInterpolationMatrix.txt");
	// (W->densify()).dlmwrite("InterpolationMatrix.txt");
	// Ar->mmwrite("SparseOperator.txt");
	// (Ar->densify()).dlmwrite("Operator.txt");
	// throw -1;

	// push the interpolation matrix 
	Ws.push_back(W);

	// recurse on the restricted matrix
	amg_setup(*Ar, Ws, As);
}


// Algebraic Multigrid V cycle 
// v1 	- number of presmoothing steps
// v2	- number of postsmoothing steps
// Ws 	- interpolation matrices
// As 	- restricted operator matrices
void amgv(const SparseMatrix & A, const Vector & b, Vector & x, unsigned int level,
		  const std::vector<SparseMatrix *> & Ws, const std::vector<SparseMatrix *> & As, 
		  unsigned int v1=1, unsigned int v2=1){

	

	// if the system is small enough, solve directly
	if (A.cols() < 40){
		// std::cout << "L(" << level << "): direct solve" << std::endl ;
		Matrix Ad = A.densify();
		Matrix Q, R, U;
		qr_householder(Ad, U, R, Q);
		x = unitary_solve(Q, b);
		x = upper_triangular_solve(R, x);
		return;
	}


	// presmoothing
	// std::cout << "L(" << level << "): presmooth " ;
	// std::cout << " resid = " << norm_2(b-A*x) << std::endl;
	if (v1>0) gauss_seidel(A, b, x, v1);
	// std::cout << " resid = " << norm_2(b-A*x) << std::endl;

	// get residual
	Vector r = b - A*x;

	// restrict residual
	SparseMatrix * W = Ws[level];
	Vector rr = W->Tmult(r);

	// get restricted Operator A
	SparseMatrix * Ar = As[level];

	// solve restricted system for error
	Vector er(Ar->cols()); er.fill(0);
	amgv(*Ar, rr, er, level+1, Ws, As, v1, v2);

	// std::cout << "residual on the error: " << norm_2(rr-(*Ar)*er) << std::endl;

	// interpolate the error
	Vector e = (*W)*er;
	// std::cout << "interpolation matrix: " << (*W) << std::endl;
	// std::cout << "restricted error: " << er << std::endl;
	// std::cout << "interpolated error: " << e << std::endl;

	// std::cout << "residual on the interpolated error: " << norm_2(r-A*e) << std::endl;


	// correct for the error
	x += e;

	// postsmoothing
	// std::cout << "L(" << level << "): postsmooth " ;
	if (v2>0) gauss_seidel(A, b, x, v2);
	// std::cout << " resid = " << norm_2(b-A*x) << std::endl;

	// std::cout << "residual on the smoothed solution: " << norm_2(b-A*x) << std::endl;


	return;
}

// algebraic multigrid
unsigned int amg(const SparseMatrix & A, const Vector & b, Vector & x, unsigned int max_iters, double res_thresh=1.0e-15){

	// *********** SETUP PHASE **************
	std::vector<SparseMatrix *> Ws, As;
	amg_setup(A, Ws, As);
	// ********** END SETUP *******************

	// std::cout << "level: 0 " ;
	// std::cout << "A nnz/total: " << A.nnz() << "/" << A.rows()*A.cols() << std::endl;
	// for (auto i=0; i<Ws.size(); i++){
	// 	std::cout << "level: " << i+1 << " " ;
	// 	SparseMatrix & Ar = *As[i];
	// 	std::cout << "A_r nnz/total: " << Ar.nnz() << "/" << Ar.rows()*Ar.cols() << std::endl;
	// }


		
	// *********** SOLUTION PHASE **************
	unsigned int it=0;
	double resid = norm_2(b-A*x);
	while(resid > res_thresh && it < max_iters){
		std::cout << "iterating..." << it << std::endl;
		amgv(A, b, x, 0, Ws, As, 2, 2);
		resid = norm_2(b-A*x);
		it++;
	}
	// ********** END SOLUTION *******************

	// delete interpolation matrices and restricted operators
	for (auto i=0; i<Ws.size(); i++) delete Ws[i];
	for (auto i=0; i<As.size(); i++) delete As[i];

	return it;
	
}

// generalized complex eigenvalue decomposition

// schur decomposition

// estimate the condition number

// }
#endif