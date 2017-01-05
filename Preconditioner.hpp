#ifndef _PRECONDITIONER_H
#define _PRECONDITIONER_H

#include "Matrix.hpp"
#include "SparseMatrix.hpp"

Vector diagonal_solve(const SparseMatrix & A, const Vector & b);
Vector lower_triangular_solve(const SparseMatrix & A, const Vector & b);
Vector upper_triangular_solve(const SparseMatrix & A, const Vector & b);
Vector diag(const SparseMatrix & A);
SparseMatrix spdiag(const Vector & x);
SparseMatrix strictly_lower(const SparseMatrix & A);
SparseMatrix strictly_upper(const SparseMatrix & A);

class Preconditioner
{
public:
	virtual Vector solve(const Vector & b) const = 0;
};


#include "LinearSolvers.hpp"


// jacobi preconditioner
// P = D
class JacobiPreconditioner : public Preconditioner
{
public:
	JacobiPreconditioner(const SparseMatrix & A){
		m_A = &A;
	}

	Vector solve(const Vector & b) const{
		return diagonal_solve(*m_A, b);
	}
protected:

	const SparseMatrix * m_A;
};


// gauss-seidel preconditioner
// P = (D+L)
class GSPreconditioner : public Preconditioner
{
public:
	GSPreconditioner(const SparseMatrix & A){
		m_A = &A;
	}

	Vector solve(const Vector & b) const{
		return lower_triangular_solve(*m_A, b);
	}
protected:

	const SparseMatrix * m_A;
};



// symmetric gauss-seidel preconditioner
// P = (D+L)*inv(D)*(D+L)'
class SGSPreconditioner : public Preconditioner
{
public:
	SGSPreconditioner(const SparseMatrix & A){
		m_A = &A;
		SparseMatrix D = spdiag(diag(*m_A));
		m_D = new SparseMatrix(D);
	}

	~SGSPreconditioner(){delete m_D;};

	Vector solve(const Vector & b) const{
		Vector x = lower_triangular_solve(*m_A,b);
		Vector y = (*m_D)*x;
		return upper_triangular_solve(*m_A, y);
	}
protected:

	const SparseMatrix * m_A;
	SparseMatrix * m_D;
};


// successive overrelaxation preconditioner
// P = (D/w+L) = (D+L) + (1-w)/w *D
class SORPreconditioner : public Preconditioner
{
public:
	SORPreconditioner(const SparseMatrix & A, double w){
		m_A = &A;
		m_w = w;
		SparseMatrix P = spdiag(diag(*m_A))/w + strictly_lower(*m_A);
		m_P = new SparseMatrix(P);
	}

	~SORPreconditioner(){delete m_P;};

	Vector solve(const Vector & b) const{
		return lower_triangular_solve(*m_P, b);
	}
protected:

	const SparseMatrix * m_A;
	SparseMatrix * m_P;
	double m_w;
};



// symmetric gauss-seidel preconditioner
// P = w/(2-w) (D/w+L)*inv(D)*(D/w+L)'
class SSORPreconditioner : public Preconditioner
{
public:
	SSORPreconditioner(const SparseMatrix & A, double w){
		m_A = &A;
		m_w = w;
		SparseMatrix L = spdiag(diag(*m_A))/w + strictly_lower(*m_A);
		m_L = new SparseMatrix(L);
		SparseMatrix R = spdiag(diag(*m_A))/w + strictly_upper(*m_A);
		m_R = new SparseMatrix(R);
		SparseMatrix D = spdiag(diag(*m_A));
		m_D = new SparseMatrix(D);
	}

	~SSORPreconditioner(){
		delete m_L;
		delete m_R;
		delete m_D;
	};

	Vector solve(const Vector & b) const{
		Vector x = lower_triangular_solve(*m_L,b);
		Vector y = (*m_D)*x;
		return (2-m_w)/m_w*upper_triangular_solve(*m_R, y);
	}
protected:

	const SparseMatrix * m_A;
	SparseMatrix * m_L;
	SparseMatrix * m_R;
	SparseMatrix * m_D;
	double m_w;
};




// incomplete cholesky preconditioner
class ICPreconditioner : public Preconditioner
{
public:
	ICPreconditioner(const SparseMatrix & A){
		m_A = &A;
	}

	Vector solve(const Vector & b) const{

	}
protected:

	const SparseMatrix * m_A;
};

#endif