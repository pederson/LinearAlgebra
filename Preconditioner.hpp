#ifndef _PRECONDITIONER_H
#define _PRECONDITIONER_H

#include "Matrix.hpp"
#include "SparseMatrix.hpp"

Vector diagonal_solve(const SparseMatrix & A, const Vector & b);
Vector lower_triangular_solve(const SparseMatrix & A, const Vector & b);
Vector upper_triangular_solve(const SparseMatrix & A, const Vector & b);
Vector diag(const SparseMatrix & A);
SparseMatrix spdiag(const Vector & x, int j);
SparseMatrix strictly_lower(const SparseMatrix & A);
SparseMatrix strictly_upper(const SparseMatrix & A);
void icholesky(const SparseMatrix & A, SparseMatrix & R);
void ilu(const SparseMatrix & A, SparseMatrix & L, SparseMatrix & U);
void amg_setup(const SparseMatrix & A, std::vector<SparseMatrix *> & Ws, std::vector<SparseMatrix *> & As);
void amgv(const SparseMatrix & A, const Vector & b, Vector & x, unsigned int level,
		  const std::vector<SparseMatrix *> & Ws, const std::vector<SparseMatrix *> & As, 
		  unsigned int v1, unsigned int v2);

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
		m_D = spdiag(diag(*m_A));
	}

	Vector solve(const Vector & b) const{
		Vector x = lower_triangular_solve(*m_A,b);
		Vector y = m_D*x;
		return upper_triangular_solve(*m_A, y);
	}
protected:

	const SparseMatrix * m_A;
	SparseMatrix m_D;
};


// successive overrelaxation preconditioner
// P = (D/w+L) = (D+L) + (1-w)/w *D
class SORPreconditioner : public Preconditioner
{
public:
	SORPreconditioner(const SparseMatrix & A, double w){
		m_A = &A;
		m_w = w;
		m_P = spdiag(diag(*m_A))/w + strictly_lower(*m_A);
	}

	Vector solve(const Vector & b) const{
		return lower_triangular_solve(m_P, b);
	}
protected:

	const SparseMatrix * m_A;
	SparseMatrix m_P;
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
		m_L = spdiag(diag(*m_A))/w + strictly_lower(*m_A);
		m_R = spdiag(diag(*m_A))/w + strictly_upper(*m_A);
		m_D = spdiag(diag(*m_A));
	}

	Vector solve(const Vector & b) const{
		Vector x = lower_triangular_solve(m_L,b);
		Vector y = m_D*x;
		return (2-m_w)/m_w*upper_triangular_solve(m_R, y);
	}
protected:

	const SparseMatrix * m_A;
	SparseMatrix m_L;
	SparseMatrix m_R;
	SparseMatrix m_D;
	double m_w;
};




// incomplete cholesky preconditioner
class ICPreconditioner : public Preconditioner
{
public:
	ICPreconditioner(const SparseMatrix & A){
		m_A = &A;
		icholesky(A,m_L);
		m_R = ~m_L;
	}

	Vector solve(const Vector & b) const{
		Vector x = lower_triangular_solve(m_L,b);
		return upper_triangular_solve(m_R,x);
	}
protected:

	const SparseMatrix * m_A;
	SparseMatrix m_R, m_L;
};



// incomplete LU preconditioner
class ILUPreconditioner : public Preconditioner
{
public:
	ILUPreconditioner(const SparseMatrix & A){
		m_A = &A;
		ilu(A,m_L,m_U);
	}

	Vector solve(const Vector & b) const{
		Vector x = lower_triangular_solve(m_L,b);
		return upper_triangular_solve(m_U,x);
	}
protected:

	const SparseMatrix * m_A;
	SparseMatrix m_U, m_L;
};


// Algebraic multigrid preconditioner
class AMGPreconditioner : public Preconditioner
{
public:
	AMGPreconditioner(const SparseMatrix & A){
		m_A = &A;
		amg_setup(A, Ws, As);
	}

	~AMGPreconditioner(){
		for (auto i=0; i<Ws.size(); i++) delete Ws[i];
		for (auto i=0; i<As.size(); i++) delete As[i];
	}

	Vector solve(const Vector & b) const{
		Vector x(m_A->cols()); x.fill(0);
		amgv(*m_A, b, x, 0, Ws, As, 2, 2);
		return x;
	}

private:

	const SparseMatrix * m_A;
	std::vector<SparseMatrix *> Ws, As;
};

#endif