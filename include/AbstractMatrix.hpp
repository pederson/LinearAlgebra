#ifndef _ABSTRACTMATRIX_H
#define _ABSTRACTMATRIX_H

#include <math.h>
#include <iostream>
#include <stdlib.h>

// namespace libra{

// forward declare some classes
class Matrix;
class Vector;
class Matrix_Proxy;
class Vector_Proxy;
class SparseMatrix;
class SparseVector;



// *******************************************************************************


class AbstractVector
{
public:

	//virtual static double dot(const AbstractVector & av1, const AbstractVector & av2) const = 0;
	virtual void size(std::size_t & rows, std::size_t & cols) const = 0;
	virtual std::size_t length() const = 0;
	virtual double norm() const = 0;
	virtual void transpose() = 0;

protected:

};



// *******************************************************************************


class AbstractMatrix
{
public:

	//virtual AbstractMatrix & operator*(const AbstractMatrix & vct) const = 0;
	virtual Vector operator*(const Vector & vct) const = 0;
	virtual Vector Tmult(const Vector & vct) const = 0;
	//virtual SparseVector operator*(const SparseVector & vct) const = 0;
	virtual void size(std::size_t & rows, std::size_t & cols) const = 0;
	virtual std::size_t rows() const = 0;
	virtual std::size_t cols() const = 0;
	//virtual double norm() const = 0;
	virtual void transpose() = 0;

protected:

};



// ******************************************************************************

// useful functions
template <typename T> int sgn(T val) {
    return (T(0) < val) - (val < T(0));
}

// }
#endif