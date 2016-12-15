#ifndef _ABSTRACTMATRIX_H
#define _ABSTRACTMATRIX_H

#include <math.h>
#include <iostream>
#include <stdlib.h>


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

	virtual void size(std::size_t & rows, std::size_t & cols) const = 0;
	virtual std::size_t length() const = 0;
	virtual double norm() const = 0;
	// virtual void transpose() = 0;

protected:

};



// *******************************************************************************


class AbstractMatrix
{
public:

	//virtual Abstract_Matrix & operator*(const Abstract_Matrix & vct) const = 0;
	//virtual Abstract_Vector & operator*(const Abstract_Vector & vct) const = 0;
	virtual void size(std::size_t & rows, std::size_t & cols) const = 0;
	virtual std::size_t rows() const = 0;
	virtual std::size_t cols() const = 0;
	virtual double norm() const = 0;
	// virtual void transpose() = 0;

protected:

};

#endif