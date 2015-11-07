#ifndef _MATRIX_H
#define _MATRIX_H

#include <iostream>
#include <stdlib.h>


// forward declare some classes
class Matrix;
class Vector;
class Matrix_Proxy;
class Vector_Proxy;


class Abstract_Matrix
{
public:

	virtual void size(std::size_t & rows, std::size_t & cols) const = 0;
	virtual std::size_t rows() const = 0;
	virtual std::size_t cols() const = 0;
	virtual double norm() const = 0;
	// virtual Vector_Proxy & row(unsigned int i) = 0;
	// virtual Vector_Proxy & col(unsigned int j) = 0;
	// virtual void transpose() = 0;

protected:

};






class Matrix_Proxy
{
	friend class Matrix;
public:

	Matrix_Proxy(Matrix & mtx, unsigned int rowStart, unsigned int rowEnd, unsigned int colStart, unsigned int colEnd)
	: m_rowStart	(rowStart)
	, m_rowEnd		(rowEnd)
	, m_colStart	(colStart)
	, m_colEnd		(colEnd)
	, m_parent		(mtx)
	, m_mrows 		(rowEnd-rowStart+1)
	, m_ncols		(colEnd-colStart+1)
	{
	}

	// copy constructor
	Matrix_Proxy(const Matrix_Proxy & mtxp)
	: m_rowStart	(mtxp.m_rowStart)
	, m_rowEnd		(mtxp.m_rowEnd)
	, m_colStart	(mtxp.m_colStart)
	, m_colEnd		(mtxp.m_colEnd)
	, m_parent		(mtxp.m_parent)
	, m_mrows 		(mtxp.m_mrows)
	, m_ncols		(mtxp.m_ncols)
	{
	}


	Matrix_Proxy & operator=(const Matrix & mtx);


protected:

private:

	unsigned int m_mrows, m_ncols;
	unsigned int m_rowStart, m_rowEnd, m_colStart, m_colEnd;
	Matrix & m_parent;

};





class Vector_Proxy{

	friend class Vector;
	friend class Matrix;

public:

	// constructor
	Vector_Proxy(double * dataptr, std::size_t length, std::size_t stride, bool is_column=true)
	: m_length		(length)
	, m_stride		(stride)
	, m_dataptr		(dataptr)
	, m_is_column 	(is_column)
	{
	}

	// copy constructor
	Vector_Proxy(const Vector_Proxy & vctp)
	: m_length		(vctp.m_length)
	, m_stride 		(vctp.m_length)
	, m_dataptr 	(vctp.m_dataptr)
	, m_is_column 	(vctp.m_is_column)
	{
	}

	Vector_Proxy & operator=(const Vector & vct);

	Vector_Proxy & operator=(Vector_Proxy & vctp);

protected:

	std::size_t m_length;
	bool m_is_column;

	double * m_dataptr;
	unsigned int m_stride;

};





class Matrix : public Abstract_Matrix
{
public:

	Matrix()
		: m_mrows	(0)
		, m_ncols	(0)
		, m_len 	(0)
		, m_data 	(nullptr)
	{
	}

	// create and allocate a new matrix
	Matrix(unsigned int mrows, unsigned int ncols)
		: m_mrows	(mrows)
		, m_ncols	(ncols)
		, m_len 	(mrows*ncols)
		, m_data 	(m_len ? new double[m_len] : nullptr)
	{
	}

	// create a matrix initialized by existing data
	Matrix(unsigned int mrows, unsigned int ncols, const double * data)
		: m_mrows	(mrows)
		, m_ncols	(ncols)
		, m_len 	(mrows*ncols)
		, m_data 	(m_len ? new double[m_len] : nullptr)
	{
		std::copy(data, data + m_len, m_data);
	}

	// copy constructor
	Matrix(const Matrix & mtx)
		: m_mrows	(mtx.m_mrows)
		, m_ncols	(mtx.m_ncols)
		, m_len 	(mtx.m_len)
		, m_data 	(m_len ? new double[m_len] : nullptr)
	{
		std::copy(mtx.m_data, mtx.m_data + m_len, m_data);
	}

	// constructor from proxy
	Matrix(const Matrix_Proxy & mtxp)
		: m_mrows	(mtxp.m_mrows)
		, m_ncols	(mtxp.m_ncols)
		, m_len 	(mtxp.m_mrows*mtxp.m_ncols)
		, m_data 	(m_len ? new double[m_len] : nullptr)
	{
		Matrix m;
		m = mtxp;
		swap(*this, m);
	}

	// destructor
	~Matrix()
	{
		if (m_data != nullptr) delete[] m_data;
	}

	friend void swap(Matrix & m1, Matrix & m2)
	{
		using std::swap;
		swap(m1.m_mrows, m2.m_mrows);
		swap(m1.m_ncols, m2.m_ncols);
		swap(m1.m_len, m2.m_len);
		swap(m1.m_data, m2.m_data);

	}

	// assignment operator
	Matrix & operator=(Matrix& mtx)
	{
		swap(*this, mtx);
		return *this;
	}

	// overloaded assignment for submatrix assignment
	Matrix & operator=(const Matrix_Proxy & mtxp)
	{
		Matrix m(mtxp.m_mrows, mtxp.m_ncols);
		
		unsigned int ix=0, jx=0;
		for (auto i=mtxp.m_rowStart; i<mtxp.m_rowEnd+1; i++)
		{
			for (auto j=mtxp.m_colStart; j<mtxp.m_colEnd+1; j++)
			{
				m(ix, jx) = mtxp.m_parent(i,j);
				jx++;
			}
			ix++;
			jx=0;
		}

		swap(*this, m);
		return *this;
	}

	// // Matrix-Matrix multiplication
	// Matrix & operator*(const Matrix & A, const Matrix & B);

	// scalar multiplication
	Matrix operator*(double val)
	{
		Matrix out(*this);
		for (auto i=0; i<m_len; i++) out.m_data[i]*=val;
		return out;
	}

	// scalar division
	Matrix operator/(double val)
	{
		Matrix out(*this);
		for (auto i=0; i<m_len; i++) out.m_data[i]/=val;
		return out;
	}

	// scalar addition
	Matrix operator+(double val)
	{
		Matrix out(*this);
		for (auto i=0; i<m_len; i++) out.m_data[i]+=val;
		return out;
	}

	// scalar subtraction
	Matrix operator-(double val)
	{
		Matrix out(*this);
		for (auto i=0; i<m_len; i++) out.m_data[i]-=val;
		return out;
	}

	// create a submatrix of an existing matrix  (no new memory allocation)
	Matrix_Proxy operator()(unsigned int rowStart, unsigned int rowEnd, unsigned int colStart, unsigned int colEnd)
	{
		return Matrix_Proxy(*this, rowStart, rowEnd, colStart, colEnd);
	}

	// accessing an element
	double & operator()(unsigned int i, unsigned int j)
	{
		return m_data[m_mrows*j + i];
	}

	// accessing an element
	double operator()(unsigned int i, unsigned int j) const
	{
		return m_data[m_mrows*j + i];
	}


	void fill(double fillval)
	{
		for (auto i=0; i<m_len; i++) m_data[i] = fillval;
	}


	// print info
	void print_summary() const 
	{
		std::cout << "rows: " << m_mrows << std::endl;
		std::cout << "cols: " << m_ncols << std::endl;
		std::cout << "total: " << m_len << std::endl;
		std::cout << " " << std::endl;
	}

	
	// member functions
	void size(std::size_t & rows, std::size_t & cols) const {rows = m_mrows; cols = m_ncols; return;};
	std::size_t rows() const {return m_mrows;};
	std::size_t cols() const {return m_ncols;};
	double norm() const {return 0;};

	const double * data() const {return m_data;};

	// void transpose(); // see: http://stackoverflow.com/questions/16737298/what-is-the-fastest-way-to-transpose-a-matrix-in-c

	// Vector_Proxy & row(unsigned int i);

	// Vector_Proxy & col(unsigned int j);

	friend std::ostream & operator<<(std::ostream & os, const Matrix & mtx);

protected:

	std::size_t m_mrows, m_ncols;
	std::size_t m_len;
	double * m_data;

};

std::ostream& operator<<(std::ostream & os, const Matrix & mtx)
{
	os << std::endl;
	for (auto i=0; i<mtx.rows(); i++)
	{
		for (auto j=0; j<mtx.cols(); j++)
		{
			os << mtx(i,j) << "\t" ;
		}
		os << std::endl;
	}
	return os;
}





// // vector is derived from matrix. Has either 1 column or 1 row
// class Vector : public Abstract_Matrix{
// public:

// 	Vector();

// 	Vector(unsigned int length);

// 	Vector(unsigned int length, const double * data);

// protected:

// private:

// };


#endif