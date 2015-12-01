#ifndef _MATRIX_H
#define _MATRIX_H

#include <math.h>
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

	//virtual Abstract_Matrix & operator*(const Abstract_Matrix & vct) const = 0;
	virtual void size(std::size_t & rows, std::size_t & cols) const = 0;
	virtual std::size_t rows() const = 0;
	virtual std::size_t cols() const = 0;
	virtual double norm() const = 0;
	virtual Vector_Proxy row(unsigned int i) = 0;
	virtual const Vector_Proxy row(unsigned int i) const = 0;
	virtual Vector_Proxy col(unsigned int j) = 0;
	virtual const Vector_Proxy col(unsigned int j) const = 0;

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

	// Matrix_Proxy(const Matrix & mtx, unsigned int rowStart, unsigned int rowEnd, unsigned int colStart, unsigned int colEnd)
	// 	: m_rowStart	(rowStart)
	// 	, m_rowEnd		(rowEnd)
	// 	, m_colStart	(colStart)
	// 	, m_colEnd		(colEnd)
	// 	, m_parent		(mtx)
	// 	, m_mrows 		(rowEnd-rowStart+1)
	// 	, m_ncols		(colEnd-colStart+1)
	// {
	// }

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

	void operator+=(const Matrix & mtx);

	void operator-=(const Matrix & mtx);

	Vector operator*(const Vector & vct);


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

	static double dot(const Vector_Proxy & vp1, const Vector_Proxy & vp2)
	{
		double result = 0.0;
		if (vp1.m_length == vp2.m_length)
		{
			for (auto i=0; i< vp1.m_length; i++)
			{
				result += *(vp1.m_dataptr + i*vp1.m_stride) * *(vp2.m_dataptr + i*vp2.m_stride);
			}
		}
		else
		{
			throw "Vectors must be of same size!";
		}

		return result;
	}

	Vector_Proxy & operator=(const Vector & vct);

	Vector operator*(double val);

	Vector operator/(double val);

	Vector operator+(double val);

	Vector operator-(double val);

	void operator*=(double val);

	void operator/=(double val);

	void operator+=(double val);

	void operator-=(double val);

	double & operator()(unsigned int i)
	{
		return m_dataptr[i*m_stride];
	}

	double operator()(unsigned int i) const
	{
		return m_dataptr[i*m_stride];
	}

	void operator=(Vector_Proxy & vctp)
	{
		for (auto i=0; i<m_length; i++)
		{
			m_dataptr[i*m_stride] = vctp.m_dataptr[i*vctp.m_stride];
		}
	}

	std::size_t length() const {return m_length;}; 

	friend std::ostream & operator<<(std::ostream & os, const Vector_Proxy & vctp);

protected:

	std::size_t m_length;
	bool m_is_column;

	double * m_dataptr;
	unsigned int m_stride;

};

std::ostream& operator<<(std::ostream & os, const Vector_Proxy & vctp)
{
	os << std::endl;
	os << std::scientific;
	if (vctp.m_is_column)
	{
		for (auto i=0; i<vctp.m_length; i++)
		{
			os << *(vctp.m_dataptr + i*vctp.m_stride) << std::endl;
		}
	}
	else
	{
		for (auto i=0; i<vctp.m_length; i++)
		{
			os << *(vctp.m_dataptr + i*vctp.m_stride) << "\t" ;
		}
		os << std::endl;
	}
	

	return os;
}



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

	// Matrix transpose
	Matrix operator~() const
	{
		Matrix m(*this);
		m.transpose();
		return m;
	}

	// Matrix-Matrix multiplication
	Matrix operator*(const Matrix & A) const
	{
		if (m_ncols != A.m_mrows)
		{
			throw "Matrix dimensions do not match!";
		}

		Matrix out(m_mrows, A.m_ncols);
		for (auto i=0; i<m_mrows; i++)
		{
			for (auto j=0; j<A.m_ncols; j++)
			{
				out(i,j) = Vector_Proxy::dot(row(i), A.col(j));
			}
		}
		return out;
	}

	// Matrix-Matrix multiplication
	Matrix & operator*=(const Matrix & A)
	{
		if (m_ncols != A.m_mrows)
		{
			throw "Matrix dimensions do not match!";
		}

		Matrix out(m_mrows, A.m_ncols);
		for (auto i=0; i<m_mrows; i++)
		{
			for (auto j=0; j<A.m_ncols; j++)
			{
				out(i,j) = Vector_Proxy::dot(row(i), A.col(j));
			}
		}
		swap(out, *this);
		return *this;
	}

	// Matrix-Vector multiplication
	Vector operator*(const Vector & v) const;

	// Matrix-Matrix addition
	Matrix operator+(const Matrix & mtx) const
	{
		if (mtx.m_mrows != m_mrows || mtx.m_ncols != m_ncols)
		{
			throw "Matrix dimensions do not match!";
		}

		Matrix out(m_mrows, m_ncols);
		for (auto i=0; i<m_mrows; i++)
		{
			for (auto j=0; j<m_ncols; j++)
			{
				out(i,j) = (*this)(i,j) + mtx(i,j);
			}
		}
		return out;
	}

	// Matrix-Matrix subtraction
	Matrix operator-(const Matrix & mtx) const
	{
		if (mtx.m_mrows != m_mrows || mtx.m_ncols != m_ncols)
		{
			throw "Matrix dimensions do not match!";
		}

		Matrix out(m_mrows, m_ncols);
		for (auto i=0; i<m_mrows; i++)
		{
			for (auto j=0; j<m_ncols; j++)
			{
				out(i,j) = (*this)(i,j) - mtx(i,j);
			}
		}
		return out;
	}

	// Matrix-Matrix addition shorthand
	Matrix & operator+=(const Matrix & mtx)
	{
		if (mtx.m_mrows != m_mrows || mtx.m_ncols != m_ncols)
		{
			throw "Matrix dimensions do not match!";
		}

		for (auto i=0; i<m_mrows; i++)
		{
			for (auto j=0; j<m_ncols; j++)
			{
				(*this)(i,j) += mtx(i,j);
			}
		}
		return *this;
	}

	// Matrix-Matrix subtraction shorthand
	Matrix & operator-=(const Matrix & mtx)
	{
		if (mtx.m_mrows != m_mrows || mtx.m_ncols != m_ncols)
		{
			throw "Matrix dimensions do not match!";
		}

		for (auto i=0; i<m_mrows; i++)
		{
			for (auto j=0; j<m_ncols; j++)
			{
				(*this)(i,j) -= mtx(i,j);
			}
		}
		return *this;
	}

	// scalar multiplication
	Matrix operator*(double val) const
	{
		Matrix out(*this);
		for (auto i=0; i<m_len; i++) out.m_data[i]*=val;
		return out;
	}

	// scalar division
	Matrix operator/(double val) const
	{
		Matrix out(*this);
		for (auto i=0; i<m_len; i++) out.m_data[i]/=val;
		return out;
	}

	// scalar addition
	Matrix operator+(double val) const
	{
		Matrix out(*this);
		for (auto i=0; i<m_len; i++) out.m_data[i]+=val;
		return out;
	}

	// scalar subtraction
	Matrix operator-(double val) const
	{
		Matrix out(*this);
		for (auto i=0; i<m_len; i++) out.m_data[i]-=val;
		return out;
	}

	// scalar multiplication
	Matrix & operator*=(double val)
	{
		for (auto i=0; i<m_len; i++) m_data[i]*=val;
		return *this;
	}

	// scalar division
	Matrix & operator/=(double val)
	{
		for (auto i=0; i<m_len; i++) m_data[i]/=val;
		return *this;
	}

	// scalar addition
	Matrix & operator+=(double val)
	{
		for (auto i=0; i<m_len; i++) m_data[i]+=val;
		return *this;
	}

	// scalar subtraction
	Matrix & operator-=(double val)
	{
		for (auto i=0; i<m_len; i++) m_data[i]-=val;
		return *this;
	}

	// create a submatrix of an existing matrix  (no new memory allocation)
	Matrix_Proxy operator()(unsigned int rowStart, unsigned int rowEnd, unsigned int colStart, unsigned int colEnd)
	{
		return Matrix_Proxy(*this, rowStart, rowEnd, colStart, colEnd);
	}

	// // create a const submatrix of an existing matrix (no new memory allocation)
	// const Matrix_Proxy operator()(unsigned int rowStart, unsigned int rowEnd, unsigned int colStart, unsigned int colEnd) const
	// {
	// 	return Matrix_Proxy(*this, rowStart, rowEnd, colStart, colEnd);
	// }

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

	// row access
	Vector_Proxy row(unsigned int i)
	{
		return Vector_Proxy(&m_data[i], m_ncols, m_mrows, false);
	}

	// row access const
	const Vector_Proxy row(unsigned int i) const
	{
		return Vector_Proxy(&m_data[i], m_ncols, m_mrows, false);
	}

	// subrow access
	Vector_Proxy subrow(unsigned int row, unsigned int colStart, unsigned int colEnd)
	{
		return Vector_Proxy(&m_data[m_mrows*colStart + row], colEnd-colStart+1, m_mrows, false);
	}

	// subrow access const
	const Vector_Proxy subrow(unsigned int row, unsigned int colStart, unsigned int colEnd) const
	{
		return Vector_Proxy(&m_data[m_mrows*colStart + row], colEnd-colStart+1, m_mrows, false);
	}

	// column access
	Vector_Proxy col(unsigned int j)
	{
		return Vector_Proxy(&m_data[j*m_mrows], m_mrows, 1, true);
	}

	// column access const
	const Vector_Proxy col(unsigned int j) const
	{
		return Vector_Proxy(&m_data[j*m_mrows], m_mrows, 1, true);
	}

	// subcolumn access
	Vector_Proxy subcol(unsigned int col, unsigned int rowStart, unsigned int rowEnd)
	{
		return Vector_Proxy(&m_data[col*m_mrows + rowStart], rowEnd-rowStart+1, 1, true);
	}

	// subcolumn access const
	const Vector_Proxy subcol(unsigned int col, unsigned int rowStart, unsigned int rowEnd) const
	{
		return Vector_Proxy(&m_data[col*m_mrows + rowStart], rowEnd-rowStart+1, 1, true);
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

	void transpose();

	friend std::ostream & operator<<(std::ostream & os, const Matrix & mtx);

protected:

	std::size_t m_mrows, m_ncols;
	std::size_t m_len;
	double * m_data;

};

std::ostream& operator<<(std::ostream & os, const Matrix & mtx)
{
	os << std::endl;
	os << std::scientific;
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



// vector is derived from matrix. Has either 1 column or 1 row
class Vector : public Matrix{
public:

	Vector()
		: Matrix(0, 0)
		, m_is_column (true)
	{
	}

	// create and allocate a new matrix
	Vector(unsigned int length)
		: Matrix(length, 1)
		, m_is_column (true)
	{
	}

	// create a matrix initialized by existing data
	Vector(unsigned int length, const double * data)
		: Matrix(length, 1, data)
		, m_is_column (true)
	{
		//std::copy(data, data + m_len, m_data);
	}


	// copy constructor
	Vector(const Vector & vct)
		: Matrix(vct.m_mrows, vct.m_ncols, vct.m_data)
		, m_is_column (vct.m_is_column)
	{
		//std::copy(vct.m_data, vct.m_data + m_len, m_data);
	}

	// construct via a matrix with 1 row or column
	Vector(const Matrix & mtx)
		: Matrix(mtx)
		, m_is_column	(true)
	{
		if (mtx.rows() != 1 && mtx.cols() != 1)
		{
			std::cout << "Cannot convert matrix to vector!" << std::endl;
			throw -1;
		}
	}

	// constructor from proxy
	Vector(const Vector_Proxy & vctp)
		: Matrix((vctp.m_is_column? vctp.m_length : 1), (vctp.m_is_column? 1 : vctp.m_length))
		, m_is_column (vctp.m_is_column)
	{
		for (auto i=0; i<vctp.m_length; i++)
		{
			m_data[i] = *(vctp.m_dataptr + i*vctp.m_stride);
		}
	}

	// destructor
	~Vector()
	{
	}

	friend void swap(Vector & v1, Vector & v2)
	{
		using std::swap;
		swap(v1.m_mrows, v2.m_mrows);
		swap(v1.m_ncols, v2.m_ncols);
		swap(v1.m_len, v2.m_len);
		swap(v1.m_data, v2.m_data);
		swap(v1.m_is_column, v2.m_is_column);

	}

	// assignment operator
	Vector & operator=(Vector& vct)
	{
		// implicit copy and swap
		swap(*this, vct);
		return *this;
	}

	// overloaded assignment for submatrix assignment
	Vector & operator=(const Vector_Proxy & vctp)
	{

		// copy and swap
		Vector m(vctp);
		swap(*this, m);
		return *this;
	}

	// overloaded assignment converting matrix
	Vector & operator=(Matrix mtx)
	{
		// copy and swap
		Vector m(mtx);
		swap(*this, m);
		return *this;
	}

	Vector operator~()
	{
		Vector out(*this);
		out.transpose();
		return out;
	}

	// Vector-Matrix multiplication
	// this returns a matrix to enable outer products
	Matrix operator*(const Matrix & A) const
	{

		if (m_ncols != A.rows())
		{
			std::cout << "Matrix dimensions do not match!" << std::endl;
			throw -1;
		}

		Matrix out(m_mrows, A.cols());
		//Vector_Proxy mine(m_data, m_len, 1, m_is_column);
		for (auto i=0; i<m_mrows; i++)
		{
			for (auto j=0; j<A.cols(); j++)
			{
				Vector_Proxy myrow(&m_data[i], m_ncols, 1, false);
				out(i,j) = Vector_Proxy::dot(myrow, A.col(j));
			}
		}

		return out;
	}

	// Vector operator*(const Vector & vct) const
	// {
	// 	if (m_ncols != vct.rows())
	// 	{
	// 		std::cout << "Matrix dimensions do not match!" << std::endl;
	// 		throw -1;
	// 	}

	// 	Vector out(m_mrows);
	// 	for (auto i=0; i<m_mrows; i++)
	// 	{
	// 		Vector_Proxy myrow(&m_data[i], m_ncols, 1, false);
	// 		out(i) = Vector_Proxy::dot(myrow, vct.col(0));
	// 	}
	// 	return out;

	// }

	// Vector-Vector addition
	Vector operator+(const Vector & vct) const
	{
		if (m_len != vct.m_len)
		{
			std::cout << "Vector dimensions do not match!" << std::endl;
			throw -1;
		}

		Vector out(m_len);
		for (auto i=0; i<m_len; i++)
		{
			out(i) = m_data[i] + vct.m_data[i];
		}

		return out;
	}

	// Vector-Vector subtraction
	Vector operator-(const Vector & vct) const
	{
		if (m_len != vct.m_len)
		{
			std::cout << "Vector dimensions do not match!" << std::endl;
			throw -1;
		}

		Vector out(m_len);
		for (auto i=0; i<m_len; i++)
		{
			out(i) = m_data[i] - vct.m_data[i];
		}

		return out;
	}

	// Vector-Vector addition shorthand
	Vector & operator+=(const Vector & vct)
	{
		if (m_len != vct.m_len)
		{
			std::cout << "Vector dimensions do not match!" << std::endl;
			throw -1;
		}

		for (auto i=0; i<m_len; i++)
		{
			m_data[i] += vct.m_data[i];
		}

		return *this;
	}

	// Vector-Vector subtraction shorthand
	Vector & operator-=(const Vector & vct)
	{
		if (m_len != vct.m_len)
		{
			std::cout << "Vector dimensions do not match!" << std::endl;
			throw -1;
		}

		for (auto i=0; i<m_len; i++)
		{
			m_data[i] -= vct.m_data[i];
		}

		return *this;
	}

	// scalar multiplication
	Vector operator*(double val) const
	{
		Vector out(*this);
		for (auto i=0; i<m_len; i++) out.m_data[i]*=val;
		return out;
	}

	// scalar division
	Vector operator/(double val) const
	{
		Vector out(*this);
		for (auto i=0; i<m_len; i++) out.m_data[i]/=val;
		return out;
	}

	// scalar addition
	Vector operator+(double val) const
	{
		Vector out(*this);
		for (auto i=0; i<m_len; i++) out.m_data[i]+=val;
		return out;
	}

	// scalar subtraction
	Vector operator-(double val) const
	{
		Vector out(*this);
		for (auto i=0; i<m_len; i++) out.m_data[i]-=val;
		return out;
	}

	// scalar multiplication
	Vector & operator*=(double val)
	{
		for (auto i=0; i<m_len; i++) m_data[i]*=val;
		return *this;
	}

	// scalar division
	Vector & operator/=(double val)
	{
		for (auto i=0; i<m_len; i++) m_data[i]/=val;
		return *this;
	}

	// scalar addition
	Vector & operator+=(double val)
	{
		for (auto i=0; i<m_len; i++) m_data[i]+=val;
		return *this;
	}

	// scalar subtraction
	Vector & operator-=(double val)
	{
		for (auto i=0; i<m_len; i++) m_data[i]-=val;
		return *this;
	}

	double & operator()(unsigned int i)
	{
		return m_data[i];
	}

	double operator()(unsigned int i) const
	{
		return m_data[i];
	}

	// create a subvector of an existing vector  (no new memory allocation)
	Vector_Proxy operator()(unsigned int indStart, unsigned int indEnd)
	{
		return Vector_Proxy(m_data+indStart, indEnd-indStart+1, 1, m_is_column);
	}

	static double dot(const Vector & v1, const Vector & v2)
	{
		double result = 0.0;
		if (v1.m_len == v2.m_len)
		{
			for (auto i=0; i< v1.m_len; i++)
			{
				result += v1.m_data[i] * v2.m_data[i];
			}
		}
		else
		{
			throw "Vectors must be of same size!";
		}

		return result;
	}

	double length() const {return m_len;};

	double norm() const 
	{
		double nm=0; 
		for (auto i=0; i<m_len; i++)
		{
			nm+=m_data[i]*m_data[i];
		}
		nm = sqrt(nm);
		return nm;
	}

	void transpose() {m_is_column = !m_is_column; std::swap(m_mrows, m_ncols);};

protected:

	bool m_is_column;

};









// global operators - Matrix
Matrix operator*(double val, const Matrix & mtx)
{
	Matrix out(mtx);
	for (auto i=0; i<mtx.rows(); i++)
	{
		for (auto j=0; j<mtx.cols(); j++)
		{ 
			out(i,j)*=val;
		}
	}
	return out;
}

Matrix operator+(double val, const Matrix & mtx)
{
	Matrix out(mtx);
	for (auto i=0; i<mtx.rows(); i++)
	{
		for (auto j=0; j<mtx.cols(); j++)
		{ 
			out(i,j)+=val;
		}
	}
	return out;
}

Matrix operator-(double val, const Matrix & mtx)
{
	Matrix out(mtx.rows(),mtx.cols());
	out.fill(val);
	for (auto i=0; i<mtx.rows(); i++)
	{
		for (auto j=0; j<mtx.cols(); j++)
		{
			out(i,j) -= mtx(i,j);
		}
	} 
	return out;
}






// global operators - Vector
Vector operator*(double val, const Vector & vct)
{
	Vector out(vct);
	for (auto i=0; i<vct.length(); i++) out(i)*=val;
	return out;
}

Vector operator+(double val, const Vector & vct)
{
	Vector out(vct);
	for (auto i=0; i<vct.length(); i++) out(i)+=val;
	return out;
}

Vector operator-(double val, const Vector & vct)
{
	Vector out(vct.length());
	out.fill(val);
	for (auto i=0; i<vct.length(); i++) out(i) -= vct(i);
	return out;
}







// global operators - Vector_Proxy
Vector operator*(double val, const Vector_Proxy & vct)
{
	Vector out(vct);
	for (auto i=0; i<vct.length(); i++) out(i)*=val;
	return out;
}

Vector operator+(double val, const Vector_Proxy & vct)
{
	Vector out(vct);
	for (auto i=0; i<vct.length(); i++) out(i)+=val;
	return out;
}

Vector operator-(double val, const Vector_Proxy & vct)
{
	Vector out(vct.length());
	out.fill(val);
	for (auto i=0; i<vct.length(); i++) out(i) -= vct(i);
	return out;
}


// things that can't be defined until everything has been defined


// for a better method,
// see: http://stackoverflow.com/questions/16737298/what-is-the-fastest-way-to-transpose-a-matrix-in-c
void Matrix::transpose()
{
	Matrix m(m_ncols, m_mrows);
	Vector c;
	for (auto i=0; i<m_ncols; i++)
	{
		c = col(i);
		c.transpose();

		m.row(i) = c;
	}

	swap(*this, m);
}



Matrix_Proxy & Matrix_Proxy::operator=(const Matrix & mtx)
{
	unsigned int ix=0, jx=0;
	for (auto i=m_rowStart; i<m_rowEnd+1; i++)
	{
		for (auto j=m_colStart; j<m_colEnd+1; j++)
		{
			m_parent(i,j) = mtx(ix, jx);
			jx++;
		}
		ix++;
		jx=0;
	}

	return *this;
}

void Matrix_Proxy::operator+=(const Matrix & mtx)
{
	unsigned int ix=0, jx=0;
	for (auto i=m_rowStart; i<m_rowEnd+1; i++)
	{
		for (auto j=m_colStart; j<m_colEnd+1; j++)
		{
			m_parent(i,j) += mtx(ix, jx);
			jx++;
		}
		ix++;
		jx=0;
	}
	return;
}

void Matrix_Proxy::operator-=(const Matrix & mtx)
{
	unsigned int ix=0, jx=0;
	for (auto i=m_rowStart; i<m_rowEnd+1; i++)
	{
		for (auto j=m_colStart; j<m_colEnd+1; j++)
		{
			m_parent(i,j) -= mtx(ix, jx);
			jx++;
		}
		ix++;
		jx=0;
	}
	return;
}


Vector Matrix_Proxy::operator*(const Vector & vct)
{
	unsigned int m = m_parent.rows();
	Vector out(m_rowEnd-m_rowStart+1);
	out.fill(0);
	for (auto i=0; i<m_colEnd-m_colStart+1; i++)
	{
		out += m_parent.subcol(m_colStart + i, m_rowStart, m_rowEnd)*vct(i);
	}
	return out;
}








Vector_Proxy & Vector_Proxy::operator=(const Vector & vct)
{
	for (auto i=0; i<m_length; i++)
	{
		*(m_dataptr + i*m_stride) = vct(i);
	}

	return *this;
}

Vector Vector_Proxy::operator*(double val)
{
	Vector out(*this);
	for (auto i=0; i<m_length; i++) out(i)*=val;
	return out;
}

Vector Vector_Proxy::operator/(double val)
{
	Vector out(*this);
	for (auto i=0; i<m_length; i++) out(i)/=val;
	return out;
}

Vector Vector_Proxy::operator+(double val)
{
	Vector out(*this);
	for (auto i=0; i<m_length; i++) out(i)+=val;
	return out;
}

Vector Vector_Proxy::operator-(double val)
{
	Vector out(*this);
	for (auto i=0; i<m_length; i++) out(i)-=val;
	return out;
}


void Vector_Proxy::operator*=(double val)
{
	for (auto i=0; i<m_length; i++) m_dataptr[i*m_stride] *= val;
}

void Vector_Proxy::operator/=(double val)
{
	for (auto i=0; i<m_length; i++) m_dataptr[i*m_stride] /= val;
}

void Vector_Proxy::operator+=(double val)
{
	for (auto i=0; i<m_length; i++) m_dataptr[i*m_stride] += val;
}

void Vector_Proxy::operator-=(double val)
{
	for (auto i=0; i<m_length; i++) m_dataptr[i*m_stride] -= val;
}



Vector Matrix::operator*(const Vector & v) const
{
	if (m_ncols != v.rows())
	{
		throw "Matrix dimensions do not match!";
	}

	Vector out(m_mrows);

	out.fill(0);
	for (auto i=0; i<m_ncols; i++)
	{
		const Vector c = col(i);
		// out += col(i)*v(i);
		out += c*v(i);
	}

	return out;
}



#endif