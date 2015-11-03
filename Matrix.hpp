#ifndef _MATRIX_H
#define _MATRIX_H

class Matrix{
public:

	Matrix();

	// create and allocate a new matrix
	Matrix(unsigned int mrows, unsigned int ncols);

	// create a matrix initialized by existing data
	Matrix(unsigned int mrows, unsigned int ncols, const double * data);

	// copy constructor
	Matrix(const Matrix & mtx);

	// assignment operator
	Matrix & operator=(const Matrix& mtx);

	// overloaded assignment for submatrix assignment
	Matrix & operator=(const Matrix_Proxy & mtxp);

	// create a submatrix of an existing matrix  (no new memory allocation)
	Matrix_Proxy & operator()(unsigned int rowStart, unsigned int rowEnd, unsigned int colStart, unsigned int colEnd);

	// accessing an element
	double & operator()(unsigned int i, unsigned int j);

	// Matrix-Matrix multiplication
	Matrix & operator*(const Matrix & A, const Matrix & B);

	// member functions
	void size(unsigned int & rows, unsigned int & cols) const {rows = m_mrows; cols = m_ncols; return;};
	unsigned int rows() const {return m_mrows;};
	unsigned int cols() const {return m_ncols;};
	const double * data() const {return m_data;};
	void transpose(); // see: http://stackoverflow.com/questions/16737298/what-is-the-fastest-way-to-transpose-a-matrix-in-c

protected:

	unsigned int m_mrows, m_ncols;
	double * m_data;

};

ostream& operator<<(ostream& os, const Matrix & mtx){
	os << 
}


class Matrix_Proxy{

	friend class Matrix;

public:

	Matrix_Proxy(Matrix & mtx);

	Matrix_Proxy & operator=(Matrix & mtx);

protected:

private:

	unsigned int m_mrows, m_ncols;
	unsigned int m_rowStart, m_rowEnd, m_colStart, m_colEnd;
	Matrix & m_parent;

};


// vector is derived from matrix. Has either 1 column or 1 row
class Vector : public Matrix{
public:

	Vector();

	Vector(unsigned int length);

	Vector(unsigned int length, const double * data);

protected:

private:

};


#endif