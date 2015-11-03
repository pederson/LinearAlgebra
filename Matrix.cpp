#include "Matrix.hpp"


using namespace std;

Matrix::Matrix(){
	m_data = nullptr;
	m_mrows = 0;
	m_ncols = 0;
}

// create and allocate a new matrix
Matrix::Matrix(unsigned int mrows, unsigned int ncols){
	m_mrows = mrows;
	m_ncols = ncols;
	m_data = new double[mrows*ncols];
}

// create a matrix initialized by existing data
Matrix::Matrix(unsigned int mrows, unsigned int ncols, const double * data{
	m_mrows = mrows;
	m_ncols = ncols;
	m_data = new double[mrows*ncols];
	//for (auto i=0; i<mrows*ncols; i++) m_data[i] = data[i];
	memcpy(m_data, data, m_mrows*m_ncols*sizeof(double));
}

// copy constructor
Matrix::Matrix(const Matrix & mtx){
	m_data = new double[mtx.m_mrows*mtx.m_ncols];
	memcpy(m_data, mtx.m_data, m_mrows*m_ncols*sizeof(double));
}

// assignment operator
Matrix & Matrix::operator=(const Matrix& mtx){

	return *this;
}

// overloaded assignment for submatrix assignment
Matrix & Matrix::operator=(const Matrix_Proxy & mtxp){
	m_data = new double[mtxp.m_mrows*mtxp.m_ncols];
	unsigned int cnt = 0;
	for (auto i=mtxp.m_rowStart; i<=mtxp.m_rowEnd; i++){
		for (auto j=mtxp.m_colStart; j<=mtxp.m_colEnd; j++){
			m_data[cnt++] = mtxp.m_parent(i, j);
		}
	}

	return *this;
}

// create a submatrix of an existing matrix  (no new memory allocation)
Matrix_Proxy & Matrix::operator()(unsigned int rowStart, unsigned int rowEnd, unsigned int colStart, unsigned int colEnd){
	Matrix_Proxy mprx(*this);
	mprx.m_rowStart = rowStart;
	mprx.m_rowEnd = rowEnd;
	mprx.m_colStart = colStart;
	mprx.m_colEnd = colEnd;
}

// accessing an element
double & Matrix::operator()(unsigned int i, unsigned int j){
	return data[m_mrows*j + i];
}

// Matrix-Matrix multiplication
Matrix & Matrix::operator*(const Matrix & A, const Matrix & B){
	if (A.cols() != B.rows()){
		throw "Matrix dimensions must match!";
	}

	Matrix result(A.rows(), B.cols());
	for (auto i=0; i<result.m_mrows; i++){
		for (auto j=0; j<result.m_ncols; j++){
			result(i, j) = 
		}
	}
}


void Matrix::transpose(); // see: http://stackoverflow.com/questions/16737298/what-is-the-fastest-way-to-transpose-a-matrix-in-c



Matrix_Proxy::Matrix_Proxy(Matrix & mtx){
	m_parent = mtx;
}

Matrix_Proxy & Matrix_Proxy::operator=(Matrix & mtx){

}

// Vector::Vector();

// Vector::Vector(unsigned int length);

// Vector::Vector(unsigned int length, const double * data);


int main(unsigned int argc, char * argv[]){

	Matrix A(10, 10);

}