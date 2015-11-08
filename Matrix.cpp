#include "Matrix.hpp"


using namespace std;


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


// // create a submatrix of an existing matrix  (no new memory allocation)
// Matrix_Proxy & Matrix::operator()(unsigned int rowStart, unsigned int rowEnd, unsigned int colStart, unsigned int colEnd){
// 	Matrix_Proxy mprx(*this);
// 	mprx.m_rowStart = rowStart;
// 	mprx.m_rowEnd = rowEnd;
// 	mprx.m_colStart = colStart;
// 	mprx.m_colEnd = colEnd;
// }

// // accessing an element
// double & Matrix::operator()(unsigned int i, unsigned int j){
// 	return data[m_mrows*j + i];
// }

// // Matrix-Matrix multiplication
// Matrix & Matrix::operator*(const Matrix & A, const Matrix & B){
// 	if (A.cols() != B.rows()){
// 		throw "Matrix dimensions must match!";
// 	}

// 	Matrix result(A.rows(), B.cols());
// 	for (auto i=0; i<result.m_mrows; i++){
// 		for (auto j=0; j<result.m_ncols; j++){
// 			result(i, j) = 
// 		}
// 	}
// }


// void Matrix::transpose(); // see: http://stackoverflow.com/questions/16737298/what-is-the-fastest-way-to-transpose-a-matrix-in-c



// Matrix_Proxy::Matrix_Proxy(Matrix & mtx){
// 	m_parent = mtx;
// }

// Matrix_Proxy & Matrix_Proxy::operator=(Matrix & mtx){

// }

// // Vector::Vector();

// // Vector::Vector(unsigned int length);

// // Vector::Vector(unsigned int length, const double * data);



#include <typeinfo>

int main(int argc, char * argv[]){

	unsigned int nrows=20, ncols = 5;

	Matrix A(nrows, ncols);
	for (auto i=0; i<nrows; i++)
	{
		for (auto j=0; j<ncols; j++)
		{
			A(i, j) = i*j;
		}
	}

	// assignment operator
	Matrix B = A;
	cout << "B = " << B << endl;

	A.print_summary();

	// submatrix extraction
	Matrix S = A(4, 8, 3, 4);
	cout << "S = " << S << endl;

	// submatrix assignment
	Matrix C(5, 2);
	C.fill(-10);
	B(0, 4, 2, 3) = C;
	B(15, 19, 0, 1) = S;
	cout << "B = " << B << endl;

	// submatrix assignment via matrix proxy
	// B(0, 1, 0, 1) = S(3, 4, 0, 1);
	// cout << "B = " << B << endl;

	// scalar multiplication, subtraction, addition, division
	Matrix D = (((C*0.5)/2.0) + 10) - 2.5;
	cout << "D = " << D << endl;
	cout << "C = " << C << endl;

	// vector proxy
	cout << "B(:, 3) = " << B.col(3) << endl;
	cout << "B(2, :) = " << B.row(2) << endl;

	// dot products
	cout << "B(:,1)' * B(:,2) = " << Vector_Proxy::dot(B.col(1), B.col(2)) << endl;

	// matrix-matrix mult

	return 0;
}