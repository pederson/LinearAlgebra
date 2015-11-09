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

Vector Matrix::operator*(const Vector & v)
{
	if (m_ncols != v.rows())
	{
		throw "Matrix dimensions do not match!";
	}

	Vector out(m_mrows);
	out.fill(0);
	for (auto i=0; i<m_ncols; i++)
	{
		out += col(i)*v(i);
	}

	return out;
}



#include <typeinfo>

int main(int argc, char * argv[]){

	unsigned int nrows=20, ncols = 5;

	//**************** MATRIX TESTS *********************//
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
	Matrix E = B(0, 2, 0, 4);
	Matrix F = B(14, 18, 2, 4);
	Matrix G = E*F;
	cout << "E = " << E << endl;
	cout << "F = " << F << endl;
	cout << "G = " << G << endl;

	//***************************************************//




	//**************** VECTOR TESTS *********************//
	Vector v1(5);
	for (auto i=0; i<5; i++)
	{
		v1(i) = i;
	}

	// assignment operator
	Vector v2 = v1;
	cout << "v2 = " << v2 << endl;

	v2.print_summary();

	// // subvector extraction
	Vector v3 = v2(1, 2);
	cout << "v3 = " << v3 << endl;

	// // subvector assignment
	Vector v4(2);
	v4.fill(-10);
	v2(0, 1) = v4;
	cout << "v2 = " << v2 << endl;

	// // subvector assignment via vector proxy
	// // B(0, 1, 0, 1) = S(3, 4, 0, 1);
	// // cout << "B = " << B << endl;

	// scalar multiplication, subtraction, addition, division
	Vector v5 = (((v2*0.5)/2.0) + 10) - 2.5;
	cout << "v2 = " << v2 << endl;
	cout << "v5 = " << v5 << endl;

	// vector proxy
	cout << " v2(1:2) = " << v2(1,2) << endl;
	cout << "v2(2, :) = " << v2.row(2) << endl;

	// dot products
	cout << "v2' * v5= " << Vector::dot(v2, v5) << endl;

	// vector-matrix mult
	v5.transpose();
	Vector v6 = v5*F;
	cout << "v5 = " << v5 << endl;
	cout << "F = " << F << endl;
	cout << "v6 = " << v6 << endl;

	// matrix-vector mult
	v6.transpose();
	cout << "F = " << F << endl;
	cout << "v6 = " << v6 << endl;
	Vector v7 = F*v6;
	cout << "v7 = " << v7 << endl;


	//***************************************************//


	return 0;
}