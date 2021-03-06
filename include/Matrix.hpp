#ifndef _MATRIX_H
#define _MATRIX_H

#include <math.h>
#include <iostream>
#include <fstream>
#include <istream>
#include <sstream>
#include <string>
#include <stdlib.h>
#include <random>
#include <chrono>

#include <array>
#include <utility>


#include "Tensor.hpp"
#include "VectorTools.hpp"
#include "AbstractMatrix.hpp"




namespace libra{
	namespace matrix{

		template <typename MatrixType>
		struct MatrixRow : public vector::VectorFunctors<MatrixRow<MatrixType>>,
						   public vector::VectorAssignment<MatrixRow<MatrixType>>{
		private:
			typedef std::remove_reference_t<
						 decltype(std::declval<MatrixType>()(0, 0))
						 > 		
						 matrix_value_type;
		public:

			typedef MatrixRow<MatrixType> 					SelfType;
			typedef vector::VectorAssignment<SelfType> 		AssignmentType;


			using AssignmentType::operator=;


			MatrixRow(MatrixType & m, size_type row)
			: mMat(&m), mRow(row) {};

			size_type size() const {return mMat->cols();};

			void resize(size_type l) {};


			template <bool is_const>
			class matrix_row_iterator{
			public:
				typedef std::remove_reference_t<decltype(std::declval<MatrixType>()(0, 0))> 		orig_value_type;

				typedef matrix_row_iterator					self_type;
				typedef std::ptrdiff_t 						difference_type;
			    typedef typename std::conditional<is_const, 
			    								  typename std::add_const<orig_value_type>::type, 
			    								  orig_value_type>::type 						
			    								  			value_type;
			    typedef value_type &			  			reference;
			    typedef value_type *						pointer;
			    typedef std::random_access_iterator_tag		iterator_category;

				// construction
				matrix_row_iterator(MatrixType * m, size_type row, size_type col)
				: mMat(m), mRow(row), mCol(col){};

				// copy assignment
				matrix_row_iterator & operator=(const matrix_row_iterator & cit){
					matrix_row_iterator i(cit);
					std::swap(i,*this);
					return *this;
				}

				// rvalue dereferencing
				// pointer operator->() {return *mMat(mRow, mCol);};
				// reference operator*(){return mMat->operator()(mRow, mCol);};

				pointer operator->() const {return mMat->operator()(mRow, mCol);};
				reference operator*() const {return mMat->operator()(mRow, mCol);};

				// increment operators
				self_type operator++(){
					mCol++;
					return *this;
				}
				self_type operator++(int blah){
					mCol++;
					return *this;
				}

				// decrement operators
				self_type operator--(){
					mCol--;
					return *this;
				}
				self_type operator--(int blah){
					mCol--;
					return *this;
				}

				// scalar arithmetic operators
				self_type operator+(int n){
					mCol = mCol + n;
					return *this;
				}
				self_type operator-(int n){
					mCol = mCol - n;
					return *this;
				}
				int operator-(const self_type & b) const {
					return mCol - b.mCol;
				}

				// equivalence operators
				bool operator!=(const self_type & leaf) const {return mCol != leaf.mCol;};
				bool operator==(const self_type & leaf) const {return mCol == leaf.mCol;};

				// relational operators
				bool operator>(const self_type & leaf) const {return mCol > leaf.mCol;};
				bool operator>=(const self_type & leaf) const {return mCol >= leaf.mCol;};
				bool operator<(const self_type & leaf) const {return mCol < leaf.mCol;};
				bool operator<=(const self_type & leaf) const {return mCol <= leaf.mCol;};


				// compound assignment operators
				self_type operator+=(int n){
					mCol += n;
					return *this;
				}
				self_type operator-=(int n){
					mCol -= n;
					return *this;
				}


				// offset dereference operator
				// reference operator[](int n){
				// 	return *mMat(mRow, n);
				// }

				reference operator[](int n) const {
					return mMat->operator()(mRow, n);
				}
			private:
				size_type mCol;
				size_type mRow;
				typename std::conditional<is_const, typename std::add_const<MatrixType>::type, MatrixType>::type * mMat;
			};


			typedef matrix_row_iterator<true> const_iterator;
			typedef matrix_row_iterator<false> iterator;


			iterator begin() {return iterator(mMat, mRow, 0);};
			iterator end()	 {return iterator(mMat, mRow, size());};

			const_iterator cbegin() const {return const_iterator(mMat, mRow, 0);};
			const_iterator cend() const	 {return const_iterator(mMat, mRow, size());};


			matrix_value_type & operator()(int n) {return mMat->operator()(mRow, n);};
		private:
			MatrixType * 	mMat;
			size_type 		mRow;
		};







		template <typename MatrixType>
		struct MatrixCol : public vector::VectorFunctors<MatrixRow<MatrixType>>,
						   public vector::VectorAssignment<MatrixCol<MatrixType>>{
		private:
			typedef std::remove_reference_t<
						 decltype(std::declval<MatrixType>()(0, 0))
						 > 		
						 matrix_value_type;
		public:

			typedef MatrixCol<MatrixType> 					SelfType;
			typedef vector::VectorAssignment<SelfType> 		AssignmentType;


			using AssignmentType::operator=;


			MatrixCol(MatrixType & m, size_type col)
			: mMat(&m), mCol(col) {};

			size_type size() const {return mMat->rows();};

			void resize(size_type l) {};

			template <bool is_const>
			class matrix_col_iterator{
			public:
				typedef std::remove_reference_t<decltype(std::declval<MatrixType>()(0, 0))> 		orig_value_type;

				typedef matrix_col_iterator					self_type;
				typedef std::ptrdiff_t 						difference_type;
			    typedef typename std::conditional<is_const, 
			    								  typename std::add_const<orig_value_type>::type, 
			    								  orig_value_type>::type 						
			    								  			value_type;
			    typedef value_type &			  			reference;
			    typedef value_type *						pointer;
			    typedef std::random_access_iterator_tag		iterator_category;

				// construction
				matrix_col_iterator(MatrixType * m, size_type row, size_type col)
				: mMat(m), mRow(row), mCol(col){};

				// copy assignment
				matrix_col_iterator & operator=(const matrix_col_iterator & cit){
					matrix_col_iterator i(cit);
					std::swap(i,*this);
					return *this;
				}

				// rvalue dereferencing
				// pointer operator->() {return *mMat(mRow, mCol);};
				// reference operator*(){return mMat->operator()(mRow, mCol);};

				pointer operator->() const {return mMat->operator()(mRow, mCol);};
				reference operator*() const {return mMat->operator()(mRow, mCol);};

				// increment operators
				self_type operator++(){
					mRow++;
					return *this;
				}
				self_type operator++(int blah){
					mRow++;
					return *this;
				}

				// decrement operators
				self_type operator--(){
					mRow--;
					return *this;
				}
				self_type operator--(int blah){
					mRow--;
					return *this;
				}

				// scalar arithmetic operators
				self_type operator+(int n){
					mRow = mRow + n;
					return *this;
				}
				self_type operator-(int n){
					mRow = mRow - n;
					return *this;
				}
				int operator-(const self_type & b) const {
					return mRow - b.mRow;
				}

				// equivalence operators
				bool operator!=(const self_type & leaf) const {return mRow != leaf.mRow;};
				bool operator==(const self_type & leaf) const {return mRow == leaf.mRow;};

				// relational operators
				bool operator>(const self_type & leaf) const {return mRow > leaf.mRow;};
				bool operator>=(const self_type & leaf) const {return mRow >= leaf.mRow;};
				bool operator<(const self_type & leaf) const {return mRow < leaf.mRow;};
				bool operator<=(const self_type & leaf) const {return mRow <= leaf.mRow;};


				// compound assignment operators
				self_type operator+=(int n){
					mRow += n;
					return *this;
				}
				self_type operator-=(int n){
					mRow -= n;
					return *this;
				}


				// offset dereference operator
				// reference operator[](int n){
				// 	return *mMat(mRow, n);
				// }

				reference operator[](int n) const {
					return mMat->operator()(n, mCol);
				}
			private:
				size_type mCol;
				size_type mRow;
				typename std::conditional<is_const, typename std::add_const<MatrixType>::type, MatrixType>::type * mMat;
			};


			typedef matrix_col_iterator<true> const_iterator;
			typedef matrix_col_iterator<false> iterator;


			iterator begin() {return iterator(mMat, 0, mCol);};
			iterator end()	 {return iterator(mMat, size(), mCol);};

			const_iterator cbegin() const {return const_iterator(mMat, 0, mCol);};
			const_iterator cend() const	 {return const_iterator(mMat, size(), mCol);};

			matrix_value_type & operator()(int n) {return mMat->operator()(n, mCol);};
		private:
			MatrixType * 	mMat;
			size_type 		mCol;
		};







		// fixed-length resize policy -- do nothing
		template <size_type rows_at_compile, size_type cols_at_compile>
		struct MatrixResizePolicy{
			template <typename MatrixType>
			static void resize(MatrixType & v, size_type r, size_type c) {};
		};
		// dynamic columns size resize policy
		template <size_type rows_at_compile>
		struct MatrixResizePolicy<rows_at_compile, dynamic_size>{
			template <typename MatrixType>
			static void resize(MatrixType & v, size_type r, size_type c) {
				for (auto i=0; i<r; i++)v[i].resize(c);};
		};

		// dynamic columns size resize policy
		template <>
		struct MatrixResizePolicy<dynamic_size, dynamic_size>{
			template <typename MatrixType>
			static void resize(MatrixType & v, size_type r, size_type c) {
				v.resize(r);
				for (auto i=0; i<r; i++)v[i].resize(c);};
		};



		// matrix class definition
		template <typename scalar_type, size_type rows_at_compile, size_type cols_at_compile>
		class Matrix : public Table<scalar_type, 2, rows_at_compile, cols_at_compile>{
		public:
			typedef Table<scalar_type, 2, rows_at_compile, cols_at_compile> 	BaseType;
			typedef Matrix<scalar_type, rows_at_compile, cols_at_compile> 		SelfType;

			// inherit the base class constructors
			using BaseType::BaseType;

			// this nasty-looking code simply allows the use of vector and array braces initialization
			// An input braces-initializer is converted to a row-major ordered matrix
			// e.g. Matrix<int, 2, 3> mat = {1,2,3,4,5,6} becomes:		|1, 2, 3|
			// 															|4, 5, 6|
			// template <typename... Args>
		 //    Matrix(Args &&... args) : BaseType({(std::forward<Args>(args))...}) {};


			// explicitly define an initializer list constructor
			// An input braces-initializer is converted to a row-major ordered matrix
			// e.g. Matrix<int, 2, 3> mat = {{1,2,3},{4,5,6}} becomes:	|1, 2, 3|
			// 															|4, 5, 6|
			Matrix(std::initializer_list<std::initializer_list<scalar_type>> il){
				auto it = std::begin(il);
				MatrixResizePolicy<rows_at_compile, cols_at_compile>::resize(*this, il.size(), (*it).size());
				
				auto itrme = std::begin(*this);
				while (it != std::end(il)){
					auto it2 = std::begin(*it);
					auto itcme = std::begin(*itrme);
					while (it2 != std::end(*it)){
						(*itcme) = *it2;
						it2++;
						itcme++;
					}
					// (*itme) = *it;
					it++;
					itrme++;
				}
			}

			Matrix(){};

			Matrix(size_type r, size_type c){
				MatrixResizePolicy<rows_at_compile, cols_at_compile>::resize(*this, r, c);
			}


			scalar_type & operator()(size_type i, size_type j){return (*this)[i][j];};
			const scalar_type & operator()(size_type i, size_type j) const {return (*this)[i][j];};



			template <typename T = size_type>
			typename std::enable_if<rows_at_compile != dynamic_size, T>::type 
			rows() const {return rows_at_compile;};

			template <typename T = size_type>
			typename std::enable_if<rows_at_compile == dynamic_size, T>::type 
			rows() const {return this->size();};


			template <typename T = size_type>
			typename std::enable_if<cols_at_compile != dynamic_size, T>::type 
			cols() const {return cols_at_compile;};

			template <typename T = size_type>
			typename std::enable_if<cols_at_compile == dynamic_size, T>::type 
			cols() const {return (rows() > 0 ? (*this)[0].size() : 0);};


			// matrix-vector multiplication
			template <typename VectorT, typename VectorT2>
			void vmult(const VectorT & v, VectorT2 && result) const {
				static_assert(type_traits::is_traversable_vector<VectorT>::value, "Matrix-Vector input must be a vector type!");
				static_assert(type_traits::is_assignable_vector<VectorT2>::value, "Matrix-Vector input must be a vector type!");

				// std::cout << "at begin of matmult" << std::endl;
				scalar_type rowsum;
				auto rit = this->cbegin();
				auto resit = result.begin();//std::begin(result);

				// std::cout << "at middle of matmult" << std::endl;

				while (rit != this->cend() && resit != result.end()){
					rowsum = vector::inner_product(v, (*rit));
					(*resit) = rowsum;
					rit++; resit++;
				}
			};


			// (matrix-transpose)-vector multiplication
			template <typename VectorT, typename VectorT2>
			void Tmult(const VectorT & v, VectorT2 && result) {
				static_assert(type_traits::is_traversable_vector<VectorT>::value, "Matrix-Vector input must be a vector type!");
				static_assert(type_traits::is_assignable_vector<VectorT2>::value, "Matrix-Vector input must be a vector type!");

				// std::cout << "at begin of matmult" << std::endl;
				scalar_type rowsum;
				auto rit = this->cbegin();
				auto resit = result.begin();//std::begin(result);

				// std::cout << "at middle of matmult" << std::endl;
				for (auto i=0; i<this->cols(); i++){
					rowsum = vector::inner_product(v, this->col(i));
					(*resit) = rowsum;
					resit++;
				}
			};




			MatrixRow<SelfType> row(size_type r) {return MatrixRow<SelfType>(*this, r);};
			// MatrixRow<SelfType> row(size_type r) const {return MatrixRow<SelfType>(*this, r);};
			MatrixCol<SelfType> col(size_type c) {return MatrixCol<SelfType>(*this, c);};
			// MatrixCol<SelfType> col(size_type c) const {return MatrixCol<SelfType>(*this, c);};
	
		private:
		};


		template <bool Header = false, typename MatrixT>
		// template <typename scalar_type, size_type rows_at_compile, size_type cols_at_compile>
		void write(const MatrixT & m, std::string dlm = " ", std::ostream & os = std::cout, std::size_t ntabs = 0){

			// static_assert(type_traits::is_matrix<MatrixT>::value, "Input must be a matrix type!");

			for (auto i=0; i<ntabs; i++) os << "\t" ;
			if (Header) os << "<Matrix>" << std::endl;
			os << std::scientific;

			for (auto r = 0; r<m.rows(); r++){
				for (auto i=0; i<ntabs+1; i++) os << "\t" ;

				// first item of the row
				os << m(r,0) ;

				// remaining items in the row get a delimiter
				for (auto c = 1; c<m.cols(); c++){
					os << dlm << m(r,c) ;
				}
				os << std::endl;
			}

			for (auto i=0; i<ntabs; i++) os << "\t" ;
			if (Header) os << "</Matrix>" << std::endl;
			
			return;
		};






		

	} // end namespace matrix
} // end namespace libra








































// using namespace libra;
// // *******************************************************************************
// namespace libra{

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




// *******************************************************************************


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
			std::cout << "Vectors must be of same size to dot!" << std::endl;
			throw -1;
		}

		return result;
	}

	Vector_Proxy & operator=(const Vector & vct);

	Vector_Proxy & operator-=(const Vector & vct);

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

	void operator-=(Vector_Proxy & vctp)
	{
		for (auto i=0; i<m_length; i++)
		{
			m_dataptr[i*m_stride] -= vctp.m_dataptr[i*vctp.m_stride];
		}
	}

	void fill(double d){
		for (auto i=0; i<m_length; i++)
		{
			m_dataptr[i*m_stride] = d;
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




// *******************************************************************************

class Matrix : public AbstractMatrix
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

	// Matrix assignment
	Matrix & operator=(const Matrix & A)
	{

		Matrix out(A);
		swap(*this, out);
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

	// Matrix-Vector transpose multiplication
	Vector Tmult(const Vector & v) const;

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

	// Matrix-SparseMatrix addition
	Matrix operator+(const SparseMatrix & mtx) const;

	// Matrix-SparseMatrix subtraction
	Matrix operator-(const SparseMatrix & mtx) const;

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
	//double norm() const {return 0;};

	const double * data() const {return m_data;};

	void transpose();

	void dlmwrite(std::string filename, std::string delimiter=",") const{
		std::ofstream file(filename);
		const char * dlm = delimiter.c_str();

		//file << std::endl;
		file << std::scientific;
		for (auto i=0; i<m_mrows; i++)
		{
			file << m_data[m_mrows*0 + i];
			for (auto j=1; j<m_ncols; j++)
			{
				file << dlm << m_data[m_mrows*j + i];
			}
			file << std::endl;
		}
	}

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



// *******************************************************************************

// vector is derived from matrix. Has either 1 column or 1 row
class Vector : public Matrix, public AbstractVector{
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


	// class const_iterator;
	// class iterator{
	// public:
	// 	friend class const_iterator;
	// 	typedef iterator 							self_type;
	// 	typedef std::ptrdiff_t 						difference_type;
	//     typedef typename IteratorT::value_type 		value_type;
	//     typedef typename IteratorT::reference 		reference;
	//     typedef typename IteratorT::pointer 		pointer;
	//     typedef typename IteratorT::iterator_category	iterator_category;

	// 	// construction
	// 	iterator(VectorView & sv, IteratorT it)
	// 	: mSV(sv)
	// 	, mIt(it){};

	// 	// copy assignment
	// 	iterator & operator=(const iterator & cit){
	// 		iterator i(cit);
	// 		std::swap(i,*this);
	// 		return *this;
	// 	}

	// 	// rvalue dereferencing
	// 	pointer operator->() {return mIt.operator->();};
	// 	reference operator*(){ return *mIt;};

	// 	// increment operators
	// 	self_type operator++(){
	// 		mIt++;
	// 		return *this;
	// 	}
	// 	self_type operator++(int blah){
	// 		mIt++;
	// 		return *this;
	// 	}

	// 	// decrement operators
	// 	self_type operator--(){
	// 		mIt--;
	// 		return *this;
	// 	}
	// 	self_type operator--(int blah){
	// 		mIt--;
	// 		return *this;
	// 	}

	// 	// scalar arithmetic operators
	// 	self_type operator+(int n){
	// 		mIt = mIt + n;
	// 		return *this;
	// 	}
	// 	self_type operator-(int n){
	// 		mIt = mIt - n;
	// 		return *this;
	// 	}
	// 	int operator-(const self_type & b) const {
	// 		return mIt - b.mIt;
	// 	}

	// 	// equivalence operators
	// 	bool operator!=(const self_type & leaf) const {return mIt != leaf.mIt;};
	// 	bool operator==(const self_type & leaf) const {return mIt == leaf.mIt;};

	// 	// relational operators
	// 	bool operator>(const self_type & leaf) const {return mIt > leaf.mIt;};
	// 	bool operator>=(const self_type & leaf) const {return mIt >= leaf.mIt;};
	// 	bool operator<(const self_type & leaf) const {return mIt < leaf.mIt;};
	// 	bool operator<=(const self_type & leaf) const {return mIt <= leaf.mIt;};


	// 	// compound assignment operators
	// 	self_type operator+=(int n){
	// 		mIt += n;
	// 		return *this;
	// 	}
	// 	self_type operator-=(int n){
	// 		mIt -= n;
	// 		return *this;
	// 	}


	// 	// offset dereference operator
	// 	reference operator[](int n){
	// 		return mIt[n];
	// 	}
	// private:
	// 	VectorView & mSV;
	// 	IteratorT mIt;
	// };

	// iterator begin() {return iterator(*this, mBegin);};
	// iterator end()	 {return iterator(*this, mEnd);};




	// class const_iterator{
	// public:
	// 	typedef typename VectorT::const_iterator 		subiterator_type;
		
	// 	typedef const_iterator 							self_type;
	// 	typedef typename subiterator_type::difference_type		difference_type;
	//     typedef typename subiterator_type::value_type 	value_type;
	//     typedef typename subiterator_type::reference 	reference;
	//     typedef typename subiterator_type::pointer 		pointer;
	//     typedef typename subiterator_type::iterator_category	iterator_category;

	// 	// construction
	// 	const_iterator(const VectorView & sv, subiterator_type it)
	// 	: mSV(sv)
	// 	, mIt(it){};

	// 	// conversion of iterator to const_iterator
	// 	const_iterator(const iterator & it)
	// 	: mSV(it.mSV)
	// 	, mIt(it.mIt) {};


	// 	const_iterator & operator=(const const_iterator & cit){
	// 		const_iterator i(cit);
	// 		std::swap(i,*this);
	// 		return *this;
	// 	}

	// 	// rvalue dereferencing
	// 	pointer operator->() {return mIt.operator->();};
	// 	reference operator*(){ return *mIt;};

	// 	// increment operators
	// 	self_type operator++(){
	// 		mIt++;
	// 		return *this;
	// 	}
	// 	self_type operator++(int blah){
	// 		mIt++;
	// 		return *this;
	// 	}

	// 	// decrement operators
	// 	self_type operator--(){
	// 		mIt--;
	// 		return *this;
	// 	}
	// 	self_type operator--(int blah){
	// 		mIt--;
	// 		return *this;
	// 	}

	// 	// scalar arithmetic operators
	// 	self_type operator+(int n){
	// 		mIt = mIt + n;
	// 		return *this;
	// 	}
	// 	self_type operator-(int n){
	// 		mIt = mIt - n;
	// 		return *this;
	// 	}
	// 	int operator-(const self_type & b) const {
	// 		return mIt - b.mIt;
	// 	}

	// 	// equivalence operators
	// 	bool operator!=(const self_type & leaf) const {return mIt != leaf.mIt;};
	// 	bool operator==(const self_type & leaf) const {return mIt == leaf.mIt;};

	// 	// relational operators
	// 	bool operator>(const self_type & leaf) const {return mIt > leaf.mIt;};
	// 	bool operator>=(const self_type & leaf) const {return mIt >= leaf.mIt;};
	// 	bool operator<(const self_type & leaf) const {return mIt < leaf.mIt;};
	// 	bool operator<=(const self_type & leaf) const {return mIt <= leaf.mIt;};


	// 	// compound assignment operators
	// 	self_type operator+=(int n){
	// 		mIt += n;
	// 		return *this;
	// 	}
	// 	self_type operator-=(int n){
	// 		mIt -= n;
	// 		return *this;
	// 	}


	// 	// offset dereference operator
	// 	reference operator[](int n){
	// 		return mIt[n];
	// 	}
	// private:
	// 	const VectorView & mSV;
	// 	subiterator_type mIt;
	// };

	// const_iterator cbegin() const {return const_iterator(*this, vector_const_iterator(mBegin));};
	// const_iterator cend() const	 {return const_iterator(*this, vector_const_iterator(mEnd));};

	// // these are required by std::cbegin()/cend()
	// const_iterator begin() const {return cbegin();};
	// const_iterator end() const	 {return cend();};

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

	// Vector-Vector elementwise multiply
	Vector elem_mult(const Vector & vct) const
	{
		if (m_len != vct.m_len)
		{
			std::cout << "Vector dimensions do not match!" << std::endl;
			throw -1;
		}

		Vector out(m_len);
		for (auto i=0; i<m_len; i++)
		{
			out(i) = m_data[i] * vct.m_data[i];
		}

		return out;
	}

	// Vector-Vector elementwise division
	Vector elem_div(const Vector & vct) const
	{
		if (m_len != vct.m_len)
		{
			std::cout << "Vector dimensions do not match!" << std::endl;
			throw -1;
		}

		Vector out(m_len);
		for (auto i=0; i<m_len; i++)
		{
			out(i) = m_data[i] / vct.m_data[i];
		}

		return out;
	}

	Vector operator+(const SparseVector & vct) const;

	Vector operator-(const SparseVector & vct) const;

	SparseVector get_support(const SparseVector & vct) const;

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

	void size(std::size_t & rows, std::size_t & cols) const {rows = m_mrows; cols = m_ncols; return;};

	std::size_t length() const {return m_len;};

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

// identity matrix
Matrix eye(unsigned int rows, unsigned int cols)
{
	Matrix out(rows, cols);
	unsigned int sz = std::min(rows,cols);
	out.fill(0);
	for (auto i=0; i<sz; i++)
	{
		out(i,i) = 1;
	}
	return out;
}

// hilbert matrix
// entries are H(i,j) = 1/(i+j+1) where i, j start at 0
Matrix hilb(unsigned int rows)
{
	Matrix out(rows, rows);
	for (auto i=0; i<rows; i++)
	{
		for (auto j=0; j<rows; j++)
		{
			out(i,j) = 1.0/(i+j+1);
		}
	}
	return out;
}

// random matrix uniformly distributed [0,1]
Matrix randmat(unsigned int rows, unsigned int cols)
{
	// seed
	// std::default_random_engine generator;
	unsigned int seed1 = std::chrono::system_clock::now().time_since_epoch().count();
	std::minstd_rand0 generator(seed1);
	std::uniform_real_distribution<double> distrib(0.0,1.0);

	Matrix out(rows, cols);
	for (auto i=0; i<rows; i++){
		for (auto j=0; j<cols; j++){
			out(i,j) = distrib(generator);
		}
	}

	return out;
}

// random matrix normally distributed
Matrix randmatn(unsigned int rows, unsigned int cols)
{
	// std::default_random_engine generator;
	unsigned int seed1 = std::chrono::system_clock::now().time_since_epoch().count();
	std::minstd_rand0 generator(seed1);
	std::normal_distribution<double> distrib(0.0,1.0);

	Matrix out(rows, cols);
	for (auto i=0; i<rows; i++){
		for (auto j=0; j<cols; j++){
			out(i,j) = distrib(generator);
		}
	}

	return out;
}


Matrix dlmread(std::string filename, std::string delimiter=",", unsigned int headerlines=0){
	// open the file as a stream
	std::ifstream file(filename);
	const char * dlm = delimiter.c_str();

	// read the number of rows
	unsigned int nrows = 0;
	std::string line;
	for (auto i=0; i<headerlines; i++) std::getline(file, line);
	while (getline(file,line)) nrows++;

	// read the number of columns
	file.clear();
	file.seekg(0, std::ios::beg);
	for (auto i=0; i<headerlines; i++) std::getline(file, line);
	std::getline(file,line);
	std::stringstream ss(line);
	std::string field;
	unsigned int ncols=0;
	//std::cout << "first line: " << line << std::endl;
	while(std::getline(ss,field,*dlm)) ncols++;

	// std::cout << "(" << nrows << " x " << ncols << ") " << std::endl;
	// std::cout << "delimiter: " << dlm << std::endl;

	// read the data
	Matrix out(nrows,ncols);
	file.clear();
	file.seekg(0, file.beg);
	for (auto i=0; i<headerlines; i++) std::getline(file, line);
	unsigned int i=0,j;
	while (std::getline(file,line)){
		std::stringstream ss(line);
		std::string field;
		j=0;
		while(std::getline(ss,field,*dlm)){
			std::stringstream fs(field);
			double f=0.0;
			fs >> f;
			out(i,j) = f;
			j++;
		}
		i++;
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

Vector abs(const Vector & vct){
	Vector out(vct);
	for (auto i=0; i<out.length(); i++) out(i) = fabs(vct(i));
	return out;
}

Vector sign(const Vector & vct){
	Vector out(vct);
	for (auto i=0; i<out.length(); i++) out(i) = sgn(vct(i));
	return out;
}

double norm_1(const Vector & vct){
	double res=0;
	for (auto i=0; i<vct.length(); i++) res += abs(vct(i));
	return res;
}

double norm_2(const Vector & vct){
	double res=0;
	for (auto i=0; i<vct.length(); i++) res += vct(i)*vct(i);
	return sqrt(res);
}

double norm_inf(const Vector & vct){
	double res=0;
	for (auto i=0; i<vct.length(); i++) res = std::max(res, fabs(vct(i)));
	return res;
}

// random vector uniformly distributed [0,1]
Vector randvec(unsigned int length)
{
	Matrix m = randmat(length,1);
	Vector out = m.col(0);
	return out;
}

// random vector normally distributed
Vector randvecn(unsigned int length)
{
	Matrix m=randmatn(length,1);
	Vector out = m.col(0);
	return out;
}

unsigned int argmin(const Vector & vct){
	unsigned int res=0;
	double mv = vct(0);
	for (auto i=0; i<vct.length(); i++){
		if (vct(i) < mv){
			res = i;
			mv = vct(i);
		}
	}
	return res;
}

unsigned int argmax(const Vector & vct){
	unsigned int res=0;
	double mv = vct(0);
	for (auto i=0; i<vct.length(); i++){
		if (vct(i) > mv){
			res = i;
			mv = vct(i);
		}
	}
	return res;
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

Vector operator+(const Vector_Proxy & vct, const Vector & v)
{
	Vector out(vct);
	for (auto i=0; i<vct.length(); i++) out(i) += v(i);
	return out;
}

Vector operator-(const Vector_Proxy & vct, const Vector & v)
{
	Vector out(vct);
	for (auto i=0; i<vct.length(); i++) out(i) -= v(i);
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


Vector_Proxy & Vector_Proxy::operator-=(const Vector & vct)
{
	for (auto i=0; i<m_length; i++)
	{
		*(m_dataptr+i*m_stride) -= vct(i);
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


Vector Matrix::Tmult(const Vector & v) const
{
	if (m_ncols != v.rows())
	{
		throw "Matrix dimensions do not match!";
	}

	Vector out(m_mrows);

	out.fill(0);
	for (auto i=0; i<m_mrows; i++)
	{
		Vector r = row(i);
		r.transpose();
		// out += col(i)*v(i);
		out += r*v(i);
	}

	return out;
}


// }




#endif