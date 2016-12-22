#ifndef _SPARSEMATRIX_H
#define _SPARSEMATRIX_H

#include <math.h>
#include <iostream>
#include <stdlib.h>
#include <list>
#include <vector>
#include <random>

#include "Matrix.hpp"
#include "SparseVector.hpp"


// *******************************************************************************


class SparseMatrix : public AbstractMatrix
{
public:
	// constructor
	SparseMatrix()
		: m_len 	(0)
		, m_mrows 	(0)
		, m_ncols 	(0)
	{
	}

	// create and allocate a new sparse matrix
	SparseMatrix(unsigned int mrows, unsigned int ncols)
		: m_len 	(mrows*ncols)
		, m_mrows 	(mrows)
		, m_ncols 	(ncols)
	{
		m_row_ptr.assign(mrows+1,m_data.begin());
	}

	// create a sparse vector initialized by existing data
	SparseMatrix(unsigned int mrows, unsigned int ncols, 
				 unsigned int num_vals, const unsigned int * row,
				 const unsigned int * col, const double * data)
		: m_len 	(mrows*ncols)
		, m_mrows 	(mrows)
		, m_ncols 	(ncols)
	{
		m_row_ptr.assign(mrows+1, m_data.begin());
		for (auto i=0; i<num_vals; i++) set(row[i], col[i], data[i]);
	}

	// copy constructor
	SparseMatrix(const SparseMatrix & sm)
		: m_len 	(sm.m_len)
		, m_mrows 	(sm.m_mrows)
		, m_ncols 	(sm.m_ncols)
	{
		m_row_ptr.assign(m_mrows+1, m_data.begin());
		// add each component of mtx one by one
		for (auto i=0; i<sm.m_mrows; i++){
			auto it = sm.m_row_ptr[i];
			while (it != sm.m_row_ptr[i+1] && it !=sm.m_data.end()){
				set(i,it->first, it->second);
				it++;
			}
		}
	}

	// destructor
	~SparseMatrix()
	{
	}

	friend void swap(SparseMatrix & sm1, SparseMatrix & sm2)
	{
		using std::swap;
		swap(sm1.m_len, sm2.m_len);
		swap(sm1.m_mrows, sm2.m_mrows);
		swap(sm1.m_ncols, sm2.m_ncols);
		swap(sm1.m_data, sm2.m_data);
		swap(sm1.m_row_ptr, sm2.m_row_ptr);

		std::cout << "SWAP OPERATION ISNT WORKING ON SPARSE MATRICES" << std::endl;
		throw -1;
	}


	// SparseMatrix-Vector product
	Vector operator*(const Vector & vct) const{
		Vector out(m_mrows);
		out.fill(0);
		
		for (auto i=0; i<vct.rows(); i++){
			double rowsum = 0.0;
			auto it = m_row_ptr[i];
			while (it != m_row_ptr[i+1]){
				rowsum += it->second*vct(it->first);
				it++;
			}
			out(i) = rowsum;
		}
		return out;
	}


	// SparseMatrix-SparseMatrix addition
	SparseMatrix operator+(const SparseMatrix & mtx) const
	{
		if (m_mrows != mtx.m_mrows || m_ncols != mtx.m_ncols)
		{
			std::cout << "SparseMatrix dimensions do not match!" << std::endl;
			throw -1;
		}

		SparseMatrix out(*this);

		// add each component of mtx one by one
		for (auto i=0; i<mtx.m_mrows; i++){
			auto it = mtx.m_row_ptr[i];
			while (it != mtx.m_row_ptr[i+1] && it !=mtx.m_data.end()){
				out.add(i,it->first, it->second);
				it++;
			}
		}

		return out;
	}

	// SparseMatrix-SparseMatrix subtraction
	SparseMatrix operator-(const SparseMatrix & mtx) const
	{
		if (m_mrows != mtx.m_mrows || m_ncols != mtx.m_ncols)
		{
			std::cout << "SparseMatrix dimensions do not match!" << std::endl;
			throw -1;
		}

		SparseMatrix out(*this);

		// add each component of mtx one by one
		for (auto i=0; i<mtx.m_mrows; i++){
			auto it = mtx.m_row_ptr[i];
			while (it != mtx.m_row_ptr[i+1] && it !=mtx.m_data.end()){
				out.add(i,it->first, -it->second);
				it++;
			}
		}

		return out;
	}

	// SparseMatrix-SparseMatrix addition shorthand
	SparseMatrix & operator+=(const SparseMatrix & mtx)
	{
		if (m_mrows != mtx.m_mrows || m_ncols != mtx.m_ncols)
		{
			std::cout << "SparseMatrix dimensions do not match!" << std::endl;
			throw -1;
		}


		// add each component of mtx one by one
		for (auto i=0; i<mtx.m_mrows; i++){
			auto it = mtx.m_row_ptr[i];
			while (it != mtx.m_row_ptr[i+1] && it !=mtx.m_data.end()){
				add(i,it->first, it->second);
				it++;
			}
		}

		return *this;
	}

	// SparseMatrix-SparseMatrix subtraction shorthand
	SparseMatrix & operator-=(const SparseMatrix & mtx)
	{
		if (m_mrows != mtx.m_mrows || m_ncols != mtx.m_ncols)
		{
			std::cout << "SparseMatrix dimensions do not match!" << std::endl;
			throw -1;
		}


		// add each component of mtx one by one
		for (auto i=0; i<mtx.m_mrows; i++){
			auto it = mtx.m_row_ptr[i];
			while (it != mtx.m_row_ptr[i+1] && it !=mtx.m_data.end()){
				add(i,it->first, -it->second);
				it++;
			}
		}

		return *this;
	}

	// transpose
	SparseMatrix operator~() const
	{
		SparseMatrix sv(*this);
		sv.transpose();
		return sv;
	}

	Matrix operator+(const Matrix & mtx) const;

	Matrix operator-(const Matrix & mtx) const;


	// scalar multiplication
	SparseMatrix operator*(double val) const
	{
		SparseMatrix out(*this);
		for (auto it=out.m_data.begin(); it!=out.m_data.end(); it++){
			it->second *=val;
		}
		return out;
	}

	// scalar division
	SparseMatrix operator/(double val) const
	{
		SparseMatrix out(*this);
		for (auto it=out.m_data.begin(); it!=out.m_data.end(); it++){
			it->second /=val;
		}
		return out;
	}

	// scalar addition
	SparseMatrix operator+(double val) const
	{
		SparseMatrix out(*this);
		for (auto it=out.m_data.begin(); it!=out.m_data.end(); it++){
			it->second +=val;
		}
		return out;
	}

	// scalar subtraction
	SparseMatrix operator-(double val) const
	{
		SparseMatrix out(*this);
		for (auto it=out.m_data.begin(); it!=out.m_data.end(); it++){
			it->second -=val;
		}
		return out;
	}

	// scalar multiplication
	SparseMatrix & operator*=(double val)
	{
		for (auto it=m_data.begin(); it!=m_data.end(); it++){
			it->second*=val;
		}
		return *this;
	}

	// scalar division
	SparseMatrix & operator/=(double val)
	{
		for (auto it=m_data.begin(); it!=m_data.end(); it++){
			it->second/=val;
		}
		return *this;
	}

	// scalar addition
	SparseMatrix & operator+=(double val)
	{
		for (auto it=m_data.begin(); it!=m_data.end(); it++){
			it->second+=val;
		}
		return *this;
	}

	// scalar subtraction
	SparseMatrix & operator-=(double val)
	{
		for (auto it=m_data.begin(); it!=m_data.end(); it++){
			it->second-=val;
		}
		return *this;
	}

	// // value access
	// // WARNING: will insert new element if it doesn't already exist
	// double & operator()(unsigned int i)
	// {
	// 	return m_data[i];
	// }

	// // const value access
	// double operator()(unsigned int i) const
	// {
	// 	int r = m_data.count(i);
	// 	if (r==0) return m_default_value;
	// 	return m_data.at(i);
	// }



	// access data
	const std::list<std::pair<unsigned int, double>> & data() const{return m_data;};
	const std::vector<std::list<std::pair<unsigned int, double>>::iterator> & row_ptr() const{return m_row_ptr;};

	// number of "nonzeros"
	std::size_t nnz() const {return m_data.size();};

	// return vector of "nonempty" index pairs
	std::vector<std::pair<unsigned int, unsigned int>> support() const{
		std::vector<std::pair<unsigned int, unsigned int>> v;

		for (auto i=0; i<m_mrows; i++){
			auto it = m_row_ptr[i];
			while (it != m_row_ptr[i+1] && it !=m_data.end()){
				v.push_back(std::pair<unsigned int, unsigned int>(i,it->first));
				it++;
			}
		}
		return v;
	}

	// set index i,j to val
	void set(unsigned int i, unsigned int j, double val){
		auto row_it = m_row_ptr[i];
		auto it=m_row_ptr[i];

		// find the minimum row that has the same iterator
		unsigned int minrow = i;
		for (auto ro=0; ro<=i; ro++){
			if (m_row_ptr[ro] == m_row_ptr[i]){
				minrow = ro;
				break;

			}
		}
		// std::cout << "minrow: " << minrow << std::endl;

		unsigned int ct=0;
		while(it->first < j && it != m_row_ptr[i+1]){
			it++;
			ct++;
		}
		if (it == m_row_ptr[i+1]) m_data.emplace(it, std::pair<unsigned int, double> (j,val));
		else{
			if (it->first == j){
				it->second = val;
				ct=1;
			}
			else m_data.emplace(it, std::pair<unsigned int, double> (j,val));
		}

		// adjust the row iterator if necessary
		--it;
		if (ct==0){
			for (auto ro=minrow; ro<=i; ro++) m_row_ptr[ro] = it;
		}

	}



	// add val to index i,j
	void add(unsigned int i, unsigned int j, double val){
		auto row_it = m_row_ptr[i];
		auto it=m_row_ptr[i];

		// find the minimum row that has the same iterator
		unsigned int minrow = i;
		for (auto ro=0; ro<=i; ro++){
			if (m_row_ptr[ro] == m_row_ptr[i]){
				minrow = ro;
				break;

			}
		}
		// std::cout << "adding: " << i << ", " << j << " = " << val << std::endl;
		// std::cout << "minrow: " << minrow << std::endl;

		unsigned int ct=0;
		while(it->first < j && it != m_row_ptr[i+1]){
			it++;
			ct++;
		}
		if (it == m_row_ptr[i+1]) m_data.emplace(it, std::pair<unsigned int, double> (j,val));
		else{
			if (it->first == j){
				// std::cout << "adding to value: " << it->second << std::endl;
				it->second += val;
				ct=1;
				// std::cout << "result: " << it->second << std::endl;
			}
			else m_data.emplace(it, std::pair<unsigned int, double> (j,val));
		}

		// adjust the row iterator if necessary
		--it;
		// std::cout << "ct: " << ct << std::endl;
		if (ct==0){
			for (auto ro=minrow; ro<=i; ro++) m_row_ptr[ro] = it;
		}

	}

	// retrieve a row of the matrix as a sparse row vector
	SparseVector row(unsigned int i) const{
		SparseVector out(m_ncols);
		auto it = m_row_ptr[i];
		while (it != m_row_ptr[i+1] && it !=m_data.end()){
			out(it->first) = it->second;
			it++;
		}
		out.transpose();
		return out;
	}

	// multiply matrix by a sparse vector
	Vector operator*(const SparseVector & vct) const{
		Vector out(m_mrows);
		for (auto i=0; i<m_mrows; i++){
			SparseVector r = row(i);
			out(i) = SparseVector::dot(r,vct);
		}
		return out;
	}

	// // get value at index i
	// double get(unsigned int i) const{
	// 	int r = m_data.count(i);
	// 	if (r==0) return m_default_value;
	// 	return m_data.at(i);
	// }

	void size(std::size_t & rows, std::size_t & cols) const{
		rows = m_mrows; cols = m_ncols;
	};

	std::size_t rows() const {return m_mrows;};
	std::size_t cols() const {return m_ncols;};

	// double norm() const{
	// 	double val = double(m_len-m_data.size())*m_default_value;
	// 	for (auto it=m_data.begin(); it!=m_data.end(); it++){
	// 		val += it->second;
	// 	}
	// 	return val;
	// };

	void transpose(){
		SparseMatrix m(m_ncols, m_mrows);
		for (auto i=0; i<m_mrows; i++){
			auto it = m_row_ptr[i];
			while (it != m_row_ptr[i+1]){
				m.set(it->first, i, it->second);
				// std::cout << "setting (" << it->first << ", " << i << ")=" << it->second << std::endl;
				it++;
			}
		}

		auto data = m.data();
		auto rowptr = m.row_ptr();
		for (auto i=0; i<rowptr.size()-1; i++){
			auto it = rowptr[i];
			// std::cout << rowptr[i]->first << std::endl;
			while (it != rowptr[i+1]){
				std::cout << "(" << i << "," << it->first << ") : " << it->second << std::endl;
				it++;
			}
		}

		// std::cout << std::endl;
		// for (auto it=data.begin(); it!=data.end(); it++){
		// 	std::cout << "(" << it->first << "," << it->second << ") " << std::endl;
				
		// }
		// std::cout << m.row_ptr().size() << std::endl;
		// for (auto i=0; i<rowptr.size()-1; i++){
		// 	std::cout << "(" << i << ", " << rowptr[i]->first << "): " << rowptr[i]->second << std::endl;
		// }


		// std::cout << m << std::endl;
		swap(*this, m);

		std::cout << "TRANSPOSING A SPARSE MATRIX DOESNT WORK...SWAP IS MESSED UP" << std::endl;
		throw -1;
	}

	// return a dense matrix version
	Matrix densify() const{
		Matrix out(m_mrows, m_ncols);
		out.fill(0);

		for (auto i=0; i<m_mrows; i++){
			auto it = m_row_ptr[i];
			while (it != m_row_ptr[i+1] && it !=m_data.end()){
				out(i,it->first) = it->second;
				it++;
			}
		}
		return out;
	}


protected:

	std::size_t m_mrows, m_ncols;
	std::size_t m_len;
	std::list<std::pair<unsigned int, double>> m_data;
	std::vector<std::list<std::pair<unsigned int, double>>::iterator> m_row_ptr;

};



std::ostream& operator<<(std::ostream & os, const SparseMatrix & mtx)
{
	// std::cout << "AM HERE" << std::endl;
	os << std::endl;
	os << std::scientific;
	auto data = mtx.data();
	auto rowptr = mtx.row_ptr();

	for (auto i=0; i<mtx.rows(); i++){
		auto it = rowptr[i];
		while (it != rowptr[i+1] && it !=data.end()){
			os << "(" << i << "," << it->first << ") : " << it->second << std::endl;
			it++;
		}
	}
	

	return os;
}



// SparseMatrix-Matrix addition
Matrix SparseMatrix::operator+(const Matrix & mtx) const
{
	if (m_mrows != mtx.rows() || m_ncols != mtx.cols())
	{
		std::cout << "SparseMatrix-Matrix dimensions do not match!" << std::endl;
		throw -1;
	}

	Matrix out(mtx);

	for (auto i=0; i<mtx.rows(); i++){
		auto it = m_row_ptr[i];
		while (it != m_row_ptr[i+1] && it !=m_data.end()){
			out(i,it->first) += it->second;
			it++;
		}
	}

	return out;
}

// SparseMatrix-Matrix subtraction
Matrix SparseMatrix::operator-(const Matrix & mtx) const
{
	if (m_mrows != mtx.rows() || m_ncols != mtx.cols())
	{
		std::cout << "SparseMatrix-Matrix dimensions do not match!" << std::endl;
		throw -1;
	}

	Matrix out = -1*mtx;

	for (auto i=0; i<mtx.rows(); i++){
		auto it = m_row_ptr[i];
		while (it != m_row_ptr[i+1] && it !=m_data.end()){
			out(i,it->first) += it->second;
			it++;
		}
	}

	return out;
}


// Matrix-SparseMatrix addition
Matrix Matrix::operator+(const SparseMatrix & mtx) const
{
	if (m_mrows != mtx.rows() || m_ncols != mtx.cols())
	{
		std::cout << "SparseMatrix-Matrix dimensions do not match!" << std::endl;
		throw -1;
	}

	Matrix out(*this);

	auto data = mtx.data();
	auto rowptr = mtx.row_ptr();

	for (auto i=0; i<mtx.rows(); i++){
		auto it = rowptr[i];
		while (it != rowptr[i+1] && it !=data.end()){
			out(i,it->first) += it->second;
			it++;
		}
	}

	return out;
}


// Matrix-SparseMatrix subtraction
Matrix Matrix::operator-(const SparseMatrix & mtx) const
{
	if (m_mrows != mtx.rows() || m_ncols != mtx.cols())
	{
		std::cout << "SparseMatrix-Matrix dimensions do not match!" << std::endl;
		throw -1;
	}

	Matrix out(*this);

	auto data = mtx.data();
	auto rowptr = mtx.row_ptr();

	for (auto i=0; i<mtx.rows(); i++){
		auto it = rowptr[i];
		while (it != rowptr[i+1] && it !=data.end()){
			out(i,it->first) -= it->second;
			it++;
		}
	}

	return out;
}





// global operators - SparseMatrix
SparseMatrix operator*(double val, const SparseMatrix & mtx)
{
	SparseMatrix out(mtx);
	auto data = mtx.data();
	auto rowptr = mtx.row_ptr();
	for (auto i=0; i<mtx.rows(); i++){
		auto it = rowptr[i];
		while (it != rowptr[i+1] && it !=data.end()){
			out.set(i, it->first, it->second*val);
			it++;
		}
	}
	return out;
}

SparseMatrix operator+(double val, const SparseMatrix & mtx)
{
	SparseMatrix out(mtx);
	auto data = mtx.data();
	auto rowptr = mtx.row_ptr();
	for (auto i=0; i<mtx.rows(); i++){
		auto it = rowptr[i];
		while (it != rowptr[i+1] && it !=data.end()){
			out.set(i, it->first, it->second+val);
			it++;
		}
	}
	return out;
}

SparseMatrix operator-(double val, const SparseMatrix & mtx)
{
	SparseMatrix out(mtx);
	auto data = mtx.data();
	auto rowptr = mtx.row_ptr();
	for (auto i=0; i<mtx.rows(); i++){
		auto it = rowptr[i];
		while (it != rowptr[i+1] && it !=data.end()){
			out.set(i, it->first, val-it->second);
			it++;
		}
	}
	return out;
}

// identity matrix
SparseMatrix speye(unsigned int rows, unsigned int cols)
{
	SparseMatrix out(rows, cols);
	unsigned int sz = std::min(rows,cols);
	for (auto i=0; i<sz; i++)
	{
		out.set(i,i,1);
	}
	return out;
}

// sparse random matrix uniformly distributed [0,1]
// with a fill factor of fill
SparseMatrix sprandmat(unsigned int rows, unsigned int cols, double fill=0.2)
{
	// seed
	std::default_random_engine generator;
	std::uniform_real_distribution<double> distrib(0.0,1.0);

	SparseMatrix out(rows, cols);
	for (auto i=0; i<rows; i++){
		for (auto j=0; j<cols; j++){
			if (distrib(generator) < fill) out.set(i,j, distrib(generator));
		}
	}

	return out;
}

// sparse random matrix normally distributed
SparseMatrix sprandmatn(unsigned int rows, unsigned int cols, double fill=0.2)
{
	std::default_random_engine generator;
	std::normal_distribution<double> distrib(0.0,1.0);

	SparseMatrix out(rows, cols);
	for (auto i=0; i<rows; i++){
		for (auto j=0; j<cols; j++){
			if (distrib(generator) < fill) out.set(i,j, distrib(generator));
		}
	}

	return out;
}

#endif