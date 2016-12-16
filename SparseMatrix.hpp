#ifndef _SPARSEMATRIX_H
#define _SPARSEMATRIX_H

#include <math.h>
#include <iostream>
#include <stdlib.h>
#include <list>
#include <vector>

#include "Matrix.hpp"


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
		, m_data 	(sm.m_data)
		, m_row_ptr (sm.m_row_ptr)
	{
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
		// sm1.m_row_ptr.swap(sm2.m_row_ptr);
		// sm1.m_data.swap(sm2.m_data);


		// sm1.m_row_ptr[sm1.m_mrows] = sm1.m_data.end();
		// sm2.m_row_ptr[sm2.m_mrows] = sm2.m_data.end();


	}


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


	// // SparseMatrix-SparseMatrix addition
	// SparseMatrix operator+(const SparseMatrix & vct) const
	// {
	// 	if (m_len != vct.m_len)
	// 	{
	// 		std::cout << "SparseMatrix dimensions do not match!" << std::endl;
	// 		throw -1;
	// 	}

	// 	SparseMatrix out(vct);

	// 	// first merge the two maps
	// 	out.m_data.insert(m_data.begin(), m_data.end());

	// 	// find all common elements of the two maps and adjust accordingly
	// 	for (auto it=m_data.begin(); it!=m_data.end(); it++)
	// 	{
	// 		auto fit = vct.m_data.find(it->first);
	// 		if (fit != vct.m_data.end()){
	// 			out.m_data[it->first] += m_data.at(it->first);
	// 		}
	// 	}

	// 	// then take care of the default value
	// 	out.m_default_value = m_default_value+vct.m_default_value;

	// 	return out;
	// }

	// // SparseMatrix-SparseMatrix subtraction
	// SparseMatrix operator-(const SparseMatrix & vct) const
	// {
	// 	if (m_len != vct.m_len)
	// 	{
	// 		std::cout << "SparseMatrix dimensions do not match!" << std::endl;
	// 		throw -1;
	// 	}

	// 	SparseMatrix out(vct*-1);

	// 	// first merge the two maps
	// 	out.m_data.insert(m_data.begin(), m_data.end());

	// 	// find all common elements of the two maps and adjust accordingly
	// 	for (auto it=m_data.begin(); it!=m_data.end(); it++)
	// 	{
	// 		auto fit = vct.m_data.find(it->first);
	// 		if (fit != vct.m_data.end()){
	// 			out.m_data[it->first] += m_data.at(it->first);
	// 		}
	// 	}

	// 	// then take care of the default value
	// 	out.m_default_value = m_default_value-vct.m_default_value;

	// 	return out;
	// }

	// // SparseMatrix-SparseMatrix addition shorthand
	// SparseMatrix & operator+=(const SparseMatrix & vct)
	// {
	// 	if (m_len != vct.m_len)
	// 	{
	// 		std::cout << "SparseMatrix dimensions do not match!" << std::endl;
	// 		throw -1;
	// 	}

	// 	// first merge the two maps
	// 	m_data.insert(vct.m_data.begin(), vct.m_data.end());

	// 	// find all common elements of the two maps and adjust accordingly
	// 	for (auto it=m_data.begin(); it!=m_data.end(); it++)
	// 	{
	// 		auto fit = vct.m_data.find(it->first);
	// 		if (fit != vct.m_data.end()){
	// 			m_data[it->first] += vct.m_data.at(it->first);
	// 		}
	// 	}

	// 	// then take care of the default value
	// 	m_default_value += vct.m_default_value;

	// 	return *this;
	// }

	// // SparseMatrix-SparseMatrix subtraction shorthand
	// SparseMatrix & operator-=(const SparseMatrix & vct)
	// {
	// 	if (m_len != vct.m_len)
	// 	{
	// 		std::cout << "SparseMatrix dimensions do not match!" << std::endl;
	// 		throw -1;
	// 	}

	// 	SparseMatrix nvct = vct*-1;

	// 	// first merge the two maps
	// 	m_data.insert(nvct.m_data.begin(), nvct.m_data.end());

	// 	// find all common elements of the two maps and adjust accordingly
	// 	for (auto it=m_data.begin(); it!=m_data.end(); it++)
	// 	{
	// 		auto fit = nvct.m_data.find(it->first);
	// 		if (fit != nvct.m_data.end()){
	// 			m_data[it->first] += nvct.m_data[it->first];
	// 		}
	// 	}

	// 	// then take care of the default value
	// 	m_default_value -= vct.m_default_value;

	// 	return *this;
	// }

	// // transpose
	// SparseMatrix operator~() const
	// {
	// 	SparseMatrix sv(*this);
	// 	sv.transpose();
	// 	return sv;
	// }


	// Vector operator+(const Vector & vct) const;


	// Vector operator-(const Vector & vct) const;


	// // scalar multiplication
	// SparseMatrix operator*(double val) const
	// {
	// 	SparseMatrix out(*this);
	// 	out.m_default_value *= val;
	// 	for (auto it=out.m_data.begin(); it!=out.m_data.end(); it++){
	// 		out.m_data[it->first]*=val;
	// 	}
	// 	return out;
	// }

	// // scalar division
	// SparseMatrix operator/(double val) const
	// {
	// 	SparseMatrix out(*this);
	// 	out.m_default_value /= val;
	// 	for (auto it=out.m_data.begin(); it!=out.m_data.end(); it++){
	// 		out.m_data[it->first]/=val;
	// 	}
	// 	return out;
	// }

	// // scalar addition
	// SparseMatrix operator+(double val) const
	// {
	// 	SparseMatrix out(*this);
	// 	out.m_default_value += val;
	// 	for (auto it=out.m_data.begin(); it!=out.m_data.end(); it++){
	// 		out.m_data[it->first]+=val;
	// 	}
	// 	return out;
	// }

	// // scalar subtraction
	// SparseMatrix operator-(double val) const
	// {
	// 	SparseMatrix out(*this);
	// 	out.m_default_value -= val;
	// 	for (auto it=out.m_data.begin(); it!=out.m_data.end(); it++){
	// 		out.m_data[it->first]-=val;
	// 	}
	// 	return out;
	// }

	// // scalar multiplication
	// SparseMatrix & operator*=(double val)
	// {
	// 	m_default_value *= val;
	// 	for (auto it=m_data.begin(); it!=m_data.end(); it++){
	// 		m_data[it->first]*=val;
	// 	}
	// 	return *this;
	// }

	// // scalar division
	// SparseMatrix & operator/=(double val)
	// {
	// 	m_default_value /= val;
	// 	for (auto it=m_data.begin(); it!=m_data.end(); it++){
	// 		m_data[it->first]/=val;
	// 	}
	// 	return *this;
	// }

	// // scalar addition
	// SparseMatrix & operator+=(double val)
	// {
	// 	m_default_value += val;
	// 	for (auto it=m_data.begin(); it!=m_data.end(); it++){
	// 		m_data[it->first]+=val;
	// 	}
	// 	return *this;
	// }

	// // scalar subtraction
	// SparseMatrix & operator-=(double val)
	// {
	// 	m_default_value -= val;
	// 	for (auto it=m_data.begin(); it!=m_data.end(); it++){
	// 		m_data[it->first]-=val;
	// 	}
	// 	return *this;
	// }

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

	// // dot product of two sparse vectors
	// static double dot(const SparseMatrix & v1, const SparseMatrix & v2)
	// {
	// 	if (v1.m_len != v2.m_len){
	// 		std::cout << "SparseMatrixs must be of same size to dot them!" << std::endl;
	// 		throw -1;
	// 	}

		
	// 	double def1 = v1.m_default_value;
	// 	double def2 = v2.m_default_value;

	// 	double result = v1.m_len*def1*def2;

	// 	// add the cross terms
	// 	if (def1 != 0){
	// 		for (auto it=v2.m_data.begin(); it!=v2.m_data.end(); it++){
	// 			result += def1*(it->second-def2);
	// 		}
	// 	}
	// 	if (def2 != 0){
	// 		for (auto it=v1.m_data.begin(); it!=v1.m_data.end(); it++){
	// 			result += def2*(it->second-def1);
	// 		}
	// 	}

	// 	auto il = v1.m_data.begin();
	// 	auto ir = v2.m_data.begin();
	// 	// add the intersection terms
	// 	while (il != v1.m_data.end() && ir != v2.m_data.end())
	//     {
	//         if (il->first < ir->first)
	//             ++il;
	//         else if (ir->first < il->first)
	//             ++ir;
	//         else
	//         {
	//             result += (il->second-def1)*(ir->second-def2);
	//             ++il;
	//             ++ir;
	//         }
	//     }

	// 	return result;
	// }

	// access data
	const std::list<std::pair<unsigned int, double>> & data() const{return m_data;};
	const std::vector<std::list<std::pair<unsigned int, double>>::iterator> & row_ptr() const{return m_row_ptr;};

	// number of "nonzeros"
	std::size_t nnz() const {return m_data.size();};

	// // return vector of "nonempty" indices
	// std::vector<unsigned int> support() const{
	// 	std::vector<unsigned int> v;
	// 	for(auto it = m_data.begin(); it != m_data.end(); it++) {
	// 	  v.push_back(it->first);
	// 	}
	// 	return v;
	// }

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
			if (it->first == j) it->second = val;
			else m_data.emplace(it, std::pair<unsigned int, double> (j,val));
		}

		// adjust the row iterator if necessary
		--it;
		if (ct==0){
			for (auto ro=minrow; ro<=i; ro++) m_row_ptr[ro] = it;
		}

	}

	// // add val to index i
	// void add(unsigned int i, double val){
	// 	int r = m_data.count(i);
	// 	if (r == 0) m_data[i] = m_default_value + val;
	// 	else m_data[i] += val;
	// }

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

	// std::cout << rowptr.size() << std::endl;
	// for (auto i=0; i<rowptr.size()-1; i++){
	// 	std::cout << "(" << i << ", " << rowptr[i]->first << "): " << rowptr[i]->second << std::endl;
	// }

	// std::cout << std::endl;
	// for (auto it=data.begin(); it!=data.end(); it++){
	// 	std::cout << "(" << it->first << "," << it->second << ") " << std::endl;
	// }
	// auto it=data.end();
	// std::cout << "(" << it->first << "," << it->second << ") " << std::endl;
	// it=data.begin();
	// std::cout << "(" << it->first << "," << it->second << ") " << std::endl;

	// for (auto i=0; i<mtx.rows(); i++){
	// 	for (auto it = rowptr[i]; (it !=rowptr[i+1] && it != data.end()); it++){
	// 		os << "(" << i << "," << it->first << ") : " << it->second << std::endl;
	// 	}
	// }

	
	for (auto i=0; i<mtx.rows(); i++){
		auto it = rowptr[i];
		while (it != rowptr[i+1] && it !=data.end()){
			os << "(" << i << "," << it->first << ") : " << it->second << std::endl;
			it++;
		}
	}
	

	return os;
}



// // SparseMatrix-Vector addition
// Vector SparseMatrix::operator+(const Vector & vct) const
// {
// 	if (m_len != vct.length())
// 	{
// 		std::cout << "SparseMatrix dimensions do not match!" << std::endl;
// 		throw -1;
// 	}

// 	Vector out(vct);

// 	// iterate through the map, ordered by key
// 	unsigned int i=0;
// 	for (auto it=m_data.begin(); it!=m_data.end(); it++)
// 	{
// 		while (i < it->first ){
// 			out(i) += m_default_value;
// 			i++;
// 		}
// 		out(i) += it->second;
// 		i++;
// 	}
// 	while (i < out.length()){
// 		out(i) += m_default_value;
// 		i++;
// 	}

// 	return out;
// }

// // SparseMatrix-Vector subtraction
// Vector SparseMatrix::operator-(const Vector & vct) const
// {
// 	if (m_len != vct.length())
// 	{
// 		std::cout << "SparseMatrix dimensions do not match!" << std::endl;
// 		throw -1;
// 	}

// 	Vector out(vct);

// 	// iterate through the map, ordered by key
// 	unsigned int i=0;
// 	for (auto it=m_data.begin(); it!=m_data.end(); it++)
// 	{
// 		while (i < it->first ){
// 			out(i) = m_default_value - out(i);
// 			i++;
// 		}
// 		out(i) = it->second - out(i);
// 		i++;
// 	}
// 	while (i < out.length()){
// 		out(i) = m_default_value - out(i);
// 		i++;
// 	}

// 	return out;
// }


// // Vector-SparseMatrix addition
// Vector Vector::operator+(const SparseMatrix & vct) const
// {
// 	if (m_len != vct.length())
// 	{
// 		std::cout << "SparseMatrix dimensions do not match!" << std::endl;
// 		throw -1;
// 	}

// 	Vector out(*this);
// 	const std::map<unsigned int, double> sdata = vct.data();

// 	// iterate through the map, ordered by key
// 	unsigned int i=0;
// 	double defval = vct.default_value();
// 	for (auto it=sdata.begin(); it!=sdata.end(); it++)
// 	{
// 		while (i < it->first ){
// 			out(i) += defval;
// 			i++;
// 		}
// 		out(i) += it->second;
// 		i++;
// 	}
// 	while (i < out.length()){
// 		out(i) += defval;
// 		i++;
// 	}

// 	return out;
// }


// // Vector-SparseMatrix subtraction
// Vector Vector::operator-(const SparseMatrix & vct) const
// {
// 	if (m_len != vct.length())
// 	{
// 		std::cout << "SparseMatrix dimensions do not match!" << std::endl;
// 		throw -1;
// 	}

// 	Vector out(*this);
// 	const std::map<unsigned int, double> sdata = vct.data();

// 	// iterate through the map, ordered by key
// 	unsigned int i=0;
// 	double defval = vct.default_value();
// 	for (auto it=sdata.begin(); it!=sdata.end(); it++)
// 	{
// 		while (i < it->first ){
// 			out(i) -= defval;
// 			i++;
// 		}
// 		out(i) -= it->second;
// 		i++;
// 	}
// 	while (i < out.length()){
// 		out(i) -= defval;
// 		i++;
// 	}

// 	return out;
// }





// // global operators - SparseMatrix
// SparseMatrix operator*(double val, const SparseMatrix & vct)
// {
// 	SparseMatrix out(vct.length(), val*vct.default_value());
// 	auto data = vct.data();
// 	for (auto it=data.begin(); it!=data.end(); it++){
// 		out.set(it->first, val*it->second);
// 	}
// 	return out;
// }

// SparseMatrix operator+(double val, const SparseMatrix & vct)
// {
// 	SparseMatrix out(vct.length(), val+vct.default_value());
// 	auto data = vct.data();
// 	for (auto it=data.begin(); it!=data.end(); it++){
// 		out.set(it->first, val+it->second);
// 	}
// 	return out;
// }

// SparseMatrix operator-(double val, const SparseMatrix & vct)
// {
// 	SparseMatrix out(vct.length(), val-vct.default_value());
// 	auto data = vct.data();
// 	for (auto it=data.begin(); it!=data.end(); it++){
// 		out.set(it->first, val-it->second);
// 	}
// 	return out;
// }

#endif