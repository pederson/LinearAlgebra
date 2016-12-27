#ifndef _SPARSEVECTOR_H
#define _SPARSEVECTOR_H

#include <math.h>
#include <iostream>
#include <stdlib.h>
#include <map>
#include <vector>

#include "Matrix.hpp"


// *******************************************************************************


class SparseVector : public AbstractVector
{
public:
	// constructor
	SparseVector()
		: m_len 	(0)
		, m_mrows 	(0)
		, m_ncols 	(0)
		, m_default_value	(0)
	{
	}

	// create and allocate a new sparse vector
	SparseVector(unsigned int length)
		: m_len 	(length)
		, m_mrows 	(length)
		, m_ncols 	(1)
		, m_default_value	(0)
	{
	}

	// create and allocate a new sparse vector
	SparseVector(unsigned int length, double default_value)
		: m_len 	(length)
		, m_mrows 	(length)
		, m_ncols 	(1)
		, m_default_value	(default_value)
	{
	}

	// create a sparse vector initialized by existing data
	SparseVector(unsigned int length, double default_value, 
				 unsigned int num_vals, const unsigned int * inds, const double * data)
		: m_len 	(length)
		, m_mrows 	(length)
		, m_ncols 	(1)
		, m_default_value	(default_value)
	{
		for (auto i=0; i<num_vals; i++){
			m_data[inds[i]] = data[i];
		}
	}

	// copy constructor
	SparseVector(const SparseVector & sv)
		: m_len 	(sv.m_len)
		, m_mrows 	(sv.m_mrows)
		, m_ncols 	(sv.m_ncols)
		, m_data 	(sv.m_data)
		, m_default_value	(sv.m_default_value)
	{
	}

	// destructor
	~SparseVector()
	{
	}

	friend void swap(SparseVector & sv1, SparseVector & sv2)
	{
		using std::swap;
		swap(sv1.m_len, sv2.m_len);
		swap(sv1.m_mrows, sv2.m_mrows);
		swap(sv1.m_ncols, sv2.m_ncols);
		swap(sv1.m_data, sv2.m_data);
		swap(sv1.m_default_value, sv2.m_default_value);

	}


	// SparseVector-SparseVector addition
	SparseVector operator+(const SparseVector & vct) const
	{
		if (m_len != vct.m_len)
		{
			std::cout << "SparseVector dimensions do not match!" << std::endl;
			throw -1;
		}

		SparseVector out(vct);

		// first merge the two maps
		out.m_data.insert(m_data.begin(), m_data.end());

		// find all common elements of the two maps and adjust accordingly
		for (auto it=m_data.begin(); it!=m_data.end(); it++)
		{
			auto fit = vct.m_data.find(it->first);
			if (fit != vct.m_data.end()){
				out.m_data[it->first] += m_data.at(it->first);
			}
		}

		// then take care of the default value
		out.m_default_value = m_default_value+vct.m_default_value;

		return out;
	}

	// SparseVector-SparseVector subtraction
	SparseVector operator-(const SparseVector & vct) const
	{
		if (m_len != vct.m_len)
		{
			std::cout << "SparseVector dimensions do not match!" << std::endl;
			throw -1;
		}

		SparseVector out(vct*-1);

		// first merge the two maps
		out.m_data.insert(m_data.begin(), m_data.end());

		// find all common elements of the two maps and adjust accordingly
		for (auto it=m_data.begin(); it!=m_data.end(); it++)
		{
			auto fit = vct.m_data.find(it->first);
			if (fit != vct.m_data.end()){
				out.m_data[it->first] += m_data.at(it->first);
			}
		}

		// then take care of the default value
		out.m_default_value = m_default_value-vct.m_default_value;

		return out;
	}

	// SparseVector-SparseVector addition shorthand
	SparseVector & operator+=(const SparseVector & vct)
	{
		if (m_len != vct.m_len)
		{
			std::cout << "SparseVector dimensions do not match!" << std::endl;
			throw -1;
		}

		// first merge the two maps
		m_data.insert(vct.m_data.begin(), vct.m_data.end());

		// find all common elements of the two maps and adjust accordingly
		for (auto it=m_data.begin(); it!=m_data.end(); it++)
		{
			auto fit = vct.m_data.find(it->first);
			if (fit != vct.m_data.end()){
				m_data[it->first] += vct.m_data.at(it->first);
			}
		}

		// then take care of the default value
		m_default_value += vct.m_default_value;

		return *this;
	}

	// SparseVector-SparseVector subtraction shorthand
	SparseVector & operator-=(const SparseVector & vct)
	{
		if (m_len != vct.m_len)
		{
			std::cout << "SparseVector dimensions do not match!" << std::endl;
			throw -1;
		}

		SparseVector nvct = vct*-1;

		// first merge the two maps
		m_data.insert(nvct.m_data.begin(), nvct.m_data.end());

		// find all common elements of the two maps and adjust accordingly
		for (auto it=m_data.begin(); it!=m_data.end(); it++)
		{
			auto fit = nvct.m_data.find(it->first);
			if (fit != nvct.m_data.end()){
				m_data[it->first] += nvct.m_data[it->first];
			}
		}

		// then take care of the default value
		m_default_value -= vct.m_default_value;

		return *this;
	}


	// SparseVector-SparseVector elementwise multiply
	SparseVector elem_mult(const SparseVector & vct) const{
		if (m_len != vct.m_len)
		{
			std::cout << "SparseVector dimensions do not match!" << std::endl;
			throw -1;
		}

		SparseVector out(vct);

		// first merge the two maps
		out.m_data.insert(m_data.begin(), m_data.end());

		// find all common elements of the two maps and adjust accordingly
		for (auto it=m_data.begin(); it!=m_data.end(); it++)
		{
			auto fit = vct.m_data.find(it->first);
			if (fit != vct.m_data.end()){
				out.m_data[it->first] *= m_data.at(it->first);
			}
		}

		// then take care of the default value
		out.m_default_value = m_default_value*vct.m_default_value;

		return out;
	}

	// SparseVector-SparseVector elementwise division
	SparseVector elem_div(const SparseVector & vct) const{
		if (m_len != vct.m_len)
		{
			std::cout << "SparseVector dimensions do not match!" << std::endl;
			throw -1;
		}

		SparseVector out(*this);
		auto data = vct.data();

		// first merge the two maps
		out.m_data.insert(data.begin(), data.end());

		// find all common elements of the two maps and adjust accordingly
		for (auto it=data.begin(); it!=data.end(); it++)
		{
			auto fit = vct.m_data.find(it->first);
			if (fit != vct.m_data.end()){
				out.m_data[it->first] /= data.at(it->first);
			}
		}

		// then take care of the default value
		out.m_default_value = m_default_value*vct.m_default_value;

		return out;
	}

	// transpose
	SparseVector operator~() const
	{
		SparseVector sv(*this);
		sv.transpose();
		return sv;
	}


	Vector operator+(const Vector & vct) const;


	Vector operator-(const Vector & vct) const;


	// scalar multiplication
	SparseVector operator*(double val) const
	{
		SparseVector out(*this);
		out.m_default_value *= val;
		for (auto it=out.m_data.begin(); it!=out.m_data.end(); it++){
			out.m_data[it->first]*=val;
		}
		return out;
	}

	// scalar division
	SparseVector operator/(double val) const
	{
		SparseVector out(*this);
		out.m_default_value /= val;
		for (auto it=out.m_data.begin(); it!=out.m_data.end(); it++){
			out.m_data[it->first]/=val;
		}
		return out;
	}

	// scalar addition
	SparseVector operator+(double val) const
	{
		SparseVector out(*this);
		out.m_default_value += val;
		for (auto it=out.m_data.begin(); it!=out.m_data.end(); it++){
			out.m_data[it->first]+=val;
		}
		return out;
	}

	// scalar subtraction
	SparseVector operator-(double val) const
	{
		SparseVector out(*this);
		out.m_default_value -= val;
		for (auto it=out.m_data.begin(); it!=out.m_data.end(); it++){
			out.m_data[it->first]-=val;
		}
		return out;
	}

	// scalar multiplication
	SparseVector & operator*=(double val)
	{
		m_default_value *= val;
		for (auto it=m_data.begin(); it!=m_data.end(); it++){
			m_data[it->first]*=val;
		}
		return *this;
	}

	// scalar division
	SparseVector & operator/=(double val)
	{
		m_default_value /= val;
		for (auto it=m_data.begin(); it!=m_data.end(); it++){
			m_data[it->first]/=val;
		}
		return *this;
	}

	// scalar addition
	SparseVector & operator+=(double val)
	{
		m_default_value += val;
		for (auto it=m_data.begin(); it!=m_data.end(); it++){
			m_data[it->first]+=val;
		}
		return *this;
	}

	// scalar subtraction
	SparseVector & operator-=(double val)
	{
		m_default_value -= val;
		for (auto it=m_data.begin(); it!=m_data.end(); it++){
			m_data[it->first]-=val;
		}
		return *this;
	}

	// value access
	// WARNING: will insert new element if it doesn't already exist
	double & operator()(unsigned int i)
	{
		return m_data[i];
	}

	// const value access
	double operator()(unsigned int i) const
	{
		int r = m_data.count(i);
		if (r==0) return m_default_value;
		return m_data.at(i);
	}

	// dot product of two sparse vectors
	static double dot(const SparseVector & v1, const SparseVector & v2)
	{
		if (v1.m_len != v2.m_len){
			std::cout << "SparseVectors must be of same size to dot them!" << std::endl;
			throw -1;
		}

		
		double def1 = v1.m_default_value;
		double def2 = v2.m_default_value;

		double result = v1.m_len*def1*def2;

		// add the cross terms
		if (def1 != 0){
			for (auto it=v2.m_data.begin(); it!=v2.m_data.end(); it++){
				result += def1*(it->second-def2);
			}
		}
		if (def2 != 0){
			for (auto it=v1.m_data.begin(); it!=v1.m_data.end(); it++){
				result += def2*(it->second-def1);
			}
		}

		auto il = v1.m_data.begin();
		auto ir = v2.m_data.begin();
		// add the intersection terms
		while (il != v1.m_data.end() && ir != v2.m_data.end())
	    {
	        if (il->first < ir->first)
	            ++il;
	        else if (ir->first < il->first)
	            ++ir;
	        else
	        {
	            result += (il->second-def1)*(ir->second-def2);
	            ++il;
	            ++ir;
	        }
	    }

		return result;
	}

	// access data
	const std::map<unsigned int, double> & data() const{return m_data;};

	// default value
	double default_value() const {return m_default_value;};

	// number of "nonempties"
	std::size_t support_size() const {return m_data.size();};

	// return vector of "nonempty" indices
	std::vector<unsigned int> support() const{
		std::vector<unsigned int> v;
		for(auto it = m_data.begin(); it != m_data.end(); it++) {
		  v.push_back(it->first);
		}
		return v;
	}

	// set index i to val
	void set(unsigned int i, double val){
		m_data[i] = val;
	}

	// add val to index i
	void add(unsigned int i, double val){
		int r = m_data.count(i);
		if (r == 0) m_data[i] = m_default_value + val;
		else m_data[i] += val;
	}

	// get value at index i
	double get(unsigned int i) const{
		int r = m_data.count(i);
		if (r==0) return m_default_value;
		return m_data.at(i);
	}

	void size(std::size_t & rows, std::size_t & cols) const{
		rows = m_mrows; cols = m_ncols;
	};

	std::size_t length() const{return m_len;};

	double norm() const{
		double val = double(m_len-m_data.size())*m_default_value;
		for (auto it=m_data.begin(); it!=m_data.end(); it++){
			val += it->second;
		}
		return val;
	};

	void transpose() {std::swap(m_mrows, m_ncols);};

	// return a dense vector version
	Vector densify() const{
		Vector out(m_len);
		out.fill(m_default_value);
		for (auto it=m_data.begin(); it!=m_data.end(); it++){
			out(it->first) = it->second;
		}
		return out;
	}


protected:

	std::size_t m_mrows, m_ncols;
	std::size_t m_len;
	std::map<unsigned int, double> m_data;
	double m_default_value;

};



std::ostream& operator<<(std::ostream & os, const SparseVector & vct)
{
	os << std::endl;
	os << std::scientific;
	os << "default = " << vct.default_value() << std::endl;
	auto data = vct.data();
	for (auto it=data.begin(); it!=data.end(); it++){
		os << "(" << it->first << ") : " << it->second << std::endl;
	}
	

	return os;
}



// SparseVector-Vector addition
Vector SparseVector::operator+(const Vector & vct) const
{
	if (m_len != vct.length())
	{
		std::cout << "SparseVector dimensions do not match!" << std::endl;
		throw -1;
	}

	Vector out(vct);

	// iterate through the map, ordered by key
	unsigned int i=0;
	for (auto it=m_data.begin(); it!=m_data.end(); it++)
	{
		while (i < it->first ){
			out(i) += m_default_value;
			i++;
		}
		out(i) += it->second;
		i++;
	}
	while (i < out.length()){
		out(i) += m_default_value;
		i++;
	}

	return out;
}

// SparseVector-Vector subtraction
Vector SparseVector::operator-(const Vector & vct) const
{
	if (m_len != vct.length())
	{
		std::cout << "SparseVector dimensions do not match!" << std::endl;
		throw -1;
	}

	Vector out(vct);

	// iterate through the map, ordered by key
	unsigned int i=0;
	for (auto it=m_data.begin(); it!=m_data.end(); it++)
	{
		while (i < it->first ){
			out(i) = m_default_value - out(i);
			i++;
		}
		out(i) = it->second - out(i);
		i++;
	}
	while (i < out.length()){
		out(i) = m_default_value - out(i);
		i++;
	}

	return out;
}


// Vector-SparseVector addition
Vector Vector::operator+(const SparseVector & vct) const
{
	if (m_len != vct.length())
	{
		std::cout << "SparseVector dimensions do not match!" << std::endl;
		throw -1;
	}

	Vector out(*this);
	const std::map<unsigned int, double> sdata = vct.data();

	// iterate through the map, ordered by key
	unsigned int i=0;
	double defval = vct.default_value();
	for (auto it=sdata.begin(); it!=sdata.end(); it++)
	{
		while (i < it->first ){
			out(i) += defval;
			i++;
		}
		out(i) += it->second;
		i++;
	}
	while (i < out.length()){
		out(i) += defval;
		i++;
	}

	return out;
}


// Vector-SparseVector subtraction
Vector Vector::operator-(const SparseVector & vct) const
{
	if (m_len != vct.length())
	{
		std::cout << "SparseVector dimensions do not match!" << std::endl;
		throw -1;
	}

	Vector out(*this);
	const std::map<unsigned int, double> sdata = vct.data();

	// iterate through the map, ordered by key
	unsigned int i=0;
	double defval = vct.default_value();
	for (auto it=sdata.begin(); it!=sdata.end(); it++)
	{
		while (i < it->first ){
			out(i) -= defval;
			i++;
		}
		out(i) -= it->second;
		i++;
	}
	while (i < out.length()){
		out(i) -= defval;
		i++;
	}

	return out;
}


// get support of sparse vector 
SparseVector Vector::get_support(const SparseVector & vct) const{
	std::vector<unsigned int> sup = vct.support();
	SparseVector out(m_len);
	for (auto it=sup.begin(); it!=sup.end(); it++){
		out(*it) = m_data[*it];
	}
	return out;
}




// global operators - SparseVector
SparseVector operator*(double val, const SparseVector & vct)
{
	SparseVector out(vct.length(), val*vct.default_value());
	auto data = vct.data();
	for (auto it=data.begin(); it!=data.end(); it++){
		out.set(it->first, val*it->second);
	}
	return out;
}

SparseVector operator+(double val, const SparseVector & vct)
{
	SparseVector out(vct.length(), val+vct.default_value());
	auto data = vct.data();
	for (auto it=data.begin(); it!=data.end(); it++){
		out.set(it->first, val+it->second);
	}
	return out;
}

SparseVector operator-(double val, const SparseVector & vct)
{
	SparseVector out(vct.length(), val-vct.default_value());
	auto data = vct.data();
	for (auto it=data.begin(); it!=data.end(); it++){
		out.set(it->first, val-it->second);
	}
	return out;
}

SparseVector sign(const SparseVector & vct){
	double deft = sgn(vct.default_value());
	SparseVector out(vct.length(), deft);
	auto dat = vct.data();
	for (auto it=dat.begin(); it!=dat.end(); it++){
		out(it->first) = sgn(it->second);
	}
	return out;
}

SparseVector abs(const SparseVector & vct){
	double deft = fabs(vct.default_value());
	SparseVector out(vct.length(), deft);
	auto dat = vct.data();
	for (auto it=dat.begin(); it!=dat.end(); it++){
		out(it->first) = fabs(it->second);
	}
	return out;
}

double norm_1(const SparseVector & vct){
	double res=0;
	auto dat = vct.data();
	for (auto it=dat.begin(); it!=dat.end(); it++){
		res += abs(it->second);
	}
	return res;
}

double norm_2(const SparseVector & vct){
	double res=0;
	auto dat = vct.data();
	for (auto it=dat.begin(); it!=dat.end(); it++){
		res += it->second*it->second;
	}
	return res;
}

double norm_inf(const SparseVector & vct){
	double res=0;
	auto dat = vct.data();
	for (auto it=dat.begin(); it!=dat.end(); it++){
		res = std::max(res, fabs(it->second));
	}
	return res;
}

// collection of matrix/vector generators
// sparse random vector uniformly distributed [0,1]
SparseVector sprandvec(unsigned int length, double fill)
{
	// seed
	// std::default_random_engine generator;
	unsigned int seed1 = std::chrono::system_clock::now().time_since_epoch().count();
	std::minstd_rand0 generator(seed1);
	std::uniform_real_distribution<double> distrib(0.0,1.0);

	SparseVector out(length);
	for (auto i=0; i<length; i++){
		if (distrib(generator) < fill) out(i) = distrib(generator);
	}
	return out;
}

// random vector normally distributed
SparseVector sprandvecn(unsigned int length, double fill)
{
	// seed
	// std::default_random_engine generator;
	unsigned int seed1 = std::chrono::system_clock::now().time_since_epoch().count();
	std::minstd_rand0 generator(seed1);
	std::normal_distribution<double> distrib(0.0,1.0);

	SparseVector out(length);
	for (auto i=0; i<length; i++){
		if (distrib(generator) < fill) out(i) = distrib(generator);
	}
	return out;
}

#endif