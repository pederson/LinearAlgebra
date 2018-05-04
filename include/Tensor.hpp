#ifndef _TENSOR_H
#define _TENSOR_H

#include <math.h>
#include <iostream>
#include <fstream>
#include <istream>
#include <sstream>
#include <string>
#include <stdlib.h>
#include <random>
#include <chrono>

namespace libra{

typedef std::size_t 			size_type;

// DO NOT CHANGE THIS...
// it defines the dynamic_size keyword that indicates that a dimension is not known at compile time
constexpr size_type dynamic_size = 0;;


// typedef for the basic storage type
template <typename value_type, size_type dim>
struct ArrayTypedef{
	typedef std::array<value_type, dim> type;
};
template <typename value_type>
struct ArrayTypedef<value_type, dynamic_size>{
	typedef std::vector<value_type> type;
};
template <typename value_type, size_type dim>
using ArrayType = typename ArrayTypedef<value_type, dim>::type;




// typedef for the tensor types
template <bool base, typename scalar_type, size_type dim1, size_type... dims_at_compile>
struct TableTypedef{
	typedef ArrayType<typename TableTypedef<base, scalar_type, dims_at_compile...>::type, dim1> type;
};
template <typename scalar_type, size_type dim1>
struct TableTypedef<true, scalar_type, dim1>{
	typedef ArrayType<scalar_type, dim1> type;
};
template <typename scalar_type, size_type rank, size_type... dims_at_compile>
using Table = typename TableTypedef<true, scalar_type, dims_at_compile...>::type;


template <typename scalar_type, size_type rank, size_type... dims_at_compile>
class Tensor : public Table<scalar_type, rank, dims_at_compile...>{
public:
	typedef Table<scalar_type, rank, dims_at_compile...> 	BaseType;

	// inherit the base class constructors
	using BaseType::BaseType;

	// this nasty-looking code simply allows the use of vector and array braces initializationn
	template <typename... Args>
    Tensor(Args &&... args) : BaseType({(std::forward<Args>(args))...}) {};


	// template <size_type... inds>
	// scalar_type & operator()(size_type i, inds... j){return (*this)[i].operator[](j...);}

	// template <size_type... inds>
	// scalar_type & operator()(size_type i, inds... i){return (*this)[i][j];}

};




} // end namespace libra




#endif