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



//////////////////////////////////////////////////
//////////////////////////////////////////////////




// fixed-length resize policy
template <size_type rows_at_compile, size_type cols_at_compile>
struct TensorResizePolicy{
	template <typename TensorType>
	static void resize(TensorType & v, size_type r, size_type c) {};
};
// dynamic columns size resize policy
template <size_type rows_at_compile>
struct TensorResizePolicy<rows_at_compile, dynamic_size>{
	template <typename TensorType>
	static void resize(TensorType & v, size_type r, size_type c) {
		for (auto i=0; i<r; i++)v[i].resize(c);};
};

// dynamic columns size resize policy
template <>
struct TensorResizePolicy<dynamic_size, dynamic_size>{
	template <typename TensorType>
	static void resize(TensorType & v, size_type r, size_type c) {
		v.resize(r);
		for (auto i=0; i<r; i++)v[i].resize(c);};
};


template <size_type... dims>
struct LinearMemoryLayout{
public:
	// void resize(
private:
};

namespace detail{

	// counts the number of args in the parameter pack
	template <size_type Arg1, size_type... Args>
	struct parameter_pack_size{
		static constexpr size_type value = 1+parameter_pack_size<Args...>::value;
	};
	template <size_type Arg1>
	struct parameter_pack_size<Arg1>{
		static constexpr size_type value = 1;
	};


	// counts the number of times "dynamic_size appears in the dims"
	template <size_type Arg1, size_type... Args>
	struct dynamic_size_count{
		static constexpr size_type value = (Arg1 == dynamic_size? 1 : 0)+dynamic_size_count<Args...>::value;
	};
	template <size_type Arg1>
	struct dynamic_size_count<Arg1>{
		static constexpr size_type value = (Arg1 == dynamic_size? 1 : 0);
	};

	
	// checks to make sure that the tensor dimensions are valid...
	// CANNOT have a non-dynamic size specified after a dynamic_size
	// i.e. Tensor<3, dynamic_size, 2, dynamic_size> ---> NOT OK
	// 		Tensor<3, 2, dynamic_size, dynamic_size> ---> OK
	template <bool is_dynamic, size_type Arg1>
	struct single_dimension_check{};
	template <size_type Arg1>
	struct single_dimension_check<false, Arg1>{
		static constexpr bool value = true;//(Arg1 == dynamic_size ? true : false);
	};
	template <size_type Arg1>
	struct single_dimension_check<true, Arg1>{
		static constexpr bool value = (Arg1 != dynamic_size ? false : true); 
	};
	template <bool is_dynamic, size_type Arg1, size_type... Args>
	struct dimension_check_impl{
		static constexpr bool value = dimension_check_impl<is_dynamic || Arg1 == dynamic_size, Args...>::value && single_dimension_check<is_dynamic, Arg1>::value;//(is_dynamic && Arg1 != dynamic_size) ? false : true;
	};
	template <bool is_dynamic, size_type Arg1>
	struct dimension_check_impl<is_dynamic, Arg1>{
		static constexpr bool value = single_dimension_check<is_dynamic, Arg1>::value;//(is_dynamic && Arg1 != dynamic_size) ? false : true;
	};
	template <size_type... Args>
	struct dimension_check{
		static constexpr bool value = dimension_check_impl<false, Args...>::value;
	};
}

template <typename scalar_type, typename memory_policy, typename storage_policy,
		  size_type... dims_at_compile
		  >
class Tensor{
public:
	static_assert(detail::dimension_check<dims_at_compile...>::value, "Tensor cannot have fixed size after a dynamic_size specification");
	typedef Tensor<scalar_type, memory_policy, storage_policy, dims_at_compile...> 		self_type; 
	static constexpr size_type tensor_rank = detail::parameter_pack_size<dims_at_compile...>::value;
	static constexpr size_type num_dynamic = detail::dynamic_size_count<dims_at_compile...>::value;

	// static constexpr size_type rank = detail::parameter_pack_size<dims_at_compile...>::value;
	constexpr size_type rank() const {return tensor_rank;};
	constexpr size_type ndynamic() const {return num_dynamic;};

	// typedef Table<scalar_type, rank, dims_at_compile...> 	BaseType;

	// inherit the base class constructors
	// using BaseType::BaseType;

	// this nasty-looking code simply allows the use of vector and array braces initializationn
	// template <typename... Args>
 //    Tensor(Args &&... args) : BaseType({(std::forward<Args>(args))...}) {};


	// template <size_type... inds>
	// scalar_type & operator()(size_type i, inds... j){return (*this)[i].operator[](j...);}

	// template <size_type... inds>
	// scalar_type & operator()(size_type i, inds... i){return (*this)[i][j];}

};




} // end namespace libra




#endif