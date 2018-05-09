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



	// check if the parameter pack is all of one type T
	template <typename T, typename Arg1, typename... Args>
	struct is_one_type{
		static constexpr bool value = std::is_same<Arg1, T>::value && is_one_type<T, Args...>::value;
	};
	template <typename T, typename Arg1>
	struct is_one_type<T, Arg1>{
		static constexpr bool value = std::is_same<Arg1, T>::value;
	};



	// convert dims parameter pack to an array
	template <size_type... dims>
	struct dims_to_array{
		static constexpr std::array<size_type, parameter_pack_size<dims...>::value> value{dims...}; 
	};


	// cumulative product of an array 
	// (e.g. {a[2]*a[1]*a[0], a[1]*a[0], a[0], 1} = cumulative_product(a))
	template<size_type N>
	std::array<size_type, N> cumulative_product(const std::array<size_type, N> & a){
		std::array<size_type, N> out = a;
		for (auto it = out.begin(); it!=out.end(); it++) *it = 1;
		for (auto i=1; i<out.size(); i++){
			out[i] = out[i-1]*a[i-1];
		}	
		return out;
	}


	// template <size_type N>
	// struct 
}

template <typename scalar_type, size_type... dims_at_compile>
struct LinearMemory{
protected:
	static constexpr size_type tensor_rank = detail::parameter_pack_size<dims_at_compile...>::value;
	static constexpr size_type num_dynamic = detail::dynamic_size_count<dims_at_compile...>::value;
	std::array<size_type, tensor_rank> 		mDims = detail::dims_to_array<dims_at_compile...>::value;
	std::array<size_type, tensor_rank>		mCumProd = detail::cumulative_product(mDims);

	std::vector<scalar_type> mValues;		// this is the data


	void resize(std::array<size_type, tensor_rank> & d) {
		size_type l = 1;
		for (auto it=d.begin(); it!=d.end(); it++) l *= *it;
		mValues.resize(l);
		mCumProd = detail::cumulative_product(mDims);
	};

public:

	constexpr std::array<size_type, tensor_rank> dims() const {return mDims;};
	
	// resize when num_dynamic size_types are passes.... only resizes the dynamic portions
	template <typename... Args>
	void resize(size_type d1, Args... d){
		static_assert(sizeof...(Args) == num_dynamic-1, "Tensor constructor and resize may only fully specify dynamic_size dimensions!");
		// auto mydims = detail::dims_to_array<dims...>::value;
		std::array<size_type, num_dynamic> dnew{{d1, size_type(d)...}};
		for (auto i=tensor_rank - num_dynamic; i<tensor_rank; i++) mDims[i] = dnew[i-(tensor_rank - num_dynamic)];
		resize(mDims);
	};

	// element access operator by index pack
	// must fully specify all the indices in order to use this
	template <typename... Args>
	scalar_type & operator()(size_type d1, Args... d){
		static_assert(sizeof...(Args) == tensor_rank-1, "Tensor index accessor may only fully specify indices!");
		// would prefer not to make this intermediate array
		std::array<size_type, 1+sizeof...(Args)> indpack = {d1, size_type(d)...};
		size_type ind = 0;
		for (auto i=0; i<sizeof...(Args)+1; i++) ind += mCumProd[i]*indpack[i];
		//= inner_product(mCumProd, indpack);
		return mValues[ind];
	};

	// const element access operator by index pack
	template <typename... Args>
	const scalar_type & operator()(size_type d1, Args... d) const {
		return operator()(d1, d...);
	};

	
};



template <typename scalar_type, template <typename, size_type...> typename memory_policy,
		  size_type... dims_at_compile
		  >
class Tensor : public memory_policy<scalar_type, dims_at_compile...>{
public:
	static_assert(detail::dimension_check<dims_at_compile...>::value, "Tensor cannot have fixed size after a dynamic_size specification");
	
	typedef memory_policy<scalar_type, dims_at_compile...> 		memory;
	// typedef Tensor<scalar_type, memory_policy, dims_at_compile...> 		self_type; 
	static constexpr size_type tensor_rank = detail::parameter_pack_size<dims_at_compile...>::value;
	static constexpr size_type num_dynamic = detail::dynamic_size_count<dims_at_compile...>::value;


	// static constexpr size_type rank = detail::parameter_pack_size<dims_at_compile...>::value;
	constexpr size_type rank() const {return tensor_rank;};
	constexpr size_type ndynamic() const {return num_dynamic;};
	

	// empty constructor
	Tensor() {};

	// expose a constructor with num_dynamic arguments of type size_type
	template <typename... Args>
	Tensor(size_type d1, Args... dims){
		memory::resize(d1, dims...);
	}

	// expose a constructor using initializer lists


	

	// MEMORY OPERATIONS
		// expose a resize() function
		// expose an operator[] that returns a lower-rank tensor view
		// expose an operator(i,j,k) that returns an element or a tensor view
		// expose iterator and const_iterator


protected:
	// std::array<size_type, tensor_rank> 		mDims = detail::dims_to_array<dims_at_compile...>::value;
};




} // end namespace libra




#endif