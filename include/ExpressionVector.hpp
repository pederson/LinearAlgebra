#ifndef _EXPRESSIONVECTOR_H
#define _EXPRESSIONVECTOR_H


#include "Traits.hpp"


// // *************** CRTP class extender for unary vector operators

// 	template <typename Derived>
// 	struct UnaryVectorOperators{
// 	private:	
// 		Derived & derived() {return *static_cast<Derived *>(this);};
// 		const Derived & derived() const {return *static_cast<const Derived *>(this);}
// 	public:
// 		decltype(auto) cos() const {return libra::vector::cos(derived());};
// 		decltype(auto) sin() const {return libra::vector::sin(derived());};
// 		decltype(auto) tan() const {return libra::vector::tan(derived());};
// 		decltype(auto) acos() const {return libra::vector::acos(derived());};
// 		decltype(auto) asin() const {return libra::vector::asin(derived());};
// 		decltype(auto) atan() const {return libra::vector::atan(derived());};
// 		decltype(auto) cosh() const {return libra::vector::cosh(derived());};
// 		decltype(auto) sinh() const {return libra::vector::sinh(derived());};
// 		decltype(auto) tanh() const {return libra::vector::tanh(derived());};
// 		decltype(auto) acosh() const {return libra::vector::acosh(derived());};
// 		decltype(auto) asinh() const {return libra::vector::asinh(derived());};
// 		decltype(auto) atanh() const {return libra::vector::atanh(derived());};
// 		decltype(auto) exp() const {return libra::vector::exp(derived());};
// 		decltype(auto) log() const {return libra::vector::log(derived());};
// 		decltype(auto) log10() const {return libra::vector::log10(derived());};
// 		decltype(auto) exp2() const {return libra::vector::exp2(derived());};
// 		decltype(auto) expm1() const {return libra::vector::expm1(derived());};
// 		decltype(auto) ilogb() const {return libra::vector::ilogb(derived());};
// 		decltype(auto) log1p() const {return libra::vector::log1p(derived());};
// 		decltype(auto) log2() const {return libra::vector::log2(derived());};
// 		decltype(auto) logb() const {return libra::vector::logb(derived());};
// 		decltype(auto) sqrt() const {return libra::vector::sqrt(derived());};
// 		decltype(auto) cbrt() const {return libra::vector::cbrt(derived());};
// 		decltype(auto) hypot() const {return libra::vector::hypot(derived());};
// 		decltype(auto) erf() const {return libra::vector::erf(derived());};
// 		decltype(auto) erfc() const {return libra::vector::erfc(derived());};
// 		decltype(auto) tgamma() const {return libra::vector::tgamma(derived());};
// 		decltype(auto) lgamma() const {return libra::vector::lgamma(derived());};
// 		decltype(auto) ceil() const {return libra::vector::ceil(derived());};
// 		decltype(auto) floor() const {return libra::vector::floor(derived());};
// 		decltype(auto) trunc() const {return libra::vector::trunc(derived());};
// 		decltype(auto) round() const {return libra::vector::round(derived());};
// 		decltype(auto) abs() const {return libra::vector::abs(derived());};
// 	};



// } // end namespace vector


// This is a library for Expression Vectors... meaning that it
// represents a vector of expressions. 
// An expression is defined as a Type, which can access its data 
// through a static member function ::get()
// 
// Members are lazily evaluated through a ::get<I>(args...) function, 
// which accesses the I-th expression and evaluates it with 
// the arguments "args..."


namespace libra{



namespace Detail{

	template <typename ... Ts>
	struct tuple_reverse_impl{};

	template <typename T, size_t... I>
	struct tuple_reverse_impl<T, std::index_sequence<I...>>
	{
		typedef std::tuple<typename std::tuple_element<sizeof...(I) - 1 - I, T>::type...> type;
	};

	// partial specialization for handling empty tuples:
	template <typename T>
	struct tuple_reverse_impl<T, std::index_sequence<>>
	{
		typedef T type;
	};

	template <typename T>
	struct tuple_reverse
	: tuple_reverse_impl<T, std::make_index_sequence<std::tuple_size<T>::value>>
	{ };


	template <std::size_t I>
	struct IndexWrapper{
		typedef std::size_t type;
		static constexpr type value = I;
	};

	template <std::size_t I, std::size_t Iend>
	struct CreateTupleRange
	{
	private:
		template <typename T>
		struct CreateTuple{};

		template <std::size_t ... Idx>
		struct CreateTuple<std::index_sequence<Idx...>>{
		  typedef std::tuple<IndexWrapper<Idx+I>...> type;
		};

		typedef std::make_index_sequence<Iend+1-I>  IdxSeq;
	public:
		typedef typename CreateTuple<IdxSeq>::type  type;
	};

	template <std::size_t begin, std::size_t end>
	using RangeTuple = typename CreateTupleRange<begin, end>::type;

	template <std::size_t ... Idx>
	using IndexTuple = std::tuple<IndexWrapper<Idx>...>;



	template <typename T, typename... Ts>
	struct FirstType{
		typedef T type;
	};

	template <typename T>
	struct FirstType<T>{
		typedef T type;
	};


	// this does the implementation
	template<typename AggregateTuple, int I, template<typename...> class Fctor, typename TupleType1, typename... TupleTypes>
	struct nested_for_each_tuple_type_struct {
		template <typename... Args>
		static void get(Args && ... args) {
			// set the Ith type
			typedef std::tuple_element_t<I, TupleType1>   Type1;
			typedef std::tuple<Type1>                   SingletonTuple1;
			typedef decltype(std::tuple_cat(std::declval<AggregateTuple>(),std::declval<SingletonTuple1>()))    NewAggregate;

			// loop through inner loop
			nested_for_each_tuple_type_struct<NewAggregate, 
			                                std::tuple_size<typename FirstType<TupleTypes...>::type>::value-1, 
			                                Fctor,
			                                TupleTypes...>::get(std::forward<Args>(args)...);

			// on to the next iterate
			nested_for_each_tuple_type_struct<AggregateTuple, 
			                                I - 1, 
			                                Fctor, 
			                                TupleType1, 
			                                TupleTypes...>::get(std::forward<Args>(args)...);
		}
	};

	// base case
	template<typename AggregateTuple, template<typename...> class Fctor, typename TupleType1, typename... TupleTypes>
	struct nested_for_each_tuple_type_struct<AggregateTuple, 0, Fctor, TupleType1, TupleTypes...> {
		template <typename... Args>
		static void get(Args && ... args) {
			// set the Ith type
			typedef std::tuple_element_t<0, TupleType1>   Type1;
			typedef std::tuple<Type1>                   SingletonTuple1;
			typedef decltype(std::tuple_cat(std::declval<AggregateTuple>(),std::declval<SingletonTuple1>()))    NewAggregate;

			// loop through inner loop
			nested_for_each_tuple_type_struct<NewAggregate, 
			                                std::tuple_size<typename FirstType<TupleTypes...>::type>::value-1, 
			                                Fctor,
			                                TupleTypes...>::get(std::forward<Args>(args)...);
		}
	};

	// this does the implementation
	template<typename AggregateTuple, int I, template<typename...> class Fctor, typename TupleType>
	struct nested_for_each_tuple_type_struct<AggregateTuple, I, Fctor, TupleType>{
	private:
		template <typename Tuple>
		struct ApplyFunctor{
		};

		// specialize
		template <typename... TupleArgs>
		struct ApplyFunctor<std::tuple<TupleArgs...>>{
		  template <typename... Args>
		  static void get(Args && ... args){
		    Fctor<TupleArgs...>::get(std::forward<Args>(args)...);
		  }
		};
	public:
		template <typename... Args>
		static void get(Args && ... args) {

			// set the Ith type
			typedef std::tuple_element_t<I, TupleType>   Type1;
			typedef std::tuple<Type1>                   SingletonTuple1;
			typedef decltype(std::tuple_cat(std::declval<AggregateTuple>(),std::declval<SingletonTuple1>()))    NewAggregate;

			//Call get() of Fctor
			ApplyFunctor<NewAggregate>::get(std::forward<Args>(args)...);

			//Continue to next iteration
			nested_for_each_tuple_type_struct<AggregateTuple, I - 1, Fctor, TupleType>::get(std::forward<Args>(args)...);
		}
	};


	// this does the implementation
	template<typename AggregateTuple, template<typename...> class Fctor, typename TupleType>
	struct nested_for_each_tuple_type_struct<AggregateTuple, 0, Fctor, TupleType>{
	private:
		template <typename Tuple>
		struct ApplyFunctor{
		};

		// specialize
		template <typename... TupleArgs>
		struct ApplyFunctor<std::tuple<TupleArgs...>>{
			template <typename... Args>
			static void get(Args && ... args){
			Fctor<TupleArgs...>::get(std::forward<Args>(args)...);
			}
		};
	public:
		template <typename... Args>
		static void get(Args && ... args) {

			// set the Ith type
			typedef std::tuple_element_t<0, TupleType>   Type1;
			typedef std::tuple<Type1>                   SingletonTuple1;
			typedef decltype(std::tuple_cat(std::declval<AggregateTuple>(),std::declval<SingletonTuple1>()))    NewAggregate;

			//Call get() of Fctor
			ApplyFunctor<NewAggregate>::get(std::forward<Args>(args)...);
		}
	};



	// check if a tuple has a certain type T in its template parameters
	template <typename T, typename Tuple>
	struct has_type;

	template <typename T>
	struct has_type<T, std::tuple<>> : std::false_type {};

	template <typename T, typename U, typename... Ts>
	struct has_type<T, std::tuple<U, Ts...>> : has_type<T, std::tuple<Ts...>> {};

	template <typename T, typename... Ts>
	struct has_type<T, std::tuple<T, Ts...>> : std::true_type {};

	template <typename T, typename Tuple>
	using tuple_contains_type = typename has_type<T, Tuple>::type;
} // end namespace detail


// Fctor is a templated class with a single template parameter. 
// it is expected to contain a static function get() that accepts the same 
// number of arguments as are passed into this function
//
// Fctor<T>::get(Args...) is called for each type T in the tuple 
template<template<typename...> class Fctor, typename... TupleTypes, typename... Args>
void nested_for_each_tuple_type(Args && ... args) {
  Detail::nested_for_each_tuple_type_struct<std::tuple<>, std::tuple_size<typename Detail::FirstType<TupleTypes...>::type>::value-1, Fctor, typename Detail::tuple_reverse<TupleTypes>::type ...>::get(std::forward<Args>(args)...);
};


//************************************************************
//************************************************************


namespace Detail{
  template<std::size_t Sz, bool Col, typename... Ts>
  struct is_expression_vector_helper {};

} // end namespace Detail



//***********************
// an expression vector requires:
//	- static member function ::get<I>(...)
//	- static variable ::size
//	- static variable ::isColumn
//	- static typedef 	type<I>, returning the I-th element type
template<typename T, typename _ = void>
struct is_expression_vector : std::false_type {};

template <typename T>
struct is_expression_vector<
      T,
      std::conditional_t<
          false,
          Detail::is_expression_vector_helper<
          	  T::size,
              T::isColumn,
              typename T::template type<0>
              >,
          void
          >
      > : public std::true_type {};
// need to add the static ::get<I>(...) with a generic number of inputs


template <typename T, std::size_t ... Is>
struct expression_type{
};

// helper function
template <typename T, std::size_t I>
struct expression_type<T,I>{
	static_assert(is_expression_vector<T>::value, "Requires an expression vector!");
	typedef typename T::template type<I> 	type;
};

//************************************************************
//************************************************************



// plain old expression vector, directly stores elements
template <typename ... Elements>
struct ExpressionVector{
private:
	typedef std::tuple<Elements...>			TupleType;
public:
	static constexpr std::size_t 	size 		= std::tuple_size<TupleType>::value;
	static constexpr bool 			isColumn 	= true;

	template <std::size_t I>
	using type = typename std::tuple_element<I, TupleType>::type;

	template <std::size_t I, typename ... Args>
	static decltype(auto) get(Args && ... args){
		static_assert(I < size, "ExpressionVector ERROR: Index out of range");
		typedef typename std::tuple_element<I, TupleType>::type 	El;
		return El::get(std::forward<Args>(args)...);
	}
};

//************************************************************
//************************************************************


// Define an algebra for ExpressionVector
template <typename T1, typename T2>
struct ExpressionSum{
private:
	static_assert(is_expression_vector<T1>::value && is_expression_vector<T2>::value,
				  "ExpressionSum ERROR: input types must be expression vectors!");
	static_assert(T1::size == T2::size,
				  "ExpressionSum ERROR: input expression vectors must have same size!");
	static_assert(T1::isColumn == T2::isColumn,
				  "ExpressionSum ERROR: input expression vectors must have same orientation!");

	// define the atomic SumType so that individual elements of this class are ::get<I> -able
	template <typename I, typename J>
	struct SumType{

		template <typename ... Args>
		static decltype(auto) get(Args && ... args){
			return I::get(std::forward<Args>(args)...) + J::get(std::forward<Args>(args)...);
		};
	};
public:
	static constexpr std::size_t 	size 		= T1::size;
	static constexpr bool 			isColumn 	= T1::isColumn;

	template <std::size_t I>
	using type = SumType<typename expression_type<T1, I>::type, 
						 typename expression_type<T2, I>::type>;

	template <std::size_t I, typename ... Args>
	static decltype(auto) get(Args && ... args){
		static_assert(I < size, "ExpressionVector ERROR: Index out of range");
		typedef typename expression_type<T1, I>::type 	El1;
		typedef typename expression_type<T2, I>::type 	El2;

		return SumType<El1, El2>::get(std::forward<Args>(args)...);
	}
};

//************************************************************
//************************************************************

template <typename T1, typename T2>
struct ExpressionElementwiseProduct{
private:
	static_assert(is_expression_vector<T1>::value && is_expression_vector<T2>::value,
				  "ExpressionElementwiseProduct ERROR: input types must be expression vectors!");
	static_assert(T1::size == T2::size,
				  "ExpressionElementwiseProduct ERROR: input expression vectors must have same size!");
	// static_assert(T1::isColumn == T2::isColumn,
	// 			  "ExpressionElementwiseProduct ERROR: input expression vectors must have same orientation!");

	// define the atomic ProductType so that individual elements of this class are ::get<I> -able
	template <typename I, typename J>
	struct ProductType{
		template <typename ... Args>
		static decltype(auto) get(Args && ... args){
			return I::get(std::forward<Args>(args)...) * J::get(std::forward<Args>(args)...);
		};
	};
public:
	static constexpr std::size_t 	size 		= T1::size;
	static constexpr bool 			isColumn 	= T1::isColumn;

	template <std::size_t I>
	using type = ProductType<typename expression_type<T1, I>::type, 
						 typename expression_type<T2, I>::type>;

	template <std::size_t I, typename ... Args>
	static decltype(auto) get(Args && ... args){
		static_assert(I < size, "ExpressionElementwiseProduct ERROR: Index out of range");
		typedef typename expression_type<T1, I>::type 	El1;
		typedef typename expression_type<T2, I>::type 	El2;
		return ProductType<El1, El2>::get(std::forward<Args>(args)...);
	}
};


//************************************************************
//************************************************************

// template <typename T1, typename T2>
// struct ExpressionCrossProduct{
// private:
// 	static_assert(is_expression_vector<T1>::value && is_expression_vector<T2>::value,
// 				  "ExpressionCrossProduct ERROR: input types must be expression vectors!");
// 	static_assert(T1::size == T2::size,
// 				  "ExpressionCrossProduct ERROR: input expression vectors must have same size!");
// 	// static_assert(T1::isColumn == T2::isColumn,
// 	// 			  "ExpressionCrossProduct ERROR: input expression vectors must have same orientation!");

// 	// define the atomic ProductType so that individual elements of this class are ::get<I> -able
// 	template <typename I, typename J>
// 	struct ProductType{
// 		template <typename ... Args>
// 		static decltype(auto) get(Args && ... args){
// 			return I::get(std::forward<Args>(args)...) * J::get(std::forward<Args>(args)...);
// 		};
// 	};
// public:
// 	static constexpr std::size_t 	size 		= T1::size;
// 	static constexpr bool 			isColumn 	= T1::isColumn;

// 	template <std::size_t I>
// 	using type = ProductType<typename expression_type<I, T1>::type, 
// 						 typename expression_type<I, T2>::type>

// 	template <std::size_t I, typename ... Args>
// 	static decltype(auto) get(Args && ... args){
// 		static_assert(I < size, "ExpressionCrossProduct ERROR: Index out of range");
// 		typedef typename expression_type<I, T1>::type 	El1;
// 		typedef typename expression_type<I, T2>::type 	El2;
// 		return ProductType<El1, El2>::get(std::forward<Args>(args)...);
// 	}
// };


//************************************************************
//************************************************************



template <typename T>
struct ExpressionSumElements{
private:
	static_assert(is_expression_vector<T>::value,
				  "ExpressionSumElements ERROR: input type must be expression vector!");
	
	template <typename I>
	struct atomic_sum{
		template <typename ReturnType, typename ... Args>
		static void get(ReturnType && o, Args && ... args){
			o += T::template get<I::value>(std::forward<Args>(args)...);
		}
	};
public:
	static constexpr std::size_t 	size 		= 1;
	static constexpr bool 			isColumn 	= true;


	template <typename ... Args>
	static decltype(auto) get(Args && ... args){
		typedef std::decay_t<decltype(T::template get<0>(args...))> 	ReturnType;
		ReturnType out = T::template get<0>(std::forward<Args>(args)...);
		nested_for_each_tuple_type< atomic_sum, Detail::RangeTuple<1, T::size-1> >(out, std::forward<Args>(args)...);
		return out;
	}
};

//************************************************************
//************************************************************

template <typename T1, typename T2>
using ExpressionInnerProduct = ExpressionSumElements<ExpressionElementwiseProduct<T1, T2>>;


//************************************************************
//************************************************************

// template <typename T1, typename T2>
// struct ExpressionPair{
// private:
// 	static_assert(is_expression_vector<T1>::value && is_expression_vector<T2>::value,
// 				  "ExpressionPair ERROR: input types must be expression vectors!");
// 	static_assert(T1::size == T2::size,
// 				  "ExpressionPair ERROR: input expression vectors must have same size!");
// 	static_assert(T1::isColumn == T2::isColumn,
// 				  "ExpressionPair ERROR: input expression vectors must have same orientation!");

// 	// define the atomic SumType so that individual elements of this class are ::get<I> -able
// 	template <typename I, typename J>
// 	struct SumType{
// 		template <typename ... Args>
// 		static decltype(auto) get(Args && ... args){
// 			return I::get(std::forward<Args>(args)...) + J::get(std::forward<Args>(args)...);
// 		};
// 	};
// public:
// 	static constexpr std::size_t 	size 		= T1::size;
// 	static constexpr bool 			isColumn 	= T1::isColumn;

// 	template <std::size_t I>
// 	using type = SumType<typename expression_type<I, T1>::type, 
// 						 typename expression_type<I, T2>::type>

// 	template <std::size_t I, typename ... Args>
// 	static decltype(auto) get(Args && ... args){
// 		static_assert(I < size, "ExpressionVector ERROR: Index out of range");
// 		typedef typename expression_type<I, T1>::type 	El1;
// 		typedef typename expression_type<I, T2>::type 	El2;
// 		return SumType<El1, El2>::get(std::forward<Args>(args)...);
// 	}
// };






namespace Detail{
  template<std::size_t Sz, bool Col, typename... Ts>
  struct is_expression_matrix_helper {};

} // end namespace Detail



//***********************
// an expression vector requires:
//	- static member function ::get<I,J>(...)
//	- static member function ::row<I>(...)
//	- static member function ::col<J>(...)
//	- static variable ::rows
//	- static variable ::cols
//	- static typedef 	type<I, J>, returning the (I,J)-th element type
template<typename T, typename _ = void>
struct is_expression_matrix : std::false_type {};

template <typename T>
struct is_expression_matrix<
      T,
      std::conditional_t<
          false,
          Detail::is_expression_matrix_helper<
          	  T::rows,
          	  T::cols,
              typename T::template type<0,0>
              >,
          void
          >
      > : public std::true_type {};
// need to add the static ::get<I,J>(...) with a generic number of inputs


// helper function
template <typename T, std::size_t I, std::size_t J>
struct expression_type<T,I,J>{
	static_assert(is_expression_matrix<T>::value, "Requires an expression matrix!");
	typedef typename T::template type<I,J> 	type;
};


//************************************************************
//************************************************************



// plain old expression matrix, directly stores elements in row major order 
template <std::size_t R, std::size_t C, typename ... Elements>
struct ExpressionMatrix{
private:
	typedef std::tuple<Elements...>			TupleType;
	static constexpr std::size_t 	size 		= std::tuple_size<TupleType>::value;
	static_assert(R*C == size, "ExpressionMatrix ERROR: Number of elements must match Row x Col count");

public:
	
	static constexpr std::size_t 	rows 		= R;
	static constexpr std::size_t 	cols 		= C;

	template <std::size_t I, std::size_t J>
	using type = typename std::tuple_element<I*cols+J, TupleType>::type;

	template <std::size_t I, std::size_t J, typename ... Args>
	static decltype(auto) get(Args && ... args){
		static_assert(I < rows || J < cols, "ExpressionMatrix ERROR: Index out of range");
		typedef typename std::tuple_element<I*cols+J, TupleType>::type 	El;
		return El::get(std::forward<Args>(args)...);
	}
};


//************************************************************
//************************************************************



// matrix transpose
template <typename T>
struct ExpressionMatrixTranspose{
private:
	static_assert(is_expression_matrix<T>::value, "ExpressionMatrixTranspose ERROR: Input type must be expression matrix!");
	static_assert(T::rows == T::cols, "ExpressionMatrixSymmetricPart ERROR: Input matrix must be square!");
public:
	
	static constexpr std::size_t 	rows 		= T::cols;
	static constexpr std::size_t 	cols 		= T::rows;

	template <std::size_t I, std::size_t J>
	using type = typename T::template type<J,I>;

	template <std::size_t I, std::size_t J, typename ... Args>
	static decltype(auto) get(Args && ... args){
		static_assert(I < rows && J < cols, "ExpressionMatrixTranspose ERROR: Index out of range");
		typedef typename T::template type<J,I> 	El;
		return El::get(std::forward<Args>(args)...);
	}
};


//************************************************************
//************************************************************



// matrix symmetric part
template <typename T>
struct ExpressionMatrixSymmetricPart{
private:
	static_assert(is_expression_matrix<T>::value, "ExpressionMatrixSymmetricPart ERROR: Input type must be expression matrix!");
	static_assert(T::rows == T::cols, "ExpressionMatrixSymmetricPart ERROR: Input matrix must be square!");

	// define the atomic SymmType so that individual elements of this class are ::get<I> -able
	template <typename I, typename J>
	struct SymmType{

		template <typename ... Args>
		static decltype(auto) get(Args && ... args){
			return 0.5*(I::get(std::forward<Args>(args)...) + J::get(std::forward<Args>(args)...));
		};
	};
public:
	static constexpr std::size_t 	rows 		= T::rows;
	static constexpr std::size_t 	cols 		= T::cols;

	template <std::size_t I, std::size_t J>
	using type = SymmType<typename expression_type<T, I, J>::type, 
						 typename expression_type<T, J, I>::type>;

	template <std::size_t I, std::size_t J, typename ... Args>
	static decltype(auto) get(Args && ... args){
		static_assert(I < rows && J < cols, "ExpressionMatrixSymmetricPart ERROR: Index out of range");
		typedef typename expression_type<T, I, J>::type 	El1;
		typedef typename expression_type<T, J, I>::type 	El2;

		return SymmType<El1, El2>::get(std::forward<Args>(args)...);
	}
};

//************************************************************
//************************************************************



// matrix symmetric part
template <typename T>
struct ExpressionMatrixAntisymmetricPart{
private:
	static_assert(is_expression_matrix<T>::value, "ExpressionMatrixAntisymmetricPart ERROR: Input type must be expression matrix!");
	static_assert(T::rows == T::cols, "ExpressionMatrixAntisymmetricPart ERROR: Input matrix must be square!");

	// define the atomic AntisymmType so that individual elements of this class are ::get<I> -able
	template <typename I, typename J>
	struct AntisymmType{

		template <typename ... Args>
		static decltype(auto) get(Args && ... args){
			return 0.5*(I::get(std::forward<Args>(args)...) - J::get(std::forward<Args>(args)...));
		};
	};
public:
	static constexpr std::size_t 	rows 		= T::rows;
	static constexpr std::size_t 	cols 		= T::cols;

	template <std::size_t I, std::size_t J>
	using type = AntisymmType<typename expression_type<T, I, J>::type, 
						 typename expression_type<T, J, I>::type>;

	template <std::size_t I, std::size_t J, typename ... Args>
	static decltype(auto) get(Args && ... args){
		static_assert(I < rows && J < cols, "ExpressionMatrixAntisymmetricPart ERROR: Index out of range");
		typedef typename expression_type<T, I, J>::type 	El1;
		typedef typename expression_type<T, J, I>::type 	El2;

		return AntisymmType<El1, El2>::get(std::forward<Args>(args)...);
	}
};


//************************************************************
//************************************************************


// Define an algebra for ExpressionVector
template <typename T1, typename T2>
struct ExpressionMatrixSum{
private:
	static_assert(is_expression_matrix<T1>::value && is_expression_matrix<T2>::value,
				  "ExpressionMatrixSum ERROR: input types must be expression matrices!");
	static_assert(T1::rows == T2::rows && T1::cols == T2::cols,
				  "ExpressionMatrixSum ERROR: input expression matrices must have same size!");

	// define the atomic SumType so that individual elements of this class are ::get<I> -able
	template <typename I, typename J>
	struct SumType{
		template <typename ... Args>
		static decltype(auto) get(Args && ... args){
			return I::get(std::forward<Args>(args)...) + J::get(std::forward<Args>(args)...);
		};
	};
public:
	static constexpr std::size_t 	rows 		= T1::rows;
	static constexpr std::size_t 	cols 		= T1::cols;

	template <std::size_t I, std::size_t J>
	using type = SumType<typename expression_type<T1, I, J>::type, 
						 typename expression_type<T2, I, J>::type>;

	template <std::size_t I, std::size_t J, typename ... Args>
	static decltype(auto) get(Args && ... args){
		static_assert(I < rows && J < cols, "ExpressionMatrixSum ERROR: Index out of range");
		typedef typename expression_type<T1, I, J>::type 	El1;
		typedef typename expression_type<T2, I, J>::type 	El2;

		return SumType<El1, El2>::get(std::forward<Args>(args)...);
	}
};

//************************************************************
//************************************************************


// matrix row
template <std::size_t R, typename T>
struct ExpressionMatrixRow{
private:
	static_assert(R < T::rows, "ExpressionMatrixRow ERROR: Index out of range");
	static_assert(is_expression_matrix<T>::value, "ExpressionMatrixRow ERROR: Input type must be expression matrix!");
public:
	static constexpr std::size_t 	size 		= T::cols;
	static constexpr bool 			isColumn 	= false;

	template <std::size_t I>
	using type = typename T::template type<R, I>;

	template <std::size_t I, typename ... Args>
	static decltype(auto) get(Args && ... args){
		static_assert(I < T::cols, "ExpressionMatrixRow ERROR: Index out of range");
		typedef typename T::template type<R, I> 	El;
		return El::get(std::forward<Args>(args)...);
	}
};


//************************************************************
//************************************************************


// matrix row
template <std::size_t C, typename T>
struct ExpressionMatrixCol{
private:
	static_assert(C < T::cols, "ExpressionMatrixCol ERROR: Index out of range");
	static_assert(is_expression_matrix<T>::value, "ExpressionMatrixCol ERROR: Input type must be expression matrix!");
public:
	static constexpr std::size_t 	size 		= T::rows;
	static constexpr bool 			isColumn 	= true;

	template <std::size_t I>
	using type = typename T::template type<I, C>;

	template <std::size_t I, typename ... Args>
	static decltype(auto) get(Args && ... args){
		static_assert(I < T::rows, "ExpressionMatrixCol ERROR: Index out of range");
		typedef typename T::template type<I, C> 	El;
		return El::get(std::forward<Args>(args)...);
	}
};

//************************************************************
//************************************************************


// matrix diagonal
template <int Off, typename T>
struct ExpressionMatrixDiagonal{
private:
	static_assert((Off < 0 ? (-Off < T::rows) : (Off < T::cols)), "ExpressionMatrixDiagonal ERROR: Offset out of range");
	static_assert(is_expression_matrix<T>::value, "ExpressionMatrixDiagonal ERROR: Input type must be expression matrix!");
	static constexpr std::size_t 	rowStart 	= (Off < 0 ? -Off : 0);
	static constexpr std::size_t 	colStart 	= (Off >= 0 ? Off : 0);
public:
	static constexpr std::size_t 	size 		= (Off < 0 ? T::rows+Off : T::cols - Off);
	static constexpr bool 			isColumn 	= true;

	template <std::size_t I>
	using type = typename T::template type<rowStart+I, colStart+I>;

	template <std::size_t I, typename ... Args>
	static decltype(auto) get(Args && ... args){
		static_assert(I < size, "ExpressionMatrixDiagonal ERROR: Index out of range");
		typedef typename T::template type<rowStart+I, colStart+I> 	El;
		return El::get(std::forward<Args>(args)...);
	}
};


//************************************************************
//************************************************************


// matrix transpose
template <typename M, typename V>
struct ExpressionMatrixVectorProduct{
private:
	static_assert(is_expression_matrix<M>::value, "ExpressionMatrixVectorProduct ERROR: Input type must be expression matrix!");
	static_assert(is_expression_vector<V>::value, "ExpressionMatrixVectorProduct ERROR: Input type must be expression vector!");
	static_assert(M::cols == V::size, "ExpressionMatrixVectorProduct ERROR: Matrix-Vector dimensions must match!");

	template <std::size_t I>
	using Mrow = ExpressionMatrixRow<I, M>;
public:
	static constexpr std::size_t 	size 		= M::rows;

	template <std::size_t I>
	using type = ExpressionInnerProduct<Mrow<I>, V>;

	template <std::size_t I, typename ... Args>
	static decltype(auto) get(Args && ... args){
		static_assert(I < M::rows, "ExpressionMatrixVectorProduct ERROR: Index out of range");
		typedef ExpressionInnerProduct<Mrow<I>, V>	El;
		return El::get(std::forward<Args>(args)...);
	}
};


//************************************************************
//************************************************************


template <typename T>
struct ExpressionMatrixTrace{
private:
	static_assert(is_expression_matrix<T>::value,
				  "ExpressionMatrixTrace ERROR: input type must be expression matrix!");
	static_assert(T::rows == T::cols,
				  "ExpressionMatrixTrace ERROR: input matrix must be square!");
	
	template <typename I>
	struct atomic_sum{
		template <typename ReturnType, typename ... Args>
		static void get(ReturnType && o, Args && ... args){
			o += T::template get<I::value, I::value>(std::forward<Args>(args)...);
		}
	};
public:
	static constexpr std::size_t 	size 		= 1;
	static constexpr bool 			isColumn 	= true;


	template <typename ... Args>
	static decltype(auto) get(Args && ... args){
		typedef std::decay_t<decltype(T::template get<0,0>(args...))> 	ReturnType;
		ReturnType out = T::template get<0,0>(std::forward<Args>(args)...);
		nested_for_each_tuple_type< atomic_sum, Detail::RangeTuple<1, T::rows-1> >(out, std::forward<Args>(args)...);
		return out;
	}
};


//************************************************************
//************************************************************


template <typename T1, typename T2>
struct ExpressionMatrixConcatHorz{
private:
	static_assert(is_expression_matrix<T1>::value && is_expression_matrix<T2>::value,
				  "ExpressionMatrixConcatHorz ERROR: input types must be expression matrices!");
	static_assert(T1::rows == T2::rows,
				  "ExpressionMatrixConcatHorz ERROR: input expression matrices must have same number of rows!");

public:
	static constexpr std::size_t 	rows 		= T1::rows;
	static constexpr std::size_t 	cols 		= T1::cols+T2::cols;

	template <std::size_t I, std::size_t J>
	using type = std::conditional_t<I < T1::rows, 
				 typename expression_type<T1, I, J>::type,
				 typename expression_type<T2, I-T1::rows, J>::type>;

	template <std::size_t I, std::size_t J, typename ... Args>
	static decltype(auto) get(Args && ... args){
		static_assert(I < rows && J < cols, "ExpressionMatrixConcatHorz ERROR: Index out of range");
		typedef std::conditional_t<I < T1::rows, 
				 typename expression_type<T1, I, J>::type,
				 typename expression_type<T2, I-T1::rows, J>::type> 	El;

		return El::get(std::forward<Args>(args)...);
	}
};



//************************************************************
//************************************************************


template <typename T1, typename T2>
struct ExpressionMatrixConcatVert{
private:
	static_assert(is_expression_matrix<T1>::value && is_expression_matrix<T2>::value,
				  "ExpressionMatrixConcatVert ERROR: input types must be expression matrices!");
	static_assert(T1::cols == T2::cols,
				  "ExpressionMatrixConcatVert ERROR: input expression matrices must have same number of cols!");

	// define the atomic SumType so that individual elements of this class are ::get<I> -able
	template <typename I, typename J>
	struct SumType{
		template <typename ... Args>
		static decltype(auto) get(Args && ... args){
			return I::get(std::forward<Args>(args)...) + J::get(std::forward<Args>(args)...);
		};
	};
public:
	static constexpr std::size_t 	rows 		= T1::rows+T2::rows;
	static constexpr std::size_t 	cols 		= T1::cols;

	template <std::size_t I, std::size_t J>
	using type = std::conditional_t<J < T1::cols, 
				 typename expression_type<T1, I, J>::type,
				 typename expression_type<T2, I, J-T1::cols>::type>;

	template <std::size_t I, std::size_t J, typename ... Args>
	static decltype(auto) get(Args && ... args){
		static_assert(I < rows && J < cols, "ExpressionMatrixConcatVert ERROR: Index out of range");
		typedef std::conditional_t<J < T1::cols, 
				 typename expression_type<T1, I, J>::type,
				 typename expression_type<T2, I, J-T1::cols>::type> 	El;

		return El::get(std::forward<Args>(args)...);
	}
};


//************************************************************
//************************************************************


template <typename T1, typename T2>
struct ExpressionMatrixMatrixProduct{
private:
	static_assert(is_expression_matrix<T1>::value && is_expression_matrix<T2>::value,
				  "ExpressionMatrixMatrixProduct ERROR: input types must be expression matrices!");
	static_assert(T1::cols == T2::rows,
				  "ExpressionMatrixMatrixProduct ERROR: inner matrix dimensions must match!");

public:
	static constexpr std::size_t 	rows 		= T1::rows;
	static constexpr std::size_t 	cols 		= T2::cols;

	template <std::size_t I, std::size_t J>
	using type = ExpressionInnerProduct<ExpressionMatrixRow<I, T1>, 
										ExpressionMatrixCol<J, T2>>;

	template <std::size_t I, std::size_t J, typename ... Args>
	static decltype(auto) get(Args && ... args){
		static_assert(I < rows && J < cols, "ExpressionMatrixMatrixProduct ERROR: Index out of range");
		typedef ExpressionInnerProduct<ExpressionMatrixRow<I, T1>, 
										ExpressionMatrixCol<J, T2>> El;

		return El::get(std::forward<Args>(args)...);
	}
};


//************************************************************
//************************************************************

// Block matrices into a larger sparse matrix

// calculate the matrix determinant

// vector cross product

// matrix row sum -- outputs as vector

// matrix col sum -- outputs as vector

// matrix - matrix double contraction

// matrix and vector tensor product

// get a submatrix or subvector

// identity matrix

// zero matrix 

// unit vectors

// vector projection

// subtraction

//************************************************************
//************************************************************

// Compile-Time algorithms...
// 			- LU decomposition
//			- QR factorization

} // end namespace libra

///////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////

#endif