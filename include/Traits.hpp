#ifndef _TRAITS_H
#define _TRAITS_H


#include <stdlib.h>
#include <complex>
#include <type_traits>



// These are tools for detecting type traits

namespace libra{



namespace type_traits{

//***********************
	// This compile-time helper is used to detect whether a type is complex
	// and if it is, it gives the underlying type
	template <typename T>
	struct is_complex {
		static const bool value = false;
	};

	template <typename U>
	struct is_complex<std::complex<U>> {
		static const bool value = 		true;
		typedef U 			type;
	};








//***********************
	template<typename... Ts>
	struct is_vector_helper {};


//***********************
	// an assignable vector requires the non-const begin() and end()
	// iterator accessors
	template<typename T, typename _ = void>
	struct is_assignable_vector : std::false_type {};

	template<typename T>
	struct is_assignable_vector<
	        T,
	        std::conditional_t<
	            false,
	            is_vector_helper<
	                decltype(std::declval<T>().begin()),
	                decltype(std::declval<T>().end())
	                >,
	            void
	            >
	        > : public std::true_type {};


//***********************
	// a traversable vector requires the const cbegin() and cend()
	// iterator accessors
	template<typename T, typename _ = void>
	struct is_traversable_vector : std::false_type {};

	// this works because of SFINAE...
	// is_vector_helper<> is not ever compiled if the cbegin() and cend() functions
	// do not exist
	template<typename T>
	struct is_traversable_vector<
	        T,
	        std::conditional_t<
	            false,
	            is_vector_helper<
	                decltype(std::declval<T>().cbegin()),
	                decltype(std::declval<T>().cend())
	                >,
	            void
	            >
	        > : public std::true_type {};



//***********************
	// This compile-time helper is used to detect whether a type is a
	// vector type by checking if it is traversable or assignable (or
	// both)
	template<typename T>
	struct is_vector : public std::conditional_t<
												 is_assignable_vector<T>::value ||
												 is_traversable_vector<T>::value,
												 std::true_type,
												 std::false_type
												>::type {};



//***********************
	// a sizable vector is one which provides a size() method
	// to determine the length of the vector
	template<typename T, typename _ = void>
	struct is_sizable_vector : std::false_type {};

	template<typename T>
	struct is_sizable_vector<
	        T,
	        std::conditional_t<
	            false,
	            is_vector_helper<
	                decltype(std::declval<T>().size())
	                >,
	            void
	            >
	        > : public std::true_type {};


//***********************
	// a resizable vector is one which provides a resize() method
	// to change the length of the vector
	template<typename T, typename _ = void>
	struct is_resizable_vector : std::false_type {};

	template<typename T>
	struct is_resizable_vector<
	        T,
	        std::conditional_t<
	            false,
	            is_vector_helper<
	                decltype(std::declval<T>().resize(1))
	                >,
	            void
	            >
	        > : public std::true_type {};



//***********************
	// This compile-time helper is used to detect the underlying contained
	// type of a general container that only has a cbegin() function...
	// Note that we do not require that the vector has value_type defined
	template<typename T>
	struct contained_type{
		static_assert(is_vector<T>::value, "Must be a vector type in order to deduce contained type");
		typedef typename std::remove_const<typename std::remove_reference<decltype(*std::declval<T>().cbegin())>::type>::type type;
	};

	// template <typename P>
	// struct contained_type<P*>{
	// 	typedef P 		type;
	// };




//***********************
	// a random access iterator is one which has a valid
	// operator-
	template<typename T, typename _ = void>
	struct is_random_access_iterator : std::false_type {};

	template<typename T>
	struct is_random_access_iterator<
	        T,
	        std::conditional_t<
	            false,
	            is_vector_helper<
	                decltype(std::declval<T>()-std::declval<T>())
	                >,
	            void
	            >
	        > : public std::true_type {};










//***********************
	// the resulting type from a multiplicative product of two types
	// i.e. the type resulting from T1*T2
	template <typename T1, typename T2>
	struct product_type{
		typedef decltype(std::declval<T1>() * std::declval<T2>()) type;
	};

	template <typename T1, typename T2>
	struct sum_type{
		typedef decltype(std::declval<T1>() + std::declval<T2>()) type;
	};

	template <typename T1, typename T2>
	struct difference_type{
		typedef decltype(std::declval<T1>() - std::declval<T2>()) type;
	};



//***********************


} // end namespace type_traits


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