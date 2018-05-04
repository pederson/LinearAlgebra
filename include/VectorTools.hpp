#ifndef _VECTORTOOLS_H
#define _VECTORTOOLS_H

#include <stdio.h>
#include <stdlib.h>
#include <complex>
#include <algorithm>
#include <set>
#include <iterator>



// These are tools for operating on GENERAL VECTORS
// The requirements for an object to be a general vector are:
// 		- has iterator/const_iterator of random access type that can be 
//		  accessed through std::cbegin()/cend() and std::begin()/end()


namespace libra{


// This compile-time helper is used to detect whether a type is complex
// and if it is, it gives the underlying type
template <typename T>
struct complex_detector {
	static const bool value = false;
};

template <typename U>
struct complex_detector<std::complex<U>> {
	static const bool value = 		true;
	typedef U 			type;
};

template <typename T>
struct is_complex{
	static const bool value = complex_detector<T>::value ;
	// typedef typename complex_detector<T>::type 		type;// the underlying type of the complex number (double, float, int, etc...)
};


// This compile-time helper is used to detect whether a type is a
// vector type by checking if it has a cbegin() and cend() 
template<typename T, typename _ = void>
struct is_vector : std::false_type {};

template<typename... Ts>
struct is_vector_helper {};

template<typename T>
struct is_vector<
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
// template<typename T>
// struct is_vector<
//         T,
//         std::conditional_t<
//             false,
//             is_vector_helper<
//                 typename T::value_type,
//                 typename T::size_type,
//                 typename T::allocator_type,
//                 typename T::iterator,
//                 typename T::const_iterator,
//                 decltype(std::declval<T>().size()),
//                 decltype(std::declval<T>().begin()),
//                 decltype(std::declval<T>().end()),
//                 decltype(std::declval<T>().cbegin()),
//                 decltype(std::declval<T>().cend())
//                 >,
//             void
//             >
//         > : public std::true_type {};



// This compile-time helper is used to detect the underlying contained
// type of a general vector that only has a cbegin() and cend() function...
// Note that we do not require that the vector has value_type defined
template<typename T>
struct vector_detector{
	typedef typename std::remove_const<typename std::remove_reference<decltype(std::declval<decltype(std::declval<T>().cbegin())>().operator*())>::type>::type type;
};


// write a generalized vector to an output stream
template <typename VectorT>
void write_vector(const VectorT & v, std::ostream & os = std::cout, std::size_t ntabs = 0){
	static_assert(is_vector<VectorT>::value, "A Vector type requires a cbegin() and a cend() method!");

	for (auto i=0; i<ntabs; i++) os << "\t" ;
	os << "<Vector>" << std::endl;
	os << std::scientific;
	for (auto it = std::cbegin(v); it!=std::cend(v); it++){
		for (auto i=0; i<ntabs+1; i++) os << "\t" ;
		std::cout << *it << std::endl;
	}
	for (auto i=0; i<ntabs; i++) os << "\t" ;
	os << "</Vector>" << std::endl;
	return;
} 



// use these to switch between conjugation schemes... only used for complex-valued types
template <typename T>
struct inner_product_conjugation{
	static const T & get(const T & val) {return val;};
};

template <typename U>
struct inner_product_conjugation<std::complex<U>>{
	static std::complex<U> get(const std::complex<U> & val) {return std::conj(val);};
};

template <typename T1, typename T2>
struct product_type{
	typedef decltype(std::declval<T1>() * std::declval<T2>()) type;
};

// this requires that we can return iterators from the vector
// types using std::cbegin and std::cend
template <typename Vector1, typename Vector2>
decltype(auto) inner_product(const Vector1 & v1, const Vector2 & v2){
	static_assert(is_vector<Vector1>::value, "A Vector type requires a cbegin() and a cend() method!");
	static_assert(is_vector<Vector2>::value, "A Vector type requires a cbegin() and a cend() method!");

	// FIXME: the V2Type doesnt work for some reason...
	typedef typename vector_detector<Vector1>::type 		V1Type; // type of the first vector
	// typedef typename vector_detector<Vector2>::type 		V2Type; // type of the second vector

	auto it1=std::cbegin(v1);
	auto it2=std::cbegin(v2);
	typedef decltype((*it1)*(*it2))		ProductType;
	ProductType out = inner_product_conjugation<V1Type>::get(*it1)*(*it2);
	it1++; it2++;
	while (it1 != std::cend(v1) && it2 != std::cend(v2)){
		out += inner_product_conjugation<V1Type>::get(*it1)*(*it2);
		it1++;
		it2++;
	}

	return out;
} 



// this requires that we can iterate over the vector using std::begin and std::end
// and that the ValueT is at least implicitly convertible to 
// the VectorT::value_type
template <typename VectorT, typename ValueT>
void fill(VectorT & v1, ValueT val){
	static_assert(is_vector<VectorT>::value, "A Vector type requires a cbegin() and a cend() method!");

	auto it1 = std::begin(v1);
	typedef typename std::remove_reference<decltype(*it1)>::type 		NonRefT;
	for (auto it = std::begin(v1); it!=std::end(v1); it++){
		*it = static_cast<NonRefT>(val);
	}
	return;
} 


// this is used to switch the distribution type, esp. when it's 
// a complex-valued type
template <typename T>
struct distribution_type{
	typedef T 		type;
};

template <typename U>
struct distribution_type<std::complex<U>>{
	typedef U 		type;
};

template <typename T>
struct random_getter{
	template <typename distribution_t, typename generator_t>
	static T get(distribution_t & d, generator_t & g){return d(g);};
};

template <typename U>
struct random_getter<std::complex<U>>{
	template <typename distribution_t, typename generator_t>
	static std::complex<U> get(distribution_t & d, generator_t & g){return std::complex<U>(d(g),d(g));};
};

// random vector uniformly distributed [0,1]
template <typename VectorT>
void fill_rand(VectorT & v)
{
	static_assert(is_vector<VectorT>::value, "A Vector type requires a cbegin() and a cend() method!");

	auto it = v.begin();
	typedef typename std::remove_reference<decltype(*it)>::type 		NonRefT;
	typedef typename distribution_type<NonRefT>::type 					DistribT;

	unsigned int seed1 = std::chrono::system_clock::now().time_since_epoch().count();
	std::minstd_rand0 generator(seed1);
	std::uniform_real_distribution<DistribT> distrib(0.0,1.0);

	for (auto it = std::begin(v); it!=std::end(v); it++){
		*it = random_getter<NonRefT>::get(distrib, generator);//distrib(generator);
	}

	return;
}


// random vector normally distributed
template <typename VectorT>
void fill_randn(VectorT & v)
{
	static_assert(is_vector<VectorT>::value, "A Vector type requires a cbegin() and a cend() method!");

	auto it = v.begin();
	typedef typename std::remove_reference<decltype(*it)>::type 		NonRefT;
	typedef typename distribution_type<NonRefT>::type 					DistribT;

	// if (is_complex<NonRefT>::value) std::cout << "I am complex" << std::endl;
	unsigned int seed1 = std::chrono::system_clock::now().time_since_epoch().count();
	std::minstd_rand0 generator(seed1);
	std::normal_distribution<DistribT> distrib(0.0,1.0);

	for (auto it = std::begin(v); it!=std::end(v); it++){
		*it = random_getter<NonRefT>::get(distrib, generator);//distrib(generator);
	}
	return;
}



// this requires that we can iterate over the vector using std::cbegin and std::cend
template <typename VectorT>
std::size_t length(VectorT & v1){
	static_assert(is_vector<VectorT>::value, "A Vector type requires a cbegin() and a cend() method!");

	return std::cend(v1) - std::cbegin(v1);
} 

// index of the min value
template <typename VectorT>
std::size_t argmin(const VectorT & v){
	static_assert(is_vector<VectorT>::value, "A Vector type requires a cbegin() and a cend() method!");

	std::size_t res=0;
	auto mv = *v.begin();
	for (auto it=std::cbegin(v); it!=std::cend(v); it++){
		if (*it < mv){
			res = it - std::cbegin(v);
			mv = *it;
		}
	}
	return res;
}

// the min value
template <typename VectorT>
decltype(auto) min(const VectorT & v){
	static_assert(is_vector<VectorT>::value, "A Vector type requires a cbegin() and a cend() method!");

	std::size_t res=0;
	auto mv = *v.begin();
	for (auto it=std::cbegin(v); it!=std::cend(v); it++){
		if (*it < mv){
			res = it - std::cbegin(v);
			mv = *it;
		}
	}
	return mv;
}

// index of the max value
template <typename VectorT>
std::size_t argmax(const VectorT & v){
	static_assert(is_vector<VectorT>::value, "A Vector type requires a cbegin() and a cend() method!");

	std::size_t res=0;
	auto mv = *v.begin();
	for (auto it=std::cbegin(v); it!=std::cend(v); it++){
		if (*it > mv){
			res = it - std::cbegin(v);
			mv = *it;
		}
	}
	return res;
}

// the max value
template <typename VectorT>
decltype(auto) max(const VectorT & v){
	static_assert(is_vector<VectorT>::value, "A Vector type requires a cbegin() and a cend() method!");

	std::size_t res=0;
	auto mv = *v.begin();
	for (auto it=std::cbegin(v); it!=std::cend(v); it++){
		if (*it > mv){
			res = it - std::cbegin(v);
			mv = *it;
		}
	}
	return mv;
}




// A norm is specified as f(Sum(Rule(v_i))) where
// 			v_i 	= i-th element of the vector v
//			Rule 	= some function applied to each vector element (e.g. abs() or square)
//			Sum  	= the summation over the whole vector
//			f 		= some function applied to the sum (e.g. the square root)
template <typename NormRule, typename VectorT>
decltype(auto) norm(const VectorT & v){
	static_assert(is_vector<VectorT>::value, "A Vector type requires a cbegin() and a cend() method!");

	auto it = std::cbegin(v);
	auto val = NormRule::element_rule(*it);
	it++;
	while (it != std::cend(v)){
		val += NormRule::element_rule(*it);
		it++;
	}
	return NormRule::function_rule(val);
}

template <std::size_t p>
struct PNormRule {
	template <typename T>
	static T element_rule(T & val){
		return std::pow(val, p);
	};

	template <typename T>
	static T function_rule(T & val){
		return std::pow(val, 1.0/static_cast<double>(p));
	}
};


template <> 
struct PNormRule<1> {
	template <typename T>
	static T element_rule(T & val){
		return std::abs(val);
	};

	template <typename T>
	static T function_rule(T & val){
		return val;
	}
};

// some common norm functions defined for convenience
template <typename VectorT>
decltype(auto) norm_inf(const VectorT & v){
	static_assert(is_vector<VectorT>::value, "A Vector type requires a cbegin() and a cend() method!");

	auto it = std::cbegin(v);
	auto val = std::abs(*it);
	it++;
	while (it != std::cend(v)){
		val = std::max(val, std::abs(*it));
		it++;
	}
	return val;
}

// some convenient declarations 
template <typename VectorT>
decltype(auto) norm_1(const VectorT & v){
	return norm<PNormRule<1>>(v);
}

template <typename VectorT>
decltype(auto) norm_2(const VectorT & v){
	return norm<PNormRule<2>>(v);
}

template <typename VectorT>
decltype(auto) norm_3(const VectorT & v){
	return norm<PNormRule<3>>(v);
}

template <typename VectorT>
decltype(auto) norm_4(const VectorT & v){
	return norm<PNormRule<4>>(v);
}

template <typename VectorT>
decltype(auto) norm_5(const VectorT & v){
	return norm<PNormRule<5>>(v);
}


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