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


template <typename VectorT>
void write_vector(const VectorT & v, std::ostream & os = std::cout, std::size_t ntabs = 0){
	for (auto i=0; i<ntabs; i++) os << "\t" ;
	os << "<Vector>" << std::endl;
	os << std::scientific;
	for (auto it = std::cbegin(v); it!=std::cend(v); it++){
		for (auto i=0; i<ntabs; i++) os << "\t" ;
		std::cout << *it << std::endl;
	}
	for (auto i=0; i<ntabs; i++) os << "\t" ;
	os << "</Vector>" << std::endl;
	return;
} 




// this requires that we can return iterators from the vector
// types using std::cbegin and std::cend
template <typename Vector1, typename Vector2>
decltype(auto) inner_product(const Vector1 & v1, const Vector2 & v2){

	// FIXME: This needs to be changed for when v1 is a complex number...
	// 		  it should get conjugated before taking the product
	auto it1=std::cbegin(v1);
	auto it2=std::cbegin(v2);
	auto out = (*it1)*(*it2);
	it1++; it2++;
	while (it1 != std::cend(v1) && it2 != std::cend(v2)){
		out += (*it1)*(*it2);
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
	auto it1 = std::begin(v1);
	typedef typename std::remove_reference<decltype(*it1)>::type 		NonRefT;
	for (auto it = std::begin(v1); it!=std::end(v1); it++){
		*it = static_cast<NonRefT>(val);
	}
	return;
} 


// random vector uniformly distributed [0,1]
template <typename VectorT>
void fill_rand(VectorT & v)
{
	auto it = v.begin();
	typedef typename std::remove_reference<decltype(*it)>::type 		NonRefT;
	
	unsigned int seed1 = std::chrono::system_clock::now().time_since_epoch().count();
	std::minstd_rand0 generator(seed1);
	std::uniform_real_distribution<NonRefT> distrib(0.0,1.0);

	for (auto it = std::begin(v); it!=std::end(v); it++){
		*it = distrib(generator);
	}
	return;
}

// random vector normally distributed
template <typename VectorT>
void fill_randn(VectorT & v)
{

	auto it = v.begin();
	typedef typename std::remove_reference<decltype(*it)>::type 		NonRefT;

	unsigned int seed1 = std::chrono::system_clock::now().time_since_epoch().count();
	std::minstd_rand0 generator(seed1);
	std::normal_distribution<NonRefT> distrib(0.0,1.0);

	for (auto it = std::begin(v); it!=std::end(v); it++){
		*it = distrib(generator);
	}
	return;
}



// this requires that we can iterate over the vector using std::cbegin and std::cend
template <typename VectorT>
std::size_t length(VectorT & v1){
	// std::size_t ct = 0;
	// for (auto it = std::cbegin(v1); it!=std::cend(v1); it++){
	// 	ct++;
	// }
	// return ct;
	return std::cend(v1) - std::cbegin(v1);
} 


template <typename VectorT>
std::size_t argmin(const VectorT & v){
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


template <typename VectorT>
decltype(auto) min(const VectorT & v){
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

template <typename VectorT>
std::size_t argmax(const VectorT & v){
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


template <typename VectorT>
decltype(auto) max(const VectorT & v){
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
	auto it = std::cbegin(v);
	auto val = std::abs(*it);
	it++;
	while (it != std::cend(v)){
		val = std::max(val, std::abs(*it));
		it++;
	}
	return val;
}

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