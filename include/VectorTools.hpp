#ifndef _VECTORTOOLS_H
#define _VECTORTOOLS_H

#include <stdio.h>
#include <stdlib.h>
#include <complex>
#include <algorithm>
#include <set>
#include <iterator>
#include <iomanip>
#include <cstring>


#include "Traits.hpp"



// These are tools for operating on GENERAL VECTORS
// The requirements for an object to be a general vector are:
// 		- has iterator/const_iterator of random access type that can be 
//		  accessed through std::cbegin()/cend() and/or std::begin()/end
// 		  or by the equivalent member functions
// 		- can determine the size either through a size() function
//		  or by subtracting random access iterators

namespace libra{

namespace vector{

//***********************
	// write a generalized vector to an output stream
	template <bool Header = false, bool AsColumn = false, typename VectorT>
	void write(const VectorT & v, std::string dlm = " ", std::ostream & os = std::cout, std::size_t ntabs = 0){
		static_assert(type_traits::is_vector<VectorT>::value, "A Vector type requires a cbegin() and a cend() method!");

		for (auto i=0; i<ntabs; i++) os << "\t" ;
		if (Header) os << "<Vector>";
		if (Header && AsColumn) os << std::endl;
		os << std::scientific << std::setprecision(10);

		// first item
		auto it = v.cbegin();
		if (AsColumn) for (auto i=0; i<ntabs+1; i++) os << "\t" ;
		os << *it ;
		if (AsColumn) os << std::endl;
		it++;

		// remaining items get a delimiter
		while (it!= v.cend()){
			if (!AsColumn) os << dlm;
			if (AsColumn) for (auto i=0; i<ntabs+1; i++) os << "\t" ;
			os << *it ;
			if (AsColumn) os << std::endl;
			it++;
		}
		for (auto i=0; i<ntabs; i++) os << "\t" ;
		if (Header) os << "</Vector>" << std::endl;
		return;
	} 


	// the read from a stream into a general vector
	// returns the length of the read
	template <bool count_only=false, typename VectorT>
	std::size_t read_to(std::istream & is, VectorT && v){
		static_assert(type_traits::is_vector<VectorT>::value, "A Vector type requires a begin() and a end() method!");

		// underlying data type of VectorT
		typedef std::remove_reference_t<decltype(*v.begin())> T;
		// determine if type T has a defined stream operator>>


		// save the beginning of the stream
		std::streampos beg0 = is.tellg();
		std::streampos isbegin = is.tellg();

		// determine if there is a <Vector> header
		char ch;
		std::string head;
		is.get(ch);
		if (ch=='<'){
			is >> head;
			if (!strcmp(head.c_str(),"Vector>")){
				// go past the carrot
				isbegin = is.tellg();
			}
		}
		is.seekg(isbegin);


		// determine delimiter
		char dlm;
		T val;
		is >> val;
		// std::cout << "first value: " << val << std::endl;
		is.get(ch); 
		dlm = ch;
		is.seekg(isbegin);


		// std::cout << "delimiter is determined to be: " << dlm << std::endl;


		// start going through the stream
		std::size_t ct=0;
		std::string line;
		std::stringstream ss;
		if (count_only){
			while (getline(is, line, dlm)){
				ss << line;
				ss >> val;
				ct++;
			}
		}
		else{
			auto it=v.begin();
			while (getline(is, line, dlm)){
				ss << line;
				ss >> *it;
				// proceed
				it++;
				ct++;
			}
		}

		// std::cout << "looped through " << ct << " values" << std::endl;

		// if (is.good()) std::cout << "istream is still GOOD" << std::endl;
		// else std::cout << "istream is BAD" << std::endl;
		is.clear();
		is.seekg(beg0);
		// if (is.good()) std::cout << "istream is still GOOD" << std::endl;
		// else std::cout << "istream is BAD" << std::endl;
		// return the count
		return ct;
	}




//***********************
	// use these to switch between conjugation schemes... only used for complex-valued types
	template <typename T>
	struct inner_product_conjugation{
		static const T & get(const T & val) {return val;};
	};

	template <typename U>
	struct inner_product_conjugation<std::complex<U>>{
		static std::complex<U> get(const std::complex<U> & val) {return val;};//std::conj(val);};
	};

	// this requires that we can return iterators from the vector
	// types using std::cbegin and std::cend
	template <typename Vector1, typename Vector2>
	typename type_traits::product_type<
				typename type_traits::contained_type<Vector1>::type, 
				typename type_traits::contained_type<Vector2>::type
				>::type
	inner_product(const Vector1 & v1, const Vector2 & v2){
		static_assert(type_traits::is_vector<Vector1>::value, "A Vector type requires a cbegin() and a cend() method!");
		static_assert(type_traits::is_vector<Vector2>::value, "A Vector type requires a cbegin() and a cend() method!");

		// FIXME: the V2Type doesnt work for some reason...
		typedef typename type_traits::contained_type<Vector1>::type 		V1Type; // type of the first vector
		// typedef typename type_traits::contained_type<Vector2>::type 		V2Type; // type of the second vector

		auto it1=v1.cbegin();
		auto it2=v2.cbegin();
		typedef decltype((*it1)*(*it2))		ProductType;
		ProductType out = inner_product_conjugation<V1Type>::get(*it1)*(*it2);
		it1++; it2++;
		while (it1 != v1.cend() && it2 != v2.cend()){
			out += inner_product_conjugation<V1Type>::get(*it1)*(*it2);
			it1++;
			it2++;
		}

		return out;
	} 




//***********************
	// this requires that we can iterate over the vector using std::begin and std::end
	// and that the ValueT is at least implicitly convertible to 
	// the VectorT::value_type
	template <typename VectorT, typename ValueT>
	void fill(VectorT && v1, ValueT val){
		static_assert(type_traits::is_vector<VectorT>::value, "A Vector type requires a cbegin() and a cend() method!");

		auto it1 = v1.begin();
		typedef typename std::remove_reference<decltype(*it1)>::type 		NonRefT;
		for (auto it = v1.begin(); it!=v1.end(); it++){
			*it = static_cast<NonRefT>(val);
		}
		return;
	} 


//***********************
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
		static_assert(type_traits::is_vector<VectorT>::value, "A Vector type requires a cbegin() and a cend() method!");

		auto it = v.begin();
		typedef typename std::remove_reference<decltype(*it)>::type 		NonRefT;
		typedef typename distribution_type<NonRefT>::type 					DistribT;

		unsigned int seed1 = std::chrono::system_clock::now().time_since_epoch().count();
		std::minstd_rand0 generator(seed1);
		std::uniform_real_distribution<DistribT> distrib(0.0,1.0);

		for (auto it = v.begin(); it!=v.end(); it++){
			*it = random_getter<NonRefT>::get(distrib, generator);//distrib(generator);
		}

		return;
	}


	// random vector normally distributed
	template <typename VectorT>
	void fill_randn(VectorT & v)
	{
		static_assert(type_traits::is_vector<VectorT>::value, "A Vector type requires a cbegin() and a cend() method!");

		auto it = v.begin();
		typedef typename std::remove_reference<decltype(*it)>::type 		NonRefT;
		typedef typename distribution_type<NonRefT>::type 					DistribT;

		unsigned int seed1 = std::chrono::system_clock::now().time_since_epoch().count();
		std::minstd_rand0 generator(seed1);
		std::normal_distribution<DistribT> distrib(0.0,1.0);

		for (auto it = v.begin(); it!=v.end(); it++){
			*it = random_getter<NonRefT>::get(distrib, generator);//distrib(generator);
		}
		return;
	}



//***********************
	// // this requires that we can iterate over the vector using std::cbegin and std::cend
	// template <typename VectorT, typename T = size_type>
	// typename std::enable_if<type_traits::is_random_access_iterator<decltype(std::declval<VectorT>().cbegin())>::value, T>::type
	// length(const VectorT & v1){
	// 	static_assert(type_traits::is_vector<VectorT>::value, "A Vector type requires a cbegin() and a cend() method!");

	// 	return v1.cend() - v1.cbegin();
	// } 

	// this requires that we can iterate over the vector using std::cbegin and std::cend
	template <typename VectorT, typename T = size_type>
	typename std::enable_if<type_traits::is_sizable_vector<VectorT>::value, T>::type
	length(const VectorT & v1){
		static_assert(type_traits::is_vector<VectorT>::value, "A Vector type requires a cbegin() and a cend() method!");

		return v1.size();
	} 


//***********************
	// index of the min value
	template <typename VectorT>
	std::size_t argmin(const VectorT & v){
		static_assert(type_traits::is_vector<VectorT>::value, "A Vector type requires a cbegin() and a cend() method!");

		std::size_t res=0;
		auto mv = *v.begin();
		for (auto it=v.cbegin(); it!=v.cend(); it++){
			if (*it < mv){
				res = it - v.cbegin();
				mv = *it;
			}
		}
		return res;
	}

	// the min value
	template <typename VectorT>
	decltype(auto) min(const VectorT & v){
		static_assert(type_traits::is_vector<VectorT>::value, "A Vector type requires a cbegin() and a cend() method!");

		std::size_t res=0;
		auto mv = *v.begin();
		for (auto it=v.cbegin(); it!=v.cend(); it++){
			if (*it < mv){
				res = it - v.cbegin();
				mv = *it;
			}
		}
		return mv;
	}

	// the min value with comparator
	template <typename VectorT, typename Compare>
	decltype(auto) min(const VectorT & v, Compare comp){
		static_assert(type_traits::is_vector<VectorT>::value, "A Vector type requires a cbegin() and a cend() method!");

		std::size_t res=0, pos=0;
		auto mv = *v.cbegin();
		for (auto it=v.cbegin(); it!=v.cend(); it++){
			if (comp(*it, mv)){
				res = pos;//it - v.cbegin();
				mv = *it;
			}
			pos++;
		}
		return mv;
	}

	// index of the max value
	template <typename VectorT>
	std::size_t argmax(const VectorT & v){
		static_assert(type_traits::is_vector<VectorT>::value, "A Vector type requires a cbegin() and a cend() method!");

		std::size_t res=0, pos=0;
		auto mv = *v.cbegin();
		for (auto it=v.cbegin(); it!=v.cend(); it++){
			if (*it > mv){
				res = pos;
				mv = *it;
			}
			pos++;
		}
		return res;
	}

	// the max value
	template <typename VectorT>
	decltype(auto) max(const VectorT & v){
		static_assert(type_traits::is_vector<VectorT>::value, "A Vector type requires a cbegin() and a cend() method!");

		std::size_t res=0, pos=0;
		auto mv = *v.cbegin();
		for (auto it=v.cbegin(); it!=v.cend(); it++){
			if (*it > mv){
				res = pos;//it - v.cbegin();
				mv = *it;
			}
			pos++;
		}
		return mv;
	}


	// the max value with comparator
	template <typename VectorT, typename Compare>
	decltype(auto) max(const VectorT & v, Compare comp){
		static_assert(type_traits::is_vector<VectorT>::value, "A Vector type requires a cbegin() and a cend() method!");

		std::size_t res=0, pos=0;
		auto mv = *v.cbegin();
		for (auto it=v.cbegin(); it!=v.cend(); it++){
			if (!comp(*it, mv)){
				res = pos;//it - v.cbegin();
				mv = *it;
			}
			pos++;
		}
		return mv;
	}




//***********************
	// A norm is specified as f(Sum(Rule(v_i))) where
	// 			v_i 	= i-th element of the vector v
	//			Rule 	= some function applied to each vector element (e.g. abs() or square)
	//			Sum  	= the summation over the whole vector
	//			f 		= some function applied to the sum (e.g. the square root)
	template <typename NormRule, typename VectorT>
	decltype(auto) norm(const VectorT & v){
		static_assert(type_traits::is_vector<VectorT>::value, "A Vector type requires a cbegin() and a cend() method!");

		auto it = v.cbegin();
		auto val = NormRule::element_rule(*it);
		it++;
		while (it != v.cend()){
			val += NormRule::element_rule(*it);
			it++;
		}
		return NormRule::function_rule(val);
	}

	template <std::size_t p>
	struct PNormRule {
		template <typename T>
		static T element_rule(const T & val){
			return std::pow(val, p);
		};

		template <typename T>
		static T function_rule(const T & val){
			return std::pow(val, 1.0/static_cast<double>(p));
		}
	};


	template <> 
	struct PNormRule<1> {
		template <typename T>
		static T element_rule(const T & val){
			return std::abs(val);
		};

		template <typename T>
		static T function_rule(const T & val){
			return val;
		}
	};

	// some common norm functions defined for convenience
	template <typename VectorT>
	decltype(auto) norm_inf(const VectorT & v){
		static_assert(type_traits::is_vector<VectorT>::value, "A Vector type requires a cbegin() and a cend() method!");

		auto it = v.cbegin();;
		auto val = std::abs(*it);
		it++;
		while (it != v.cend()){
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

//***********************
	// use CRTP to extend vector functors defined here to any
	// class that inherits this class
	template <typename Derived>
	struct VectorFunctors {
	private:
		// static_assert(type_traits::is_vector<Derived>::value, "to inherit Vector Assignment, the vector class must be assignable!");
		
		Derived & derived() {return *static_cast<Derived *>(this);};
		const Derived & derived() const {return *static_cast<const Derived *>(this);};

	public:


	/////// CONST FUNCTORS ///////////
	// These are enabled for all vector types
		template <bool Header = false, bool AsColumn = false>
		void write(std::string dlm = " ", std::ostream & os = std::cout, std::size_t ntabs = 0) const {
			vector::write<Header, AsColumn>(derived(), dlm, os, ntabs);
		}

		template <typename VectorType>
		decltype(auto) inner_product(const VectorType & v) const { return vector::inner_product(derived(), v);}

		decltype(auto) norm_inf() const { return vector::norm_inf(derived());};
		decltype(auto) norm_1() const { return vector::norm_1(derived());};
		decltype(auto) norm_2() const { return vector::norm_2(derived());};
		decltype(auto) norm_3() const { return vector::norm_3(derived());};
		decltype(auto) norm_4() const { return vector::norm_4(derived());};
		decltype(auto) norm_5() const { return vector::norm_5(derived());};

		decltype(auto) max() const { return vector::max(derived());};
		template <typename Compare>
		decltype(auto) max(Compare comp) const { return vector::max(derived(), comp);};
		decltype(auto) min() const { return vector::min(derived());};
		template <typename Compare>
		decltype(auto) min(Compare comp) const { return vector::min(derived(), comp);};
		decltype(auto) argmax() const { return vector::argmax(derived());};
		decltype(auto) argmin() const { return vector::argmin(derived());};

		// length()

	/////// NON-CONST FUNCTORS ///////
	// These are enabled only for vector types that 
	// are assignable

		template <typename ScalarType,
				  typename T = Derived>
		typename std::enable_if<type_traits::is_assignable_vector<T>::value, void>::type 
		fill(ScalarType s){return vector::fill(derived(), s);};

		template <typename T = Derived>
		typename std::enable_if<type_traits::is_assignable_vector<T>::value, void>::type 
		fill_rand() {return vector::fill_rand(derived());};

		template <typename T = Derived>
		typename std::enable_if<type_traits::is_assignable_vector<T>::value, void>::type 
		fill_randn() {return vector::fill_randn(derived());};

		template <typename T = Derived>
		typename std::enable_if<type_traits::is_assignable_vector<T>::value
							 && type_traits::is_resizable_vector<T>::value, void>::type 
		read(std::string filename) {
			std::ifstream fs(filename.c_str(), std::ifstream::in);
			if (!fs.is_open()){
				std::cerr << "libra::vector Error attempting to open file: " << filename << std::endl;
				throw -1;
			}

			// read vector size
			std::size_t s = vector::read_to<true>(fs, derived());

			// resize self
			derived().resize(s);
			
			// read to vector
			vector::read_to(fs, derived());
			return;
		};

	};


//***********************



	// use CRTP to extend generalized vector assignment to any
	// class that inherits this class
	template <typename Derived>
	struct VectorAssignment {
	private:
		// static_assert(type_traits::is_assignable_vector<Derived>::value, "to inherit Vector Assignment, the vector class must be assignable!");
		
		Derived & derived() {return *static_cast<Derived *>(this);};
		const Derived & derived() const {return *static_cast<Derived *>(this);};

	public:


		template <typename T = Derived>
		typename std::enable_if<type_traits::is_resizable_vector<T>::value, void>::type
		resize_derived(size_type l) {derived().resize(l);};

		template <typename T = Derived>
		typename std::enable_if<!type_traits::is_resizable_vector<T>::value, void>::type
		resize_derived(size_type l) {};

		// assignment operator
		template <typename VectorT, 
				  typename T = typename std::enable_if<!std::is_same<VectorT, Derived>::value, void>::type>
		Derived & operator=(const VectorT & v){
			static_assert(type_traits::is_traversable_vector<VectorT>::value, "Vector Assignment requires a traversable vector input to operator=!");
			// std::cout << "am here non-class assignment..." << std::endl;
			resize_derived(vector::length(v));
			auto it1 = derived().begin();
			auto it2 = v.cbegin();
			while (it1 != derived().end() && it2 != v.cend()){
				(*it1) = (*it2);
				it1++;
				it2++;
			}
			return derived();
		}


		// move assignment operator
		template <typename VectorT, typename T = typename std::enable_if<!std::is_same<VectorT, Derived>::value, void>::type>
		Derived & operator=(const VectorT && v){
			static_assert(type_traits::is_traversable_vector<VectorT>::value, "Vector Assignment requires a traversable vector input to operator=!");
			// std::cout << "am here non-class move assignment..." << std::endl;
			resize_derived(vector::length(v));
			auto it1 = derived().begin();
			auto it2 = v.cbegin();
			while (it1 != derived().end() && it2 != v.cend()){
				(*it1) = (*it2);
				it1++;
				it2++;
			}
			return derived();
		}


		// add assignment operator
		template <typename VectorT, 
				  typename T = typename std::enable_if<!std::is_same<VectorT, Derived>::value, void>::type>
		Derived & operator+=(const VectorT & v){
			static_assert(type_traits::is_traversable_vector<VectorT>::value, "Vector Assignment requires a traversable vector input to operator=!");
			// std::cout << "am here non-class assignment..." << std::endl;
			// resize_derived(vector::length(v));
			auto it1 = derived().begin();
			auto it2 = v.cbegin();
			while (it1 != derived().end() && it2 != v.cend()){
				(*it1) += (*it2);
				it1++;
				it2++;
			}
			return derived();
		}

		// subtract assignment operator
		template <typename VectorT, 
				  typename T = typename std::enable_if<!std::is_same<VectorT, Derived>::value, void>::type>
		Derived & operator-=(const VectorT & v){
			static_assert(type_traits::is_traversable_vector<VectorT>::value, "Vector Assignment requires a traversable vector input to operator=!");
			// std::cout << "am here non-class assignment..." << std::endl;
			// resize_derived(vector::length(v));
			auto it1 = derived().begin();
			auto it2 = v.cbegin();
			while (it1 != derived().end() && it2 != v.cend()){
				(*it1) -= (*it2);
				it1++;
				it2++;
			}
			return derived();
		}
	};





} // end namespace vector

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