#ifndef _UNARYVECTOREXPRESSION_H
#define _UNARYVECTOREXPRESSION_H

#include <stdio.h>
#include <stdlib.h>
#include <complex>
#include <algorithm>
#include <set>
#include <iterator>

#include <utility>

#include "Traits.hpp"



// These are Unary Vector Expression templates for GENERAL VECTORS


namespace libra{
namespace vector{


	// functor is required to have a method get() that can operate on 
	// the stored type of VectorType
	template <typename VectorType, typename Functor>
	struct UnaryVectorExpression {
		static_assert(type_traits::is_traversable_vector<VectorType>::value, "UnaryVectorExpression requires traversable vector types!");
		typedef decltype(std::declval<VectorType>().cbegin()) 	iterator_type;

		UnaryVectorExpression(const VectorType & vec, Functor func) : v(vec), f(func) {};
		// UnaryVectorExpression(const VectorType && vec, Functor func) : v(vec), f(func) {};

		size_type size() const {return v.size();};

		struct const_iterator{
			typedef const_iterator					self_type;
			typedef typename type_traits::contained_type<VectorType>::type 	value_type;
			typedef std::forward_iterator_tag		iterator_category;

			const_iterator(iterator_type vit, const Functor * fn) : it(vit), fun(fn) {};

			// copy assignment
			const_iterator & operator=(const const_iterator & cit){
				const_iterator i(cit);
				std::swap(i,*this);
				return *this;
			}

			// dereferencing returns the functor acting on dereferenced iterator value
			value_type operator*() const {return fun->get(*it);};


			// increment operators
			self_type operator++(){
				it++;
				return *this;
			}
			self_type operator++(int blah){
				it++;
				return *this;
			}

			// equivalence operators
			bool operator!=(const self_type & leaf) const {return it != leaf.it;};
			bool operator==(const self_type & leaf) const {return it == leaf.it;};


		private:
			iterator_type 		it;
			const Functor * 	fun;
		};

		const_iterator cbegin() const {return const_iterator(v.cbegin(), &f);};
		const_iterator cend()	const {return const_iterator(v.cend(), &f);};
	private:
		const VectorType & v;
		Functor 	 f;
	};


// *************** C++ Math Library Functors *************
	// trigonometric functions
	struct CosFunctor{
		template <typename T>
		T get(const T & val) const {return std::cos(val);};
	};

	struct SinFunctor{
		template <typename T>
		T get(const T & val) const {return std::sin(val);};
	};

	struct TanFunctor{
		template <typename T>
		T get(const T & val) const {return std::tan(val);};
	};

	struct AcosFunctor{
		template <typename T>
		T get(const T & val) const {return std::acos(val);};
	};

	struct AsinFunctor{
		template <typename T>
		T get(const T & val) const {return std::asin(val);};
	};

	struct AtanFunctor{
		template <typename T>
		T get(const T & val) const {return std::atan(val);};
	};

	// hyperbolic trigonometric functions
	struct CoshFunctor{
		template <typename T>
		T get(const T & val) const {return std::cosh(val);};
	};

	struct SinhFunctor{
		template <typename T>
		T get(const T & val) const {return std::sinh(val);};
	};

	struct TanhFunctor{
		template <typename T>
		T get(const T & val) const {return std::tanh(val);};
	};

	struct AcoshFunctor{
		template <typename T>
		T get(const T & val) const {return std::acosh(val);};
	};

	struct AsinhFunctor{
		template <typename T>
		T get(const T & val) const {return std::asinh(val);};
	};

	struct AtanhFunctor{
		template <typename T>
		T get(const T & val) const {return std::atanh(val);};
	};

	// exponential and logarithmic functions
	struct ExpFunctor{
		template <typename T>
		T get(const T & val) const {return std::exp(val);};
	};

	struct LogFunctor{
		template <typename T>
		T get(const T & val) const {return std::log(val);};
	};

	struct Log10Functor{
		template <typename T>
		T get(const T & val) const {return std::log10(val);};
	};

	struct Exp2Functor{
		template <typename T>
		T get(const T & val) const {return std::exp2(val);};
	};

	struct Expm1Functor{
		template <typename T>
		T get(const T & val) const {return std::expm1(val);};
	};

	struct IlogbFunctor{
		template <typename T>
		T get(const T & val) const {return std::ilogb(val);};
	};

	struct Log1pFunctor{
		template <typename T>
		T get(const T & val) const {return std::log1p(val);};
	};

	struct Log2Functor{
		template <typename T>
		T get(const T & val) const {return std::log2(val);};
	};

	struct LogbFunctor{
		template <typename T>
		T get(const T & val) const {return std::logb(val);};
	};

	// power functions
	// struct PowFunctor{
	// 	template <typename T, typename T2>
	// 	T get(const T & val, const const T2 & exponent) {return std::pow(val, exponent);};
	// };

	struct SqrtFunctor{
		template <typename T>
		T get(const T & val) const {return std::sqrt(val);};
	};

	struct CbrtFunctor{
		template <typename T>
		T get(const T & val) const {return std::cbrt(val);};
	};

	struct HypotFunctor{
		template <typename T>
		T get(const T & val) const {return std::hypot(val);};
	};

	// error and gamma functions
	struct ErfFunctor{
		template <typename T>
		T get(const T & val) const {return std::erf(val);};
	};

	struct ErfcFunctor{
		template <typename T>
		T get(const T & val) const {return std::erfc(val);};
	};

	struct TgammaFunctor{
		template <typename T>
		T get(const T & val) const {return std::tgamma(val);};
	};

	struct LgammaFunctor{
		template <typename T>
		T get(const T & val) const {return std::lgamma(val);};
	};

	// rounding and remainder functions
	struct CeilFunctor{
		template <typename T>
		T get(const T & val) const {return std::ceil(val);};
	};

	struct FloorFunctor{
		template <typename T>
		T get(const T & val) const {return std::floor(val);};
	};

	// struct FmodFunctor{
	// 	template <typename T, typename T2>
	// 	T get(const T & val, const const T2 & denom) {return std::fmod(val, denom);};
	// };

	struct TruncFunctor{
		template <typename T>
		T get(const T & val) const {return std::trunc(val);};
	};

	struct RoundFunctor{
		template <typename T>
		T get(const T & val) const {return std::round(val);};
	};
	////////// TODO: THERE ARE MORE IN ROUNDING AND REMAINDER TO ADD LATER...

	// floating point manipulation functions
	// struct CopySignFunctor{
	// 	template <typename T, typename T2>
	// 	T get(const T & val, const const T2 & t) {return std::copysign(val, t);};
	// };

	// minimum, maximum, difference functions 
	// ... I don't know if these are appropriate for this location

	// other functions
	struct AbsFunctor{
		template <typename T>
		T get(const T & val) const {return std::abs(val);};
	};


	// c++17 defines a bunch of new special mathematical functions
	// (Legendre polynomials, etc...) but I'm not sure if they're 
	// useful here yet... also need to check for c++17 specification
	// at compile time


// *************** Custom Functors **********************
	// vector + scalar and scalar + vector
	template <typename ScalarType>
	struct ScalarAdditionFunctor{
		ScalarAdditionFunctor(ScalarType cval) : c(cval) {};

		template <typename T>
		T get(const T & val) const {return val + c;};

	private:
		ScalarType c;
	};

	// vector - scalar and scalar - vector


	// scalar * vector and vector * scalar
	template <typename ScalarType>
	struct ScalarMultiplicationFunctor{
		ScalarMultiplicationFunctor(ScalarType cval) : c(cval) {};

		template <typename T>
		T get(const T & val) const {return val*c;};

	private:
		ScalarType c;
	};





// *************** C++ Math Library Overloads ************
	// overload for vector types to return an expression template

	template <typename VectorType>
	UnaryVectorExpression<VectorType, CosFunctor> cos(VectorType & v){
		return UnaryVectorExpression<VectorType,CosFunctor>(v, CosFunctor()); 
	}

	template <typename VectorType>
	UnaryVectorExpression<VectorType, SinFunctor> sin(VectorType & v){
		return UnaryVectorExpression<VectorType,SinFunctor>(v, SinFunctor()); 
	}

	template <typename VectorType>
	UnaryVectorExpression<VectorType, TanFunctor> tan(VectorType & v){
		return UnaryVectorExpression<VectorType,TanFunctor>(v, TanFunctor()); 
	}

	template <typename VectorType>
	UnaryVectorExpression<VectorType, AcosFunctor> acos(VectorType & v){
		return UnaryVectorExpression<VectorType, AcosFunctor>(v, AcosFunctor());
	}

	template <typename VectorType>
	UnaryVectorExpression<VectorType, AsinFunctor> asin(VectorType & v){
		return UnaryVectorExpression<VectorType, AsinFunctor>(v, AsinFunctor());
	}

	template <typename VectorType>
	UnaryVectorExpression<VectorType, AtanFunctor> atan(VectorType & v){
		return UnaryVectorExpression<VectorType, AtanFunctor>(v, AtanFunctor());
	}

	template <typename VectorType>
	UnaryVectorExpression<VectorType, CoshFunctor> cosh(VectorType & v){
		return UnaryVectorExpression<VectorType, CoshFunctor>(v, CoshFunctor());
	}

	template <typename VectorType>
	UnaryVectorExpression<VectorType, SinhFunctor> sinh(VectorType & v){
		return UnaryVectorExpression<VectorType, SinhFunctor>(v, SinhFunctor());
	}

	template <typename VectorType>
	UnaryVectorExpression<VectorType, TanhFunctor> tanh(VectorType & v){
		return UnaryVectorExpression<VectorType, TanhFunctor>(v, TanhFunctor());
	}

	template <typename VectorType>
	UnaryVectorExpression<VectorType, AcoshFunctor> acosh(VectorType & v){
		return UnaryVectorExpression<VectorType, AcoshFunctor>(v, AcoshFunctor());
	}

	template <typename VectorType>
	UnaryVectorExpression<VectorType, AsinhFunctor> asinh(VectorType & v){
		return UnaryVectorExpression<VectorType, AsinhFunctor>(v, AsinhFunctor());
	}

	template <typename VectorType>
	UnaryVectorExpression<VectorType, AtanhFunctor> atanh(VectorType & v){
		return UnaryVectorExpression<VectorType, AtanhFunctor>(v, AtanhFunctor());
	}

	template <typename VectorType>
	UnaryVectorExpression<VectorType, ExpFunctor> exp(VectorType & v){
		return UnaryVectorExpression<VectorType,ExpFunctor>(v, ExpFunctor()); 
	}

	template <typename VectorType>
	UnaryVectorExpression<VectorType, LogFunctor> log(VectorType & v){
		return UnaryVectorExpression<VectorType,LogFunctor>(v, LogFunctor()); 
	}

	template <typename VectorType>
	UnaryVectorExpression<VectorType, Log10Functor> log10(VectorType & v){
		return UnaryVectorExpression<VectorType, Log10Functor>(v, Log10Functor());
	}

	template <typename VectorType>
	UnaryVectorExpression<VectorType, Exp2Functor> exp2(VectorType & v){
		return UnaryVectorExpression<VectorType, Exp2Functor>(v, Exp2Functor());
	}

	template <typename VectorType>
	UnaryVectorExpression<VectorType, Expm1Functor> expm1(VectorType & v){
		return UnaryVectorExpression<VectorType, Expm1Functor>(v, Expm1Functor());
	}

	template <typename VectorType>
	UnaryVectorExpression<VectorType, IlogbFunctor> ilogb(VectorType & v){
		return UnaryVectorExpression<VectorType, IlogbFunctor>(v, IlogbFunctor());
	}

	template <typename VectorType>
	UnaryVectorExpression<VectorType, Log1pFunctor> log1p(VectorType & v){
		return UnaryVectorExpression<VectorType, Log1pFunctor>(v, Log1pFunctor());
	}

	template <typename VectorType>
	UnaryVectorExpression<VectorType, Log2Functor> log2(VectorType & v){
		return UnaryVectorExpression<VectorType, Log2Functor>(v, Log2Functor());
	}

	template <typename VectorType>
	UnaryVectorExpression<VectorType, LogbFunctor> logb(VectorType & v){
		return UnaryVectorExpression<VectorType, LogbFunctor>(v, LogbFunctor());
	}

	// template <typename VectorType>
	// UnaryVectorExpression<VectorType, PowFunctor> pow(VectorType & v){
	// 	return UnaryVectorExpression<VectorType,PowFunctor>(v, PowFunctor()); 
	// }

	template <typename VectorType>
	UnaryVectorExpression<VectorType, SqrtFunctor> sqrt(VectorType & v){
		return UnaryVectorExpression<VectorType, SqrtFunctor>(v, SqrtFunctor());
	}

	template <typename VectorType>
	UnaryVectorExpression<VectorType, CbrtFunctor> cbrt(VectorType & v){
		return UnaryVectorExpression<VectorType, CbrtFunctor>(v, CbrtFunctor());
	}

	template <typename VectorType>
	UnaryVectorExpression<VectorType, HypotFunctor> hypot(VectorType & v){
		return UnaryVectorExpression<VectorType, HypotFunctor>(v, HypotFunctor());
	}

	template <typename VectorType>
	UnaryVectorExpression<VectorType, ErfFunctor> erf(VectorType & v){
		return UnaryVectorExpression<VectorType,ErfFunctor>(v, ErfFunctor()); 
	}

	template <typename VectorType>
	UnaryVectorExpression<VectorType, ErfcFunctor> erfc(VectorType & v){
		return UnaryVectorExpression<VectorType, ErfcFunctor>(v, ErfcFunctor());
	}

	template <typename VectorType>
	UnaryVectorExpression<VectorType, TgammaFunctor> tgamma(VectorType & v){
		return UnaryVectorExpression<VectorType, TgammaFunctor>(v, TgammaFunctor());
	}

	template <typename VectorType>
	UnaryVectorExpression<VectorType, LgammaFunctor> lgamma(VectorType & v){
		return UnaryVectorExpression<VectorType, LgammaFunctor>(v, LgammaFunctor());
	}

	template <typename VectorType>
	UnaryVectorExpression<VectorType, CeilFunctor> ceil(VectorType & v){
		return UnaryVectorExpression<VectorType, CeilFunctor>(v, CeilFunctor());
	}

	template <typename VectorType>
	UnaryVectorExpression<VectorType, FloorFunctor> floor(VectorType & v){
		return UnaryVectorExpression<VectorType, FloorFunctor>(v, FloorFunctor());
	}

	// template <typename VectorType>
	// UnaryVectorExpression<VectorType, FmodFunctor> fmod(VectorType & v){
	// 	return UnaryVectorExpression<VectorType, FmodFunctor>(v, FmodFunctor());
	// }

	template <typename VectorType>
	UnaryVectorExpression<VectorType, TruncFunctor> trunc(VectorType & v){
		return UnaryVectorExpression<VectorType, TruncFunctor>(v, TruncFunctor());
	}

	template <typename VectorType>
	UnaryVectorExpression<VectorType, RoundFunctor> round(VectorType & v){
		return UnaryVectorExpression<VectorType, RoundFunctor>(v, RoundFunctor());
	}

	// template <typename VectorType>
	// UnaryVectorExpression<VectorType, CopySignFunctor> copysign(VectorType & v){
	// 	return UnaryVectorExpression<VectorType, CopySignFunctor>(v, CopySignFunctor());
	// }

	template <typename VectorType>
	UnaryVectorExpression<VectorType, AbsFunctor> abs(VectorType & v){
		return UnaryVectorExpression<VectorType, AbsFunctor>(v, AbsFunctor());
	}


// *************** Custom Overloads ************

	// template <typename VectorType, typename ScalarType>
	// UnaryVectorExpression<VectorType, ScalarAdditionFunctor<ScalarType>> operator+(VectorType & v, ScalarType s){
	// 	return UnaryVectorExpression<VectorType, ScalarAdditionFunctor<ScalarType>>(v, ScalarAdditionFunctor<ScalarType>(s));
	// }

	template <typename VectorType, typename ScalarType,
			  typename T1 = typename std::enable_if<type_traits::is_vector<VectorType>::value, void>::type,
			  typename T2 = typename std::enable_if<!type_traits::is_vector<ScalarType>::value, void>::type
			  >
	const UnaryVectorExpression<VectorType, ScalarAdditionFunctor<ScalarType>> operator+(ScalarType s, const VectorType & v){
		return UnaryVectorExpression<VectorType, ScalarAdditionFunctor<ScalarType>>(v, ScalarAdditionFunctor<ScalarType>(s));
	}

	// template <typename VectorType, typename ScalarType>
	// UnaryVectorExpression<VectorType, ScalarAdditionFunctor<ScalarType>> operator+(ScalarType s, VectorType && v){
	// 	return UnaryVectorExpression<VectorType, ScalarAdditionFunctor<ScalarType>>(v, ScalarAdditionFunctor<ScalarType>(s));
	// }

	// template <typename VectorType, typename ScalarType>
	// UnaryVectorExpression<VectorType, ScalarMultiplicationFunctor<ScalarType>> operator*(VectorType & v, ScalarType s){
	// 	return UnaryVectorExpression<VectorType, ScalarMultiplicationFunctor<ScalarType>>(v, ScalarMultiplicationFunctor<ScalarType>(s));
	// }

	template <typename VectorType, typename ScalarType,
			  typename T1 = typename std::enable_if<type_traits::is_vector<VectorType>::value, void>::type,
			  typename T2 = typename std::enable_if<!type_traits::is_vector<ScalarType>::value, void>::type
			  >
	const UnaryVectorExpression<VectorType, ScalarMultiplicationFunctor<ScalarType>> operator*(ScalarType s, const VectorType & v){
		return UnaryVectorExpression<VectorType, ScalarMultiplicationFunctor<ScalarType>>(v, ScalarMultiplicationFunctor<ScalarType>(s));
	}


// *************** CRTP class extender for unary vector operators

	template <typename Derived>
	struct UnaryVectorOperators{
		Derived & derived() {return *static_cast<Derived *>(this);};
		const Derived & derived() const {return *static_cast<Derived *>(this);};

		UnaryVectorExpression<Derived, CosFunctor> cos(){
			return libra::vector::cos(derived());
		}

		
		UnaryVectorExpression<Derived, SinFunctor> sin(){
			return libra::vector::sin(derived());
		}

		
		UnaryVectorExpression<Derived, TanFunctor> tan(){
			return libra::vector::tan(derived());
		}

		
		UnaryVectorExpression<Derived, AcosFunctor> acos(){
			return libra::vector::acos(derived());
		}

		
		UnaryVectorExpression<Derived, AsinFunctor> asin(){
			return libra::vector::asin(derived());
		}

		
		UnaryVectorExpression<Derived, AtanFunctor> atan(){
			return libra::vector::atan(derived());
		}

		
		UnaryVectorExpression<Derived, CoshFunctor> cosh(){
			return libra::vector::cosh(derived());
		}

		
		UnaryVectorExpression<Derived, SinhFunctor> sinh(){
			return libra::vector::sinh(derived());
		}

		
		UnaryVectorExpression<Derived, TanhFunctor> tanh(){
			return libra::vector::tanh(derived());
		}

		
		UnaryVectorExpression<Derived, AcoshFunctor> acosh(){
			return libra::vector::acosh(derived());
		}

		
		UnaryVectorExpression<Derived, AsinhFunctor> asinh(){
			return libra::vector::asinh(derived());
		}

		
		UnaryVectorExpression<Derived, AtanhFunctor> atanh(){
			return libra::vector::atanh(derived());
		}

		
		UnaryVectorExpression<Derived, ExpFunctor> exp(){
			return libra::vector::exp(derived());
		}

		
		UnaryVectorExpression<Derived, LogFunctor> log(){
			return libra::vector::log(derived());
		}

		
		UnaryVectorExpression<Derived, Log10Functor> log10(){
			return libra::vector::log10(derived());
		}

		
		UnaryVectorExpression<Derived, Exp2Functor> exp2(){
			return libra::vector::exp2(derived());
		}

		
		UnaryVectorExpression<Derived, Expm1Functor> expm1(){
			return libra::vector::expm1(derived());
		}

		
		UnaryVectorExpression<Derived, IlogbFunctor> ilogb(){
			return libra::vector::ilogb(derived());
		}

		
		UnaryVectorExpression<Derived, Log1pFunctor> log1p(){
			return libra::vector::log1p(derived());
		}

		
		UnaryVectorExpression<Derived, Log2Functor> log2(){
			return libra::vector::log2(derived());
		}

		
		UnaryVectorExpression<Derived, LogbFunctor> logb(){
			return libra::vector::logb(derived());
		}

		
		UnaryVectorExpression<Derived, SqrtFunctor> sqrt(){
			return libra::vector::sqrt(derived());
		}

		
		UnaryVectorExpression<Derived, CbrtFunctor> cbrt(){
			return libra::vector::cbrt(derived());
		}

		
		UnaryVectorExpression<Derived, HypotFunctor> hypot(){
			return libra::vector::hypot(derived());
		}

		
		UnaryVectorExpression<Derived, ErfFunctor> erf(){
			return libra::vector::erf(derived());
		}

		
		UnaryVectorExpression<Derived, ErfcFunctor> erfc(){
			return libra::vector::erfc(derived());
		}

		
		UnaryVectorExpression<Derived, TgammaFunctor> tgamma(){
			return libra::vector::tgamma(derived());
		}

		
		UnaryVectorExpression<Derived, LgammaFunctor> lgamma(){
			return libra::vector::lgamma(derived());
		}

		
		UnaryVectorExpression<Derived, CeilFunctor> ceil(){
			return libra::vector::ceil(derived());
		}

		
		UnaryVectorExpression<Derived, FloorFunctor> floor(){
			return libra::vector::floor(derived());
		}
		
		UnaryVectorExpression<Derived, TruncFunctor> trunc(){
			return libra::vector::trunc(derived());
		}

		
		UnaryVectorExpression<Derived, RoundFunctor> round(){
			return libra::vector::round(derived());
		}


		UnaryVectorExpression<Derived, AbsFunctor> abs(){
			return libra::vector::abs(derived());
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