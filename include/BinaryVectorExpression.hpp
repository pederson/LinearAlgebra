#ifndef _BINARYVECTOREXPRESSION_H
#define _BINARYVECTOREXPRESSION_H

#include <stdio.h>
#include <stdlib.h>
#include <complex>
#include <algorithm>
#include <set>
#include <iterator>

#include <utility>

#include "Traits.hpp"



// These are Binary Vector Expression templates for GENERAL VECTORS


namespace libra{
namespace vector{


	// functor is required to have a method get() that can operate on 
	// the stored type of VectorType
	template <typename VectorType1, typename VectorType2, typename Functor>
	struct BinaryVectorExpression {
		static_assert(type_traits::is_traversable_vector<VectorType1>::value, "BinaryVectorExpression requires traversable vector types!");
		static_assert(type_traits::is_traversable_vector<VectorType2>::value, "BinaryVectorExpression requires traversable vector types!");
		typedef decltype(std::declval<VectorType1>().cbegin()) 	iterator_type1;
		typedef decltype(std::declval<VectorType2>().cbegin()) 	iterator_type2;

		BinaryVectorExpression(const VectorType1 & vec1, const VectorType2 & vec2, Functor func) : v1(vec1), v2(vec2), f(func) {};
		// BinaryVectorExpression(const VectorType && vec, Functor func) : v(vec), f(func) {};

		size_type size() const {return v1.size();};

		struct const_iterator{
			typedef const_iterator					self_type;
			typedef typename type_traits::contained_type<VectorType1>::type 	value_type1;
			typedef typename type_traits::contained_type<VectorType2>::type 	value_type2;
			typedef decltype(std::declval<Functor>().get(std::declval<value_type1>(), std::declval<value_type2>())) 								value_type;
			typedef std::forward_iterator_tag		iterator_category;

			const_iterator(iterator_type1 vit1, iterator_type2 vit2, const Functor * fn) : it1(vit1), it2(vit2), fun(fn) {};

			// copy assignment
			const_iterator & operator=(const const_iterator & cit){
				const_iterator i(cit);
				std::swap(i,*this);
				return *this;
			}

			// dereferencing returns the functor acting on dereferenced iterator value
			value_type operator*() const {return fun->get(*it1, *it2);};


			// increment operators
			self_type operator++(){
				it1++; it2++;
				return *this;
			}
			self_type operator++(int blah){
				it1++; it2++;
				return *this;
			}

			// equivalence operators
			bool operator!=(const self_type & leaf) const {return it1 != leaf.it1 && it2 != leaf.it2;};
			bool operator==(const self_type & leaf) const {return it1 == leaf.it1 && it2 == leaf.it2;};


		private:
			iterator_type1 		it1;
			iterator_type2		it2;
			const Functor * 	fun;
		};

		const_iterator cbegin() const {return const_iterator(v1.cbegin(), v2.cbegin(), &f);};
		const_iterator cend()	const {return const_iterator(v1.cend(), v2.cend(), &f);};
	private:
		const VectorType1 & v1;
		const VectorType2 & v2;		
		Functor 	 f;
	};


// *************** C++ Math Library Functors *************
	


// *************** Custom Functors **********************
	// vector + vector
	struct VectorAdditionFunctor{
		template <typename T1, typename T2>
		typename type_traits::sum_type<T1, T2>::type get(const T1 & val1, const T2 & val2) const {return val1 + val2;};
	};

	// vector - vector
	struct VectorSubtractionFunctor{
		template <typename T1, typename T2>
		typename type_traits::difference_type<T1, T2>::type get(const T1 & val1, const T2 & val2) const {return val1 - val2;};
	};





// *************** C++ Math Library Overloads ************
	// overload for vector types to return an expression template


// *************** Custom Overloads ************

	template <typename VectorType1, typename VectorType2, 
			  typename T1 = typename std::enable_if<type_traits::is_vector<VectorType1>::value, void>::type,
			  typename T2 = typename std::enable_if<type_traits::is_vector<VectorType2>::value, void>::type
			  >
	const BinaryVectorExpression<VectorType1, VectorType2, VectorAdditionFunctor> operator+(const VectorType1 & v1, const VectorType2 & v2){
		return BinaryVectorExpression<VectorType1, VectorType2, VectorAdditionFunctor>(v1, v2, VectorAdditionFunctor());
	}

	template <typename VectorType1, typename VectorType2, 
			  typename T1 = typename std::enable_if<type_traits::is_vector<VectorType1>::value, void>::type,
			  typename T2 = typename std::enable_if<type_traits::is_vector<VectorType2>::value, void>::type
			  >
	const BinaryVectorExpression<VectorType1, VectorType2, VectorSubtractionFunctor> operator-(const VectorType1 & v1, const VectorType2 & v2){
		return BinaryVectorExpression<VectorType1, VectorType2, VectorSubtractionFunctor>(v1, v2, VectorSubtractionFunctor());
	}


// *************** CRTP class extender for unary vector operators

	template <typename Derived>
	struct BinaryVectorOperators{
		Derived & derived() {return *static_cast<Derived *>(this);};
		const Derived & derived() const {return *static_cast<Derived *>(this);};

		// BinaryVectorExpression<Derived, CosFunctor> cos(){
		// 	return libra::vector::cos(derived());
		// }
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