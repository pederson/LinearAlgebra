#ifndef _VECTOR_H
#define _VECTOR_H

#include <math.h>
#include <iostream>
#include <fstream>
#include <istream>
#include <sstream>
#include <string>
#include <stdlib.h>
#include <random>
#include <chrono>

#include "Tensor.hpp"
#include "VectorTools.hpp"
#include "UnaryVectorExpression.hpp"

namespace libra{





	// fixed-length resize policy
	template <size_type length_at_compile>
	struct VectorResizePolicy{
		template <typename VectorType>
		static void resize(VectorType & v, size_type n) {};
	};
	// dynamic size resize policy
	template <>
	struct VectorResizePolicy<dynamic_size>{
		template <typename VectorType, typename T = void>
		static typename std::enable_if<type_traits::is_resizable_vector<VectorType>::value, T>::type 
		resize(VectorType & v, size_type n) {v.resize(n);};

		template <typename VectorType, typename T = void>
		static typename std::enable_if<!type_traits::is_resizable_vector<VectorType>::value, T>::type
		resize(VectorType & v, size_type n) {};
	};





	// vector class definition
	template <typename scalar_type, size_type length_at_compile>
	class Vector : public Table<scalar_type, 1, length_at_compile>, 
				   public vector::VectorAssignment<Vector<scalar_type, length_at_compile>>,
				   public vector::UnaryVectorOperators<Vector<scalar_type, length_at_compile>>
	{
	public:
		typedef Vector<scalar_type, length_at_compile> 		SelfType;
		typedef Table<scalar_type, 1, length_at_compile> 	BaseType;
		typedef vector::VectorAssignment<SelfType> 			AssignmentType;


		// inherit the base class constructors
		// using BaseType::BaseType; 
		using BaseType::begin;
		using BaseType::end;
		using BaseType::cbegin;
		using BaseType::cend;
		using AssignmentType::operator=;


		// copy constructor for non-class types that 
		// still qualify as vectors
		template <typename VectorT, typename T = typename std::enable_if<!std::is_same<VectorT, SelfType>::value && type_traits::is_vector<VectorT>::value, void>::type>
		Vector(const VectorT & v){
			static_assert(type_traits::is_traversable_vector<VectorT>::value, "Vector Assignment requires a traversable vector input to operator=!");
			// std::cout << "am here non-class ctor..." << std::endl;
			(*this).resize(vector::length(v));
			auto it1 = (*this).begin();
			auto it2 = v.cbegin();
			while (it1 != (*this).end() && it2 != v.cend()){
				(*it1) = (*it2);
				it1++;
				it2++;
			}
		}


		// Vector(const Vector & v) {
		// 	std::cout << "am here copy ctor..."<< std::endl; 
		// 	*this = AssignmentType::operator=(v);
		// };

		// // Vector & operator=(const Vector & v) = delete;
		// Vector & operator=(const Vector && v){
		// 	std::cout << "am here move ctor..."<< std::endl;
		// };

		// decltype(auto) begin() {return BaseType::begin();};
		// decltype(auto) end() {return BaseType::end();};

		// explicitly define an initializer list constructor
		Vector(std::initializer_list<scalar_type> il){
			// std::cout << "am here ilist ctor..." << il.size() << std::endl;
			VectorResizePolicy<length_at_compile>::resize(*this, il.size());
			// std::cout << "am here after resize..." << std::endl;
			auto it = std::begin(il);
			auto itme = (*this).begin();//std::begin(*this);
			// std::cout << "got beginnings... " << (*this).size() << std::endl;
			while (it != std::end(il)){
				(*itme) = *it;
				// std::cout << "assign" << std::endl;
				it++;
				itme++;
				// std::cout << "incr" << std::endl;
			}
			// std::cout << "finished" << std::endl;
		}

		// explicitly define an empty constructor
		Vector() {};

		// vector constructor with size
		Vector(size_type l){VectorResizePolicy<length_at_compile>::resize(*this, l);};


		// void resize(size_type l) {VectorResizePolicy<length_at_compile>::resize(*this, l);};
		// 		BEWARE OF FORWARDING CONSTRUCTOR...
		// 		it can be very bad for overload resolution... basically it always gets chosen
		// // this nasty-looking code simply allows the use of initializer list for construction
		// template <typename... Args, size_type s = length_at_compile>
		// Vector(Args &&... args, typename std::enable_if<s != dynamic_size, void>::type * = 0) : BaseType{std::forward<Args>(args)...} {};

		// compound addition with arbitrary vector
		template <typename OtherVector>
		Vector & operator+=(const OtherVector & v) {
			// VectorResizePolicy<length_at_compile>::resize(out, vector::length(v));
			auto itv = std::cbegin(v);
			auto it = std::begin(*this);
			while (it != std::end(*this)){
				(*it) += (*itv);
				it++;
				itv++;
			}
			return *this;
		};

		// random access operators by index
		scalar_type & operator()(size_type i){return (*this)[i];};
		const scalar_type & operator()(size_type i) const {return (*this)[i];};

	};



// // *******************************************************************************



// *******************************************************************************



} // end namespace libra




#endif