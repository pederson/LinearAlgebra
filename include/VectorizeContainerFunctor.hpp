/** @file VectorizeContainerFunctor.hpp
 *  @brief File with VectorizeContainerFunctor class
 *
 *  This contains the VectorizeContainerFunctor class descriptor
 *
 *  @author D. Pederson
 *  @bug No known bugs. 
 */
 
#ifndef _VECTORIZECONTAINERFUNCTOR_H
#define _VECTORIZECONTAINERFUNCTOR_H

#include <type_traits>
#include <algorithm>

#include "Macros.hpp"
#include "VectorTools.hpp"

 namespace libra{



/** @class VectorizeContainerFunctor
 *  @brief VectorizeContainerFunctor class to extend vector-type iteration
 *		   to a functor of the contained object iterator
 *
 *  
 *
 */
template <typename ContainerType,
		  typename StaticMethodFunctor>
class VectorizeContainerFunctor : public vector::VectorFunctors<VectorizeContainerFunctor<ContainerType, StaticMethodFunctor>>,
						   		  public vector::VectorAssignment<VectorizeContainerFunctor<ContainerType, StaticMethodFunctor>>{
public:
	typedef VectorizeContainerFunctor<ContainerType, StaticMethodFunctor> 		SelfType;
	typedef vector::VectorAssignment<SelfType> 									AssignmentType;
	using AssignmentType::operator=;
private:
	ContainerType * mCont;

	template <bool is_const>
	class vcf_iterator{
	private:
		typedef typename std::conditional<is_const, 
						const VectorizeContainerFunctor, 
						VectorizeContainerFunctor>::type 	container_type;
		typedef typename std::conditional<is_const, 
				typename ContainerType::const_iterator, 
				typename ContainerType::iterator>::type 	iterator_type;

		container_type * mVCM;
		iterator_type 	mIt;
	public:
		typedef vcf_iterator								self_type;
		typedef typename iterator_type::difference_type		difference_type;
		typedef std::remove_reference_t<decltype(StaticMethodFunctor::get(mIt))> 	value_type;
		typedef typename std::conditional<is_const, 
						 const value_type,
						 value_type>::type &				reference;
		typedef typename std::conditional<is_const, 
						 const value_type,
						 value_type>::type * 				pointer;
		typedef typename iterator_type::iterator_category 	iterator_category;

		vcf_iterator() {};

		vcf_iterator(const vcf_iterator & v)
		: mVCM(v.mVCM), mIt(v.mIt) {};

		vcf_iterator(container_type * vcm, iterator_type it)
		: mVCM(vcm), mIt(it) {};

		// copy assignment
		vcf_iterator & operator=(const vcf_iterator & cit){
			vcf_iterator i(cit);
			std::swap(i,*this);
			return *this;
		}

		pointer operator->() const {return &StaticMethodFunctor::get(mIt);};
		reference operator*() const {return StaticMethodFunctor::get(mIt);};

		self_type operator++(){mIt++; return *this;};
		self_type operator++(int blah) {mIt++; return *this;};

		bool operator!=(const self_type & leaf) const {return mIt != leaf.mIt;};
		bool operator==(const self_type & leaf) const {return mIt == leaf.mIt;};

	};
public:
	typedef vcf_iterator<true> 		const_iterator;
	typedef vcf_iterator<false> 	iterator;

	VectorizeContainerFunctor(ContainerType & c) : mCont(&c) {};

	decltype(mCont->size()) size() const {return mCont->size();};

	iterator begin() {return iterator(this, mCont->begin());};
	iterator end()	 {return iterator(this, mCont->end());};

	const_iterator cbegin() const {return const_iterator(this, mCont->cbegin());};
	const_iterator cend() const	 {return const_iterator(this, mCont->cend());};


};




template <typename ContainerT, typename Functor>
using VCF = libra::VectorizeContainerFunctor<ContainerT, Functor>;



// This takes a member function of a value_type that is held in a container
// and adds a vectorized version of that member function to the 
// entire container. 
//
#define LIBRA_VECTORIZE_FUNCTOR(ResultName) CRTP_Vectorize_Functor_##ResultName

#define LIBRA_VECTORIZE_FUNCTOR_DEF_NOINTERFACE(ResultName, FunctionName)			\
																		\
	namespace libra{													\
		namespace detail{												\
			namespace vectorize_functor{								\
				namespace ResultName {								\
					LIBRA_FUNCTOR_PREPARE(FunctionName);				\
																		\
					template <typename Derived> 						\
					struct LIBRA_VECTORIZE_FUNCTOR(ResultName){ 		\
					private: 											\
						Derived & derived() {return *static_cast<Derived *>(this);};	\
						const Derived & derived() const {return *static_cast<const Derived *>(this);};	\
																										\
					public: 																			\
						libra::VCF<Derived, LIBRA_FUNCTOR_FOR(FunctionName)> 							\
						ResultName() 																	\
						{																				\
									static_assert(LIBRA_HAS_METHOD(FunctionName)< 						\
									  decltype(*derived().begin())>::value, 							\
									  LIBRA_STRINGIZE(FunctionName)); 									\
							return libra::VCF<Derived, LIBRA_FUNCTOR_FOR(FunctionName)>(derived());};	\
					};																					\
				}														\
			} 															\
		}																\
	}																	\
																		\
	using libra::detail::vectorize_functor::ResultName::LIBRA_VECTORIZE_FUNCTOR(ResultName);		\






// This takes a member function of a value_type that is held in a container
// and adds a vectorized version of that member function to the 
// entire container. 
//
// This version adds an interface for when e.g. the value_type is a std::pair
// and you need to call ".second" before you can access the member function
#define LIBRA_VECTORIZE_FUNCTOR_DEF_INTERFACE(ResultName, FunctionName, InterfaceName)		\
																				\
	namespace libra{															\
		namespace detail{														\
			namespace vectorize_functor{										\
				namespace ResultName {										\
					typedef InterfaceName InterfacePolicy;						\
					LIBRA_FUNCTOR_PREPARE_INTERFACE(FunctionName);				\
																				\
					template <typename Derived> 								\
					struct LIBRA_VECTORIZE_FUNCTOR(ResultName){ 				\
					private: 													\
						Derived & derived() {return *static_cast<Derived *>(this);};	\
						const Derived & derived() const {return *static_cast<const Derived *>(this);};	\
					public: 																			\
						libra::VCF<Derived, LIBRA_FUNCTOR_FOR(FunctionName)> 	\
						ResultName() 											\
						{														\
							static_assert(LIBRA_HAS_METHOD(FunctionName)< 		\
									  decltype(InterfaceName::get(*derived().begin()))>::value, \
									  LIBRA_STRINGIZE(FunctionName)); 			\
							return libra::VCF<Derived, LIBRA_FUNCTOR_FOR(FunctionName)>(derived());};	\
					};															\
				}																\
			} 																	\
		}																		\
	}																			\
																				\
	using libra::detail::vectorize_functor::ResultName::LIBRA_VECTORIZE_FUNCTOR(ResultName);		\





// define an overloaded version of VECTORIZE_FUNCTOR_DEF
#define LIBRA_VECTORIZE_FUNCTOR_DEF(...) LIBRA_GET_MACRO(__VA_ARGS__, LIBRA_VECTORIZE_FUNCTOR_DEF_INTERFACE, LIBRA_VECTORIZE_FUNCTOR_DEF_NOINTERFACE)(LIBRA_GET_FIRST(__VA_ARGS__), __VA_ARGS__)
#define LIBRA_VECTORIZE_FUNCTOR_RENAME_DEF(ResultName, ...) LIBRA_GET_MACRO(__VA_ARGS__, LIBRA_VECTORIZE_FUNCTOR_DEF_INTERFACE, LIBRA_VECTORIZE_FUNCTOR_DEF_NOINTERFACE)(ResultName, __VA_ARGS__)



} // end namespace libra
#endif