/** @file VectorizeContainerMultifunctor.hpp
 *  @brief File with VectorizeContainerMultifunctor class
 *
 *  This contains the VectorizeContainerMultifunctor class descriptor
 *
 *  @author D. Pederson
 *  @bug No known bugs. 
 */
 
#ifndef _VECTORIZECONTAINERMULTIFUNCTOR_H
#define _VECTORIZECONTAINERMULTIFUNCTOR_H

#include <type_traits>
#include <algorithm>
#include <iostream>

#include "Macros.hpp"

 namespace libra{

/** @class VectorizeContainerMultifunctor
 *  @brief VectorizeContainerMultifunctor class to extend vector-type iteration
 *		   to a functor of the contained object iterator
 *
 *  
 *
 */
template <typename ContainerType,
		  typename StaticMethodFunctor>
class VectorizeContainerMultifunctor : public vector::VectorFunctors<VectorizeContainerMultifunctor<ContainerType, StaticMethodFunctor>>,
						   		  	   public vector::VectorAssignment<VectorizeContainerMultifunctor<ContainerType, StaticMethodFunctor>>{
public:
	typedef VectorizeContainerMultifunctor<ContainerType, StaticMethodFunctor> 		SelfType;
	typedef vector::VectorAssignment<SelfType> 									AssignmentType;
	using AssignmentType::operator=;
private:
	friend class vcm_iterator;
	ContainerType * 					mCont;
	std::vector<StaticMethodFunctor> 	mFuncts;

	template <bool is_const>
	class vcm_iterator{
	private:
		typedef typename std::conditional<is_const, 
						const VectorizeContainerMultifunctor, 
						VectorizeContainerMultifunctor>::type 	container_type;
		typedef typename ContainerType::iterator 			iterator_type;

		container_type * mVSM;
		iterator_type 	mIt;
		unsigned int 	mCtr_functs;
		unsigned int 	mCtr_objs;
	public:
		typedef vcm_iterator														self_type;
		typedef typename iterator_type::difference_type								difference_type;
		typedef decltype(std::declval<StaticMethodFunctor>().operator()(mIt)) 		value_type;
		typedef typename std::conditional<is_const, 
								  std::add_lvalue_reference_t<const std::remove_reference_t<value_type>>,
								  value_type &>::type 								reference;
		typedef typename std::conditional<is_const, 
								  const std::remove_reference_t<value_type> *,
								  std::remove_reference_t<value_type> *>::type 		pointer;
		typedef typename iterator_type::iterator_category 							iterator_category;

		vcm_iterator(container_type * vcm, iterator_type it, unsigned int fctr)
		: mVSM(vcm), mIt(it), mCtr_functs(fctr), mCtr_objs(0) {};

		// copy assignment
		vcm_iterator & operator=(const vcm_iterator & cit){
			vcm_iterator i(cit);
			std::swap(i,*this);
			return *this;
		}

		pointer operator->() const {return &mVSM->mFuncts[mCtr_functs](mIt);};
		reference operator*() const {return mVSM->mFuncts[mCtr_functs](mIt);};

		self_type operator++(){
			mCtr_functs++;
			unsigned int mod = mCtr_functs/mVSM->mFuncts.size();
			mCtr_objs += mod;
			mCtr_functs -= mod*mVSM->mFuncts.size();
			for (int i=0; i<mod; i++) mIt++; // have to do it this way because iterator_type could possibly be forward_iterator only
			return *this;
		};
		self_type operator++(int blah) {
			mCtr_functs++;
			unsigned int mod = mCtr_functs/mVSM->mFuncts.size();
			mCtr_objs += mod;
			mCtr_functs -= mod*mVSM->mFuncts.size();
			for (int i=0; i<mod; i++) mIt++;
			return *this;
		};

		bool operator!=(const self_type & leaf) const {return mIt != leaf.mIt || mCtr_functs != leaf.mCtr_functs;};
		bool operator==(const self_type & leaf) const {return mIt == leaf.mIt && mCtr_functs == leaf.mCtr_functs;};
	};
public:
	typedef vcm_iterator<true> 		const_iterator;
	typedef vcm_iterator<false> 	iterator;

	VectorizeContainerMultifunctor(ContainerType & c, std::vector<StaticMethodFunctor> v) : mCont(&c), mFuncts(v){};

	decltype(mCont->size()) size() const {return mFuncts.size()*mCont->size();};

	iterator begin() {return iterator(this, mCont->begin(), 0);};
	iterator end()	 {return iterator(this, mCont->end(), 0);};

	const_iterator cbegin() const {return const_iterator(this, mCont->begin(), 0);};
	const_iterator cend() const	 {return const_iterator(this, mCont->end(), 0);};
};






template <typename ContainerT, typename Functor>
using VCM = libra::VectorizeContainerMultifunctor<ContainerT, Functor>;				


#define LIBRA_ASSERT_EXISTENCE(FunctionName) 				\
	static_assert(LIBRA_HAS_METHOD(FunctionName)< 			\
				  decltype(*derived().begin())>::value 		\
				  ,LIBRA_STRINGIZE(FunctionName)); 


#define LIBRA_ASSERT_EXISTENCE_INTERFACE(FunctionName) 		\
	static_assert(LIBRA_HAS_METHOD(FunctionName)< 			\
				  decltype(InterfacePolicy::get(*derived().begin()))>::value 		\
				  ,LIBRA_STRINGIZE(FunctionName)); 


#define LIBRA_VECTORIZE_MULTIFUNCTOR(ResultName) CRTP_Vectorize_Multifunctor_##ResultName
#define LIBRA_VECTORIZE_MULTIFUNCTOR_DEF(ResultName, FunctorType, ...)					\
																						\
	namespace libra{																	\
		namespace detail{																\
			namespace vectorize_multifunctor{											\
				namespace ResultName {													\
					LIBRA_FOR_EACH_SEP(LIBRA_FUNCTOR_PREPARE, __VA_ARGS__);				\
					std::vector<FunctorType> functs = 									\
					{LIBRA_FOR_EACH_SEQ(LIBRA_INSTANTIATE_FUNCTOR_FOR, __VA_ARGS__)};	\
																						\
																						\
					template <typename Derived> 										\
					struct LIBRA_VECTORIZE_MULTIFUNCTOR(ResultName){ 					\
					private: 															\
						Derived & derived() {return *static_cast<Derived *>(this);};	\
						const Derived & derived() const {return *static_cast<const Derived *>(this);};	\
					public: 															\
						libra::VCM<Derived, FunctorType> 								\
						ResultName() 													\
						{																\
							LIBRA_FOR_EACH_SEP(LIBRA_ASSERT_EXISTENCE, __VA_ARGS__);	\
							return libra::VCM<Derived, FunctorType>(derived(), 			\
									functs);};											\
					};																	\
				}																		\
			} 																			\
		}																				\
	}																					\
																						\
	using libra::detail::vectorize_multifunctor::ResultName::LIBRA_VECTORIZE_MULTIFUNCTOR(ResultName);	\






#define LIBRA_VECTORIZE_MULTIFUNCTOR_INTERFACE_DEF(ResultName, InterfaceName, FunctorType, ...)		\
																						\
	namespace libra{																	\
		namespace detail{																\
			namespace vectorize_multifunctor{											\
				namespace ResultName {													\
					typedef InterfaceName InterfacePolicy;								\
					LIBRA_FOR_EACH_SEP(LIBRA_FUNCTOR_PREPARE_INTERFACE, __VA_ARGS__);	\
					std::vector<FunctorType> functs = 									\
					{LIBRA_FOR_EACH_SEQ(LIBRA_INSTANTIATE_FUNCTOR_FOR, __VA_ARGS__)};	\
																						\
																						\
					template <typename Derived> 										\
					struct LIBRA_VECTORIZE_MULTIFUNCTOR(ResultName){ 					\
					private: 															\
						Derived & derived() {return *static_cast<Derived *>(this);};	\
						const Derived & derived() const {return *static_cast<const Derived *>(this);};	\
					public: 															\
						libra::VCM<Derived, FunctorType> 								\
						ResultName() 													\
						{																\
							LIBRA_FOR_EACH_SEP(LIBRA_ASSERT_EXISTENCE_INTERFACE, __VA_ARGS__);	\
							return libra::VCM<Derived, FunctorType>(derived(), 			\
									functs);};											\
					};																	\
				}																		\
			} 																			\
		}																				\
	}																					\
																						\
	using libra::detail::vectorize_multifunctor::ResultName::LIBRA_VECTORIZE_MULTIFUNCTOR(ResultName);	\




} // end namespace simbox
#endif