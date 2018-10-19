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
#include <functional>

#include "Macros.hpp"
#include "Traits.hpp"

 namespace libra{



 	//***********************
	// an assignable vector requires the non-const begin() and end()
	// iterator accessors
	template<typename StaticFunctor, typename IteratorType, typename _ = void>
	struct can_apply_static_functor : std::false_type {};

	template<typename StaticFunctor, typename IteratorType>
	struct can_apply_static_functor<
	        StaticFunctor,
	        IteratorType,
	        std::conditional_t<
	            false,
	            type_traits::is_vector_helper<
	                decltype(std::declval<StaticFunctor>().operator()(std::declval<IteratorType>()))
	                >,
	            void
	            >
	        > : public std::true_type {};



/** @class VectorizeContainerMultifunctor
 *  @brief VectorizeContainerMultifunctor class to extend vector-type iteration
 *		   to a functor of the contained object iterator
 *
 *  
 *
 */
template <typename ContainerType,
		  typename... StaticMethodFunctor>
class VectorizeContainerMultifunctor : public vector::VectorFunctors<VectorizeContainerMultifunctor<ContainerType, StaticMethodFunctor...>>,
						   		  	   public vector::VectorAssignment<VectorizeContainerMultifunctor<ContainerType, StaticMethodFunctor...>>{
public:
	typedef VectorizeContainerMultifunctor<ContainerType, StaticMethodFunctor...> 		SelfType;
	typedef vector::VectorAssignment<SelfType> 									AssignmentType;
	using AssignmentType::operator=;
private:
	typedef std::tuple<StaticMethodFunctor...> functor_tuple;
	typedef typename std::tuple_element<0, functor_tuple>::type first_functor;
	friend class vcm_iterator;

	ContainerType * 					mCont;
	// std::vector<StaticMethodFunctor> 	mFuncts;

	template <bool is_const>
	class vcm_iterator{
	private:
		typedef typename std::conditional<is_const, 
						const VectorizeContainerMultifunctor, 
						VectorizeContainerMultifunctor>::type 			container_type;
		typedef typename std::conditional<is_const, 
						typename ContainerType::const_iterator, 
						typename ContainerType::iterator>::type 		iterator_type;


		// static_assert(can_apply_static_functor<StaticMethodFunctor, iterator_type>::value, 
		// 			  "Must be able to apply static functor to iterator type");

		container_type * mVSM;
		iterator_type 	mIt;
		unsigned int 	mCtr_functs;
		unsigned int 	mCtr_objs;
	public:
		typedef vcm_iterator														self_type;
		typedef typename iterator_type::difference_type								difference_type;
		typedef std::remove_reference_t<decltype(std::declval<first_functor>().operator()(mIt))> 		value_type;
		typedef typename std::conditional<is_const, 
								  std::add_lvalue_reference_t<const std::remove_reference_t<value_type>>,
								  value_type &>::type 								reference;
		typedef typename std::conditional<is_const, 
								  const std::remove_reference_t<value_type> *,
								  std::remove_reference_t<value_type> *>::type 		pointer;
		typedef typename iterator_type::iterator_category 							iterator_category;

	private:
		typedef std::function<reference(iterator_type &)> functor_type;
		std::vector<functor_type> mFuncts;// = {StaticMethodFunctor()...}; //  instantiate into seperate generic functors
	public:
		vcm_iterator() {
			// std::cout << "constructed empty vcm_iterator" << std::endl;
		};

		vcm_iterator(vcm_iterator && v)
		: mVSM(v.mVSM), mIt(v.mIt), mCtr_functs(v.mCtr_functs), mCtr_objs(v.mCtr_objs), mFuncts(v.mFuncts) {
			// std::cout << "move-constructed vcm_iterator" << std::endl;
			// static_assert(std::is_same<bool, iterator_type>::value, "check 4 5");
			// *this = std::move(v);
		};

		vcm_iterator(const vcm_iterator & v)
		: mVSM(v.mVSM), mIt(v.mIt), mCtr_functs(v.mCtr_functs), mCtr_objs(v.mCtr_objs), mFuncts(v.mFuncts) {
			// std::cout << "copy-constructed vcm_iterator" << std::endl;
			// static_assert(std::is_same<bool, iterator_type>::value, "check 4 5");

		};

		vcm_iterator(container_type * vcm, iterator_type it, unsigned int fctr)
		: mVSM(vcm), mIt(it), mCtr_functs(fctr), mCtr_objs(0), mFuncts({StaticMethodFunctor()...}) {
			// std::cout << "constructed vcm_iterator" << std::endl;
		};

		// copy assignment
		vcm_iterator & operator=(const vcm_iterator & cit){
			if (this != &cit){
				// std::cout << "copy assign vcm_iterator" << std::endl;
				vcm_iterator i(cit);
				// std::swap(i,*this);
				std::swap(mVSM, i.mVSM); 
				std::swap(mIt, i.mIt); 
				std::swap(mCtr_functs, i.mCtr_functs); 
				std::swap(mCtr_objs, i.mCtr_objs); 
				std::swap(mFuncts, i.mFuncts);	
			}
			return *this;
		}

		// move assignment
		vcm_iterator & operator=(vcm_iterator && cit){
			if (this != &cit){
				// std::cout << "move assign vcm_iterator" << std::endl;
				vcm_iterator i(std::move(cit));
				// std::swap(i,*this);
				std::swap(mVSM, i.mVSM); 
				std::swap(mIt, i.mIt); 
				std::swap(mCtr_functs, i.mCtr_functs); 
				std::swap(mCtr_objs, i.mCtr_objs); 
				std::swap(mFuncts, i.mFuncts);			}
			// std::cout << "exiting move-assign" << std::endl;
			return *this;
		}

		pointer operator->() {return &mFuncts[mCtr_functs](mIt);};
		reference operator*()  {return mFuncts[mCtr_functs](mIt);};

		self_type operator++(){
			mCtr_functs++;
			unsigned int mod = mCtr_functs/mFuncts.size();
			mCtr_objs += mod;
			mCtr_functs -= mod*mFuncts.size();
			for (int i=0; i<mod; i++) mIt++; // have to do it this way because iterator_type could possibly be forward_iterator only
			return *this;
		};
		self_type operator++(int blah) {
			mCtr_functs++;
			unsigned int mod = mCtr_functs/mFuncts.size();
			mCtr_objs += mod;
			mCtr_functs -= mod*mFuncts.size();
			for (int i=0; i<mod; i++) mIt++;
			return *this;
		};

		bool operator!=(const self_type & leaf) const {return mIt != leaf.mIt || mCtr_functs != leaf.mCtr_functs;};
		bool operator==(const self_type & leaf) const {return mIt == leaf.mIt && mCtr_functs == leaf.mCtr_functs;};
	};
public:
	typedef vcm_iterator<true> 		const_iterator;
	typedef vcm_iterator<false> 	iterator;

	VectorizeContainerMultifunctor(ContainerType & c) : mCont(&c){};

	decltype(mCont->size()) size() const {return std::tuple_size<functor_tuple>::value*mCont->size();};

	iterator begin() {return iterator(this, mCont->begin(), 0);};
	iterator end()	 {return iterator(this, mCont->end(), 0);};

	const_iterator cbegin() const {return const_iterator(this, mCont->cbegin(), 0);};
	const_iterator cend() const	 {return const_iterator(this, mCont->cend(), 0);};
};






template <typename ContainerT, typename... Functor>
using VCM = libra::VectorizeContainerMultifunctor<ContainerT, Functor...>;				


#define LIBRA_ASSERT_EXISTENCE(FunctionName) 				\
	static_assert(LIBRA_HAS_METHOD(FunctionName)< 			\
				  decltype(*derived().begin())>::value 		\
				  ,"Iterator does not have the requested function"); 			\
	static_assert(LIBRA_HAS_METHOD(FunctionName)< 			\
				  decltype(*derived().cbegin())>::value 	\
				  ,"Const Iterator does not have the requested function"); 


#define LIBRA_ASSERT_EXISTENCE_INTERFACE(FunctionName) 		\
	static_assert(LIBRA_HAS_METHOD(FunctionName)< 			\
				  decltype(InterfacePolicy::get(*derived().begin()))>::value 		\
				  ,"Iterator does not have the requested function"); 			\
	static_assert(LIBRA_HAS_METHOD(FunctionName)< 			\
				  decltype(InterfacePolicy::get(*derived().cbegin()))>::value 		\
				  ,"Const Iterator does not have the requested function");


#define LIBRA_VECTORIZE_MULTIFUNCTOR_NO_INTERFACE(ResultName) CRTP_Vectorize_Multifunctor_Named_ ##ResultName
#define LIBRA_VECTORIZE_MULTIFUNCTOR_INTERFACE(ResultName, InterfaceName) CRTP_Vectorize_Multifunctor_Named_ ##ResultName ##_With_Interface_ ##InterfaceName
#define LIBRA_VECTORIZE_MULTIFUNCTOR(...) LIBRA_GET_MACRO(__VA_ARGS__, LIBRA_VECTORIZE_MULTIFUNCTOR_INTERFACE, LIBRA_VECTORIZE_MULTIFUNCTOR_NO_INTERFACE)(__VA_ARGS__)


#define LIBRA_VECTORIZE_MULTIFUNCTOR_DEF(ResultName, ...)								\
																						\
	namespace libra{																	\
		namespace detail{																\
			namespace vectorize_multifunctor{											\
				namespace ResultName {													\
					LIBRA_FOR_EACH_SEP(LIBRA_FUNCTOR_PREPARE, __VA_ARGS__);				\
																						\
					typedef std::tuple<LIBRA_FOR_EACH_SEQ(LIBRA_FUNCTOR_FOR, __VA_ARGS__)> functor_tuple;	\
					functor_tuple ftup;													\
																						\
					template <typename Derived> 										\
					struct LIBRA_VECTORIZE_MULTIFUNCTOR(ResultName){ 					\
					private: 															\
						Derived & derived() {return *static_cast<Derived *>(this);};	\
						const Derived & derived() const {return *static_cast<const Derived *>(this);};	\
					public: 															\
																						\
						decltype(auto) ResultName() 									\
						{																\
							LIBRA_FOR_EACH_SEP(LIBRA_ASSERT_EXISTENCE, __VA_ARGS__);	\
							return libra::VCM<Derived, LIBRA_FOR_EACH_SEQ(LIBRA_FUNCTOR_FOR, __VA_ARGS__)>(derived());	\
						};																\
					};																	\
				}																		\
			} 																			\
		}																				\
	}																					\
																						\
	using libra::detail::vectorize_multifunctor::ResultName::LIBRA_VECTORIZE_MULTIFUNCTOR(ResultName);	\




#define LIBRA_VECTORIZE_MULTIFUNCTOR_INTERFACE_DEF(ResultName, InterfaceName, ...)			\
	namespace libra{																		\
		namespace detail{																	\
			namespace vectorize_multifunctor{												\
				namespace ResultName{														\
					namespace InterfaceName_space {											\
						typedef InterfaceName InterfacePolicy;								\
						LIBRA_FOR_EACH_SEP(LIBRA_FUNCTOR_PREPARE_INTERFACE, __VA_ARGS__);	\
																							\
																							\
						template <typename Derived> 										\
						struct LIBRA_VECTORIZE_MULTIFUNCTOR(ResultName, InterfaceName){ 	\
						private: 															\
							Derived & derived() {return *static_cast<Derived *>(this);};	\
							const Derived & derived() const {return *static_cast<const Derived *>(this);};	\
						public: 															\
																							\
							decltype(auto) ResultName() 									\
							{																\
								LIBRA_FOR_EACH_SEP(LIBRA_ASSERT_EXISTENCE_INTERFACE, __VA_ARGS__);			\
								return libra::VCM<Derived, LIBRA_FOR_EACH_SEQ(LIBRA_FUNCTOR_FOR, __VA_ARGS__)>(derived());	\
							};																\
						};																	\
					}																		\
				}																			\
			} 																				\
		}																					\
	}																						\
																							\
	using libra::detail::vectorize_multifunctor::ResultName::InterfaceName_space::LIBRA_VECTORIZE_MULTIFUNCTOR(ResultName,InterfaceName);	




} // end namespace simbox
#endif