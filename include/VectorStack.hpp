/** @file VectorStack.hpp
 *  @brief File with VectorStack class
 *
 *  This contains the VectorStack class descriptor
 *
 *  @author D. Pederson
 *  @bug No known bugs. 
 */
 
#ifndef _VECTORSTACK_H
#define _VECTORSTACK_H

#include <type_traits>
#include <algorithm>

#include <cassert>

#include "Macros.hpp"
#include "VectorTools.hpp"

namespace libra{


namespace detail{
	template<typename T, typename... Ts>
	struct LastTypeOf {
	   typedef typename LastTypeOf<Ts...>::type type;
	};

	template<typename T>
	struct LastTypeOf<T> {
	  typedef T type;
	};



	template<typename T, typename... Ts>
	struct FirstTypeOf {
	   typedef T type;
	};

	template<typename T>
	struct FirstTypeOf<T> {
	  typedef T type;
	};




	// tuple runtime_get that 
	//*************ONLY WORKS IF ALL TYPES IN THE TUPLE
	// ARE THE SAME*******************
	template <typename Tuple,
	  		  typename Indices = std::make_index_sequence<std::tuple_size<Tuple>::value>>
	struct runtime_get_func_table;

	template<typename Tuple,size_t ... Indices>
	struct runtime_get_func_table<Tuple,std::index_sequence<Indices...>>{
	    using return_type 		= typename std::tuple_element<0,Tuple>::type &;
	    using get_func_ptr 		= return_type (*)(Tuple&) noexcept;
	    
	    static constexpr get_func_ptr table[std::tuple_size<Tuple>::value] = {
	        &std::get<Indices>...
	    };
	};

	template <typename Tuple, size_t ... Indices>
	constexpr typename
	runtime_get_func_table<Tuple,std::index_sequence<Indices...>>::get_func_ptr
	runtime_get_func_table<Tuple,std::index_sequence<Indices...>>::table[std::tuple_size<Tuple>::value];

	template<typename Tuple>
	constexpr
	typename std::tuple_element<0, typename std::remove_reference<Tuple>::type>::type&
	runtime_get(Tuple&& t, size_t index){
	    using tuple_type 		= typename std::remove_reference<Tuple>::type;
	    if(index>=std::tuple_size<tuple_type>::value)
	        throw std::runtime_error("Out of range");
	    return runtime_get_func_table<tuple_type>::table[index](t);
	}




	// tuple visitor... less-beautiful but safe implementation
	// of a tuple visitor to apply some functor to each tuple type
	template<std::size_t I>
    struct visit_impl
    {
        template<typename Tuple, typename F, typename ...Args>
        inline static constexpr int visit(Tuple &tuple, std::size_t idx, F fun, Args &&...args) noexcept(noexcept(fun(std::get<I - 1U>(tuple), std::forward<Args>(args)...)) && noexcept(visit_impl<I - 1U>::visit(tuple, idx, fun, std::forward<Args>(args)...)))
        {
            return (idx == (I - 1U) ? (fun(std::get<I - 1U>(tuple), std::forward<Args>(args)...), void(), 0) : visit_impl<I - 1U>::visit(tuple, idx, fun, std::forward<Args>(args)...));
        }

        template<typename R, typename Tuple, typename F, typename ...Args>
        inline static constexpr R visit(Tuple &tuple, std::size_t idx, F fun, Args &&...args) noexcept(noexcept(fun(std::get<I - 1U>(tuple), std::forward<Args>(args)...)) && noexcept(visit_impl<I - 1U>::template visit<R>(tuple, idx, fun, std::forward<Args>(args)...)))
        {
            return (idx == (I - 1U) ? fun(std::get<I - 1U>(tuple), std::forward<Args>(args)...) : visit_impl<I - 1U>::template visit<R>(tuple, idx, fun, std::forward<Args>(args)...));
        }

        // with 2 tuples
        template <typename Tuple1, typename Tuple2, typename F>
        inline static constexpr int visit(Tuple1 &tup1, Tuple2 & tup2, std::size_t idx, F fun) //noexcept(noexcept(fun(std::get<I - 1U>(tuple), std::forward<Args>(args)...)) && noexcept(visit_impl<I - 1U>::visit(tuple, idx, fun, std::forward<Args>(args)...)))
        {
            return (idx == (I - 1U) ? (fun(std::get<I - 1U>(tup1), std::get<I - 1U>(tup2)), void(), 0) : visit_impl<I - 1U>::visit(tup1, tup2, idx, fun));
        }
    };

    template<>
    struct visit_impl<0U>
    {
        template<typename Tuple, typename F, typename ...Args>
        inline static constexpr int visit(Tuple &, std::size_t, F, Args&&...) noexcept
        {
            return 0;
        }

        template<typename R, typename Tuple, typename F, typename ...Args>
        inline static constexpr R visit(Tuple &, std::size_t, F, Args&&...) noexcept(noexcept(R{}))
        {
            static_assert(std::is_default_constructible<R>::value, "Explicit return type of visit_at method must be default-constructible");
            return R{};
        }

        // with 2 tuples
        template <typename Tuple1, typename Tuple2, typename F>
        inline static constexpr int visit(Tuple1 &tup1, Tuple2 & tup2, std::size_t idx, F fun) //noexcept(noexcept(fun(std::get<I - 1U>(tuple), std::forward<Args>(args)...)) && noexcept(visit_impl<I - 1U>::visit(tuple, idx, fun, std::forward<Args>(args)...)))
        {
            return 0;
        }
    };


	template<typename Tuple, typename F, typename ...Args>
	inline constexpr void visit_at(Tuple &tuple, std::size_t idx, F fun, Args &&...args) noexcept(noexcept(detail::visit_impl<std::tuple_size<Tuple>::value>::visit(tuple, idx, fun, std::forward<Args>(args)...)))
	{
	    detail::visit_impl<std::tuple_size<Tuple>::value>::visit(tuple, idx, fun, std::forward<Args>(args)...);
	}

	template<typename R, typename Tuple, typename F, typename ...Args>
	inline constexpr R visit_at(Tuple &tuple, std::size_t idx, F fun, Args &&...args) noexcept(noexcept(detail::visit_impl<std::tuple_size<Tuple>::value>::template visit<R>(tuple, idx, fun, std::forward<Args>(args)...)))
	{
	    return detail::visit_impl<std::tuple_size<Tuple>::value>::template visit<R>(tuple, idx, fun, std::forward<Args>(args)...);
	}

	// with 2 tuples... only works with tuples of the SAME SIZE
	template <typename Tuple1, typename Tuple2, typename F>
	inline constexpr void visit_at(Tuple1 &tup1, Tuple2 & tup2, std::size_t idx, F fun) //noexcept(noexcept(detail::visit_impl<std::tuple_size<Tuple>::value>::visit(tuple, idx, fun, std::forward<Args>(args)...)))
	{
	    detail::visit_impl<std::tuple_size<Tuple1>::value>::visit(tup1, tup2, idx, fun);
	}

	// // tuple visitor... beautiful implementation
	// // of a tuple visitor to apply some functor to each tuple type
	// template <size_t I>
	// struct visit_impl
	// {
	//     template <typename T, typename F>
	//     static void visit(T& tup, size_t idx, F fun)
	//     {
	//         if (idx == I - 1) fun(std::get<I - 1>(tup));
	//         else visit_impl<I - 1>::visit(tup, idx, fun);
	//     }
	// };

	// template <>
	// struct visit_impl<0>
	// {
	//     template <typename T, typename F>
	//     static void visit(T& tup, size_t idx, F fun) { assert(false); }
	// };

	// template <typename F, typename... Ts>
	// void visit_at(std::tuple<Ts...> const& tup, size_t idx, F fun)
	// {
	//     visit_impl<sizeof...(Ts)>::visit(tup, idx, fun);
	// }

	// template <typename F, typename... Ts>
	// void visit_at(std::tuple<Ts...>& tup, size_t idx, F fun)
	// {
	//     visit_impl<sizeof...(Ts)>::visit(tup, idx, fun);
	// }
}






/** @class VectorStack
 *  @brief VectorStack class to extend vector-type iteration
 *		   to multiple vectors "stacked" one after the other
 *
 *  
 *
 */
template <typename... VectorType>
class VectorStack {
 // : public vector::VectorFunctors<VectorStack<VectorType...>>,
	// 				public vector::VectorAssignment<VectorStack<VectorType...>>{
public:
	typedef VectorStack<VectorType...> 					SelfType;
	// typedef vector::VectorAssignment<SelfType> 			AssignmentType;
	// using AssignmentType::operator=;


private:
	typedef std::tuple<VectorType * ...>				PointerTupleType;
	static constexpr std::size_t mNumVecs 			  = std::tuple_size<PointerTupleType>::value;


	PointerTupleType mVectors;
	std::tuple<VectorType ...> mDerr;
	std::vector<std::size_t> mCumSize; // cumulative size

	template <bool is_const>
	class vs_iterator{
	private:
		typedef typename std::conditional<is_const, 
						const VectorStack, 
						VectorStack>::type 	container_type;
		// typedef typename std::conditional<is_const, 
		// 		typename ContainerType::const_iterator, 
		// 		typename ContainerType::iterator>::type 	iterator_type;

		container_type * mCont;
		// PointerTupleType * mPTup;

		unsigned int mIdx;		// tracker index... starts at 0, goes to size()-1
		unsigned int mContIdx; // index of which container we are at
		std::tuple<typename VectorType::iterator ..., 
				   typename detail::LastTypeOf<VectorType...>::type ::iterator> mIters; // tuple of iterators, plus one extra iterator at the end
	public:
		typedef vs_iterator									self_type;
		typedef typename detail::FirstTypeOf<VectorType...>::type ::iterator::difference_type		difference_type;
		typedef std::remove_reference_t<decltype(*std::get<0>(std::declval<PointerTupleType>())->begin())> 	value_type;
		typedef value_type & 								reference;
		typedef value_type * 								pointer;
		typedef typename detail::FirstTypeOf<VectorType...>::type ::iterator::iterator_category	iterator_category;

		struct tupleganger{
			template <typename T1, typename T2>
			void operator()(T1 & t1, T2 & t2) {t1 = t2->begin();};
		};

		vs_iterator(container_type * cont, unsigned int idx)
		: mCont(cont), mIdx(idx) {
			std::tuple<typename VectorType::iterator ...> its;
			for (unsigned int i=0; i<mNumVecs; i++){
				detail::visit_at(its, mCont->mVectors, i, tupleganger()); //it1 = std::begin(*it2);
			}
			// end it
			std::tuple<typename detail::LastTypeOf<VectorType...>::type ::iterator> endit;
			std::get<0>(endit) = std::get<mNumVecs-1>(mCont->mVectors)->end();

			mIters = std::tuple_cat(its, endit);

			// std::cout << "Idx: " << idx ;
			if (idx == 0){
				mContIdx = 0;
				return;
			}

			for (auto i=0; i<mCont->mCumSize.size(); i++){
				if (mCont->mCumSize[i]/idx >= 1){
					mContIdx = i-1;
					break;
				}
			}
			if (mContIdx == mNumVecs-1) mContIdx = mNumVecs;

			// std::cout << " complete " << std::endl;
		};

		// copy assignment
		vs_iterator & operator=(const vs_iterator & cit){
			vs_iterator i(cit);
			std::swap(i,*this);
			return *this;
		}

		pointer operator->() const {
			return detail::visit_at<pointer>(mIters, mContIdx, [&](auto & it){return it.operator->();});
		};
		reference operator*() const {
			// reference 
			return *detail::visit_at<pointer>(mIters, mContIdx, [&](auto & it){return it.operator->();});
		};

		self_type operator++(){
			mIdx++;
			detail::visit_at(mIters, mContIdx, [&](auto & it){it++;});
			mContIdx += mIdx/mCont->mCumSize[mContIdx];
			return *this;
		};

		self_type operator++(int blah){
			mIdx++;
			detail::visit_at(mIters, mContIdx, [&](auto & it){it++;});
			mContIdx += mIdx/mCont->mCumSize[mContIdx];
			return *this;
		};

		bool operator!=(const self_type & leaf) const {return mIdx != leaf.mIdx;};
		bool operator==(const self_type & leaf) const {return mIdx == leaf.mIdx;};

	};
	
public:
	VectorStack(std::tuple<VectorType * ...> vecs) : mVectors(vecs) {
		mCumSize.resize(mNumVecs);
		detail::visit_at(mVectors, 0, [&](auto & c){mCumSize[0] = c->size();});
		for (auto i=1; i<mNumVecs; i++){
			detail::visit_at(mVectors, i, [&](auto & c){mCumSize[i] = c->size() + mCumSize[i-1];});
		}
	};

	typedef vs_iterator<true> 	const_iterator;
	typedef vs_iterator<false> 	iterator;

	std::size_t size() const {return mCumSize[mNumVecs - 1];};

	iterator begin() {return iterator(this, 0);};
	iterator end()	 {return iterator(this, mCumSize[mNumVecs-1]);};

	// const_iterator cbegin() const {return const_iterator(this, mCont->cbegin());};
	// const_iterator cend() const	 {return const_iterator(this, mCont->cend());};


};




// template <typename ContainerT, typename Functor>
// using VS = libra::VectorStack<ContainerT, Functor>;


template <typename... VectorType>
VectorStack<VectorType...> make_vector_stack(VectorType & ... vectors){
	return VectorStack<VectorType...>(std::make_tuple(&vectors...));
}


} // end namespace libra
#endif