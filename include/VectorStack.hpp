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
	template <typename T, typename... Ts>
	struct LastTypeOf {
	   typedef typename LastTypeOf<Ts...>::type type;
	};

	template <typename T>
	struct LastTypeOf<T> {
	  typedef T type;
	};



	template <typename T, typename... Ts>
	struct FirstTypeOf {
	   typedef T type;
	};

	template <typename T>
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
        	// std::cout << "idx = " << idx << " : Visiting Tuple Index =" << I << std::endl;
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
        	// std::cout << "idx = " << idx << " : Visiting Tuple Index =" << 0 << std::endl;
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
class VectorStack : public vector::VectorFunctors<VectorStack<VectorType...>>,
					public vector::VectorAssignment<VectorStack<VectorType...>>{
public:
	typedef VectorStack<VectorType...> 					SelfType;
	typedef vector::VectorAssignment<SelfType> 			AssignmentType;
	using AssignmentType::operator=;


private:
	typedef std::tuple<VectorType * ...>				PointerTupleType;
	static constexpr std::size_t mNumVecs 			  = std::tuple_size<PointerTupleType>::value;


	PointerTupleType mVectors;
	std::vector<std::size_t> mCumSize; // cumulative size

	template <bool is_const>
	class vs_iterator{
	private:
		typedef typename std::conditional<is_const, 
						const VectorStack, 
						VectorStack>::type 	container_type;

		typedef typename std::conditional<is_const, 
				typename detail::FirstTypeOf<VectorType...>::type ::const_iterator, 
				typename detail::FirstTypeOf<VectorType...>::type ::iterator>::type 	iterator_type;

		typedef typename std::conditional<is_const, 
				std::tuple<typename VectorType::const_iterator ...>,
				std::tuple<typename VectorType::iterator ...>
				>::type 			first_iterator_tuple_types;

		typedef typename std::conditional<is_const, 
				std::tuple<typename detail::LastTypeOf<VectorType...>::type ::const_iterator>,
				std::tuple<typename detail::LastTypeOf<VectorType...>::type ::iterator>
				>::type 			last_iterator_tuple_type;


		typedef typename std::conditional<is_const, 
				std::tuple<typename VectorType::const_iterator ..., typename detail::LastTypeOf<VectorType...>::type ::const_iterator>,
				std::tuple<typename VectorType::iterator ..., typename detail::LastTypeOf<VectorType...>::type ::iterator>
				>::type 			iterator_tuple_type;


		container_type * mCont;
		unsigned int mIdx;		// tracker index... starts at 0, goes to size()-1
		unsigned int mContIdx; // index of which container we are at
		iterator_tuple_type mIters; // tuple of iterators, plus one extra iterator at the end
	

		struct tupleganger{
			template <typename T1, typename T2, bool C = is_const>
			typename std::enable_if<C==false, void>::type 
			operator()(T1 & t1, T2 & t2) {
				t1 = t2->begin();
			};

			template <typename T1, typename T2, bool C = is_const>
			typename std::enable_if<C==true, void>::type 
			operator()(T1 & t1, T2 & t2) {
				t1 = t2->cbegin();
			};
		};

		struct tupleganger_end{
			template <typename T1, typename T2, bool C = is_const>
			static typename std::enable_if<C==false, void>::type 
			do_it(T1 & t1, T2 & t2) {
				t1 = t2->end();
			};

			template <typename T1, typename T2, bool C = is_const>
			static typename std::enable_if<C==true, void>::type 
			do_it(T1 & t1, T2 & t2) {
				t1 = t2->cend();
			};
		};

	public:
		typedef vs_iterator									self_type;
		typedef typename iterator_type::difference_type		difference_type;
		typedef typename std::conditional<is_const, 
						 std::remove_reference_t<decltype(*std::get<0>(std::declval<PointerTupleType>())->cbegin())>,
						 std::remove_reference_t<decltype(*std::get<0>(std::declval<PointerTupleType>())->begin())>>::type
						 value_type;
		typedef typename std::conditional<is_const, 
						 const value_type,
						 value_type>::type &				reference;
		typedef typename std::conditional<is_const, 
						 const value_type,
						 value_type>::type * 				pointer;
		typedef typename iterator_type::iterator_category	iterator_category;



		vs_iterator(container_type * cont, unsigned int idx)
		: mCont(cont), mIdx(idx) {

				first_iterator_tuple_types its;// = std::make_tuple(mCont);
				detail::visit_at(its, mCont->mVectors, 1, tupleganger());
				for (unsigned int i=0; i<mNumVecs; i++){
					detail::visit_at(its, mCont->mVectors, i, tupleganger()); //it1 = std::begin(*it2);
				}

				// end it
				last_iterator_tuple_type endit;
				tupleganger_end::do_it(std::get<0>(endit), std::get<mNumVecs-1>(mCont->mVectors));

				mIters = std::tuple_cat(its, endit);

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
		};

		// copy constructor
		vs_iterator(const vs_iterator & v)
		: mCont(v.mCont), mIdx(v.mIdx), mContIdx(v.mContIdx), mIters(v.mIters) {};


		// move constructor
		vs_iterator(vs_iterator && v)
		: mCont(v.mCont), mIdx(v.mIdx), mContIdx(v.mContIdx), mIters(v.mIters) {};

		// copy assignment
		vs_iterator & operator=(const vs_iterator & cit){
			if (cit != *this){
				vs_iterator i(cit);
				std::swap(i.mCont, mCont);
				std::swap(i.mIdx, mIdx);
				std::swap(i.mContIdx, mContIdx);
				std::swap(i.mIters, mIters);
			}
			return *this;
		}

		// move assignment
		vs_iterator & operator=(vs_iterator && cit){
			if (cit != *this){
				vs_iterator i(std::move(cit));
				std::swap(i.mCont, mCont);
				std::swap(i.mIdx, mIdx);
				std::swap(i.mContIdx, mContIdx);
				std::swap(i.mIters, mIters);
			}
			return *this;
		}

		pointer operator->() {
			return detail::visit_at<pointer>(mIters, mContIdx, [&](auto & it){return it.operator->();});
		};
		reference operator*() {
			// reference 
			return *detail::visit_at<pointer>(mIters, mContIdx, [&](auto & it){return it.operator->();});
		};


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

	const_iterator cbegin() const {return const_iterator(this, 0);};
	const_iterator cend() const	 {return const_iterator(this, mCumSize[mNumVecs-1]);};


	template <std::size_t I>
	typename std::tuple_element<I, std::tuple<VectorType...>>::type & 
	get() {
		return *std::get<I>(mVectors);
	}

	template <std::size_t I>
	const typename std::tuple_element<I, std::tuple<VectorType...>>::type & 
	get() const {
		return *std::get<I>(mVectors);
	}

};




template <typename... VectorType>
VectorStack<VectorType...> make_vector_stack(VectorType & ... vectors){
	return VectorStack<VectorType...>(std::make_tuple(&vectors...));
}


// implementation with rvalue references as arguments
template <typename... VectorType>
VectorStack<VectorType...> make_vector_stack(VectorType && ... vectors){
	return VectorStack<VectorType...>(std::make_tuple(new VectorType(vectors)...));
}


} // end namespace libra
#endif