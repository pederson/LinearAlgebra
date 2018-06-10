#ifndef _VECTORVIEW_H
#define _VECTORVIEW_H

#include <stdio.h>
#include <stdlib.h>
#include <complex>
#include <algorithm>
#include <set>
#include <iterator>


namespace libra{




template <typename VectorT, typename IteratorT>
class VectorView{
public:
	typedef typename VectorT::const_iterator 		vector_const_iterator;


	VectorView(VectorT & v, IteratorT beg, IteratorT end)
	: mVec(v), mBegin(beg), mEnd(end) {};


	size_type size() const {return mEnd - mBegin;};


	class const_iterator;
	class iterator{
	public:
		friend class const_iterator;
		typedef iterator 							self_type;
		typedef typename IteratorT::difference_type difference_type;
	    typedef typename IteratorT::value_type 		value_type;
	    typedef typename IteratorT::reference 		reference;
	    typedef typename IteratorT::pointer 		pointer;
	    typedef typename IteratorT::iterator_category	iterator_category;

		// construction
		iterator(VectorView & sv, IteratorT it)
		: mSV(sv)
		, mIt(it){};

		// copy assignment
		iterator & operator=(const iterator & cit){
			iterator i(cit);
			std::swap(i,*this);
			return *this;
		}

		// rvalue dereferencing
		pointer operator->() {return mIt.operator->();};
		reference operator*(){ return *mIt;};

		// increment operators
		self_type operator++(){
			mIt++;
			return *this;
		}
		self_type operator++(int blah){
			mIt++;
			return *this;
		}

		// decrement operators
		self_type operator--(){
			mIt--;
			return *this;
		}
		self_type operator--(int blah){
			mIt--;
			return *this;
		}

		// scalar arithmetic operators
		self_type operator+(int n){
			mIt = mIt + n;
			return *this;
		}
		self_type operator-(int n){
			mIt = mIt - n;
			return *this;
		}
		int operator-(const self_type & b) const {
			return mIt - b.mIt;
		}

		// equivalence operators
		bool operator!=(const self_type & leaf) const {return mIt != leaf.mIt;};
		bool operator==(const self_type & leaf) const {return mIt == leaf.mIt;};

		// relational operators
		bool operator>(const self_type & leaf) const {return mIt > leaf.mIt;};
		bool operator>=(const self_type & leaf) const {return mIt >= leaf.mIt;};
		bool operator<(const self_type & leaf) const {return mIt < leaf.mIt;};
		bool operator<=(const self_type & leaf) const {return mIt <= leaf.mIt;};


		// compound assignment operators
		self_type operator+=(int n){
			mIt += n;
			return *this;
		}
		self_type operator-=(int n){
			mIt -= n;
			return *this;
		}


		// offset dereference operator
		reference operator[](int n){
			return mIt[n];
		}
	private:
		VectorView & mSV;
		IteratorT mIt;
	};

	iterator begin() {return iterator(*this, mBegin);};
	iterator end()	 {return iterator(*this, mEnd);};




	class const_iterator{
	public:
		typedef typename VectorT::const_iterator 		subiterator_type;
		
		typedef const_iterator 							self_type;
		typedef typename subiterator_type::difference_type		difference_type;
	    typedef typename subiterator_type::value_type 	value_type;
	    typedef typename subiterator_type::reference 	reference;
	    typedef typename subiterator_type::pointer 		pointer;
	    typedef typename subiterator_type::iterator_category	iterator_category;

		// construction
		const_iterator(const VectorView & sv, subiterator_type it)
		: mSV(sv)
		, mIt(it){};

		// conversion of iterator to const_iterator
		const_iterator(const iterator & it)
		: mSV(it.mSV)
		, mIt(it.mIt) {};


		const_iterator & operator=(const const_iterator & cit){
			const_iterator i(cit);
			std::swap(i,*this);
			return *this;
		}

		// rvalue dereferencing
		pointer operator->() {return mIt.operator->();};
		reference operator*(){ return *mIt;};

		// increment operators
		self_type operator++(){
			mIt++;
			return *this;
		}
		self_type operator++(int blah){
			mIt++;
			return *this;
		}

		// decrement operators
		self_type operator--(){
			mIt--;
			return *this;
		}
		self_type operator--(int blah){
			mIt--;
			return *this;
		}

		// scalar arithmetic operators
		self_type operator+(int n){
			mIt = mIt + n;
			return *this;
		}
		self_type operator-(int n){
			mIt = mIt - n;
			return *this;
		}
		int operator-(const self_type & b) const {
			return mIt - b.mIt;
		}

		// equivalence operators
		bool operator!=(const self_type & leaf) const {return mIt != leaf.mIt;};
		bool operator==(const self_type & leaf) const {return mIt == leaf.mIt;};

		// relational operators
		bool operator>(const self_type & leaf) const {return mIt > leaf.mIt;};
		bool operator>=(const self_type & leaf) const {return mIt >= leaf.mIt;};
		bool operator<(const self_type & leaf) const {return mIt < leaf.mIt;};
		bool operator<=(const self_type & leaf) const {return mIt <= leaf.mIt;};


		// compound assignment operators
		self_type operator+=(int n){
			mIt += n;
			return *this;
		}
		self_type operator-=(int n){
			mIt -= n;
			return *this;
		}


		// offset dereference operator
		reference operator[](int n){
			return mIt[n];
		}
	private:
		const VectorView & mSV;
		subiterator_type mIt;
	};

	const_iterator cbegin() const {return const_iterator(*this, vector_const_iterator(mBegin));};
	const_iterator cend() const	 {return const_iterator(*this, vector_const_iterator(mEnd));};

	// these are required by std::cbegin()/cend()
	const_iterator begin() const {return cbegin();};
	const_iterator end() const	 {return cend();};

protected:
	VectorT & mVec;
	IteratorT mBegin, mEnd;
};



template <typename VectorT, typename IteratorT>
VectorView<VectorT, IteratorT> vector_view(VectorT & vec, IteratorT beg, IteratorT end){
	return VectorView<VectorT, IteratorT>(vec, beg, end);
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