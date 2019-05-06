#include <libra.h>

using namespace std;
using namespace libra;

// g++ -std=c++14 -O2 -I../ expressions_test.cpp -o expr_test

#include <typeinfo>
#include <functional>


struct Data{
private:
	std::array<double, 5> mA;
public:
	template <std::size_t I>
	double & A(){return mA[I];};
};

template <std::size_t I>
struct GetA{
	template <typename T>
	static double & get(T && o){return o.template A<I>();};
};

template <typename VecType>
struct Print{
	template <typename ... T>
	struct type{
		template <typename ... Args>
		static void get(Args && ... args){
			std::cout << VecType::template get<T::value...>(std::forward<Args>(args)...) << ", " ;
			// if (T2::value == 0) std::cout << std::endl;
		};
	};
};

int main(int argc, char * argv[]){


	//**************** GENERALIZED VECTOR TESTS *********************//
	cout << "********************************************" << endl;
	cout << "********************************************" << endl;
	cout << "******** EXPRESSION VECTOR TESTS ***********" << endl;
	cout << "********************************************" << endl;
	cout << "********************************************" << endl;
	cout << "********************************************" << endl;
	
	Data d;
	d.A<0>() = 0;
	d.A<1>() = 1;
	d.A<2>() = 2;
	d.A<3>() = 3;
	d.A<4>() = 4;

	typedef ExpressionVector<GetA<0>, GetA<1>, GetA<2>, GetA<3>, GetA<4>> 	fwd_vec;
	typedef ExpressionVector<GetA<4>, GetA<3>, GetA<2>, GetA<1>, GetA<0>> 	bwd_vec;
	typedef ExpressionVector<GetA<0>, GetA<2>, GetA<1>, GetA<4>, GetA<3>> 	swizzle_vec;

	std::cout << "Forward: " ;
	nested_for_each_tuple_type<Print<fwd_vec>::type, Detail::RangeTuple<0, 4>>(d);
	std::cout << std::endl;

	std::cout << "Forward (Reordered): " ;
	nested_for_each_tuple_type<Print<fwd_vec>::type, Detail::IndexTuple<4, 3, 2, 1, 0>>(d);
	std::cout << std::endl;

	std::cout << "Backward: " ;
	nested_for_each_tuple_type<Print<bwd_vec>::type, Detail::RangeTuple<0, 4>>(d);
	std::cout << std::endl;

	std::cout << "Swizzle: " ;
	nested_for_each_tuple_type<Print<swizzle_vec>::type, Detail::RangeTuple<0, 4>>(d);
	std::cout << std::endl;

	typedef ExpressionSum<fwd_vec, bwd_vec> 	fwd_plus_bwd_vec;

	std::cout << "Fwd + Bwd: " ;
	nested_for_each_tuple_type<Print<fwd_plus_bwd_vec>::type, Detail::RangeTuple<0, 4>>(d);
	std::cout << std::endl;

	typedef ExpressionElementwiseProduct<fwd_vec, bwd_vec> 	fwd_times_bwd_vec;

	std::cout << "Fwd * Bwd: " ;
	nested_for_each_tuple_type<Print<fwd_times_bwd_vec>::type, Detail::RangeTuple<0, 4>>(d);
	std::cout << std::endl;

	std::cout << "Sum(Fwd * Bwd): " << ExpressionSumElements<fwd_times_bwd_vec>::get(d) << std::endl;
	std::cout << "InnerProduct(Fwd, Swizzle): " << ExpressionInnerProduct<fwd_vec, swizzle_vec>::get(d) << std::endl;

	typedef ExpressionMatrix<3, 3, 
							 GetA<0>, GetA<0>, GetA<0>,
							 GetA<1>, GetA<1>, GetA<1>, 
							 GetA<2>, GetA<2>, GetA<2>> 	two_mtx;
	std::cout << "Matrix: " << std::endl;
	nested_for_each_tuple_type<Print<two_mtx>::type, Detail::RangeTuple<0,two_mtx::rows-1>, Detail::RangeTuple<0,two_mtx::cols-1>>(d);
	std::cout << std::endl;

	typedef ExpressionMatrixTranspose<two_mtx> 	two_trans;
	std::cout << "Matrix Transpose: " << std::endl;
	nested_for_each_tuple_type<Print<two_trans>::type, Detail::RangeTuple<0,two_mtx::rows-1>, Detail::RangeTuple<0,two_mtx::cols-1>>(d);
	std::cout << std::endl;

	typedef ExpressionMatrixRow<0, two_mtx> 	two_row;
	std::cout << "Matrix Row 0: " << std::endl;
	nested_for_each_tuple_type<Print<two_row>::type, Detail::RangeTuple<0,two_row::size-1>>(d);
	std::cout << std::endl;

	typedef ExpressionMatrixCol<0, two_mtx> 	two_col;
	std::cout << "Matrix Col 0: " << std::endl;
	nested_for_each_tuple_type<Print<two_col>::type, Detail::RangeTuple<0,two_col::size-1>>(d);
	std::cout << std::endl;

	typedef ExpressionVector<GetA<0>, GetA<1>, GetA<2>> 	two_vec;
	typedef ExpressionMatrixVectorProduct<two_mtx, two_vec> 	two_matvec;
	std::cout << "Matrix-Vector multiply: " << std::endl;
	nested_for_each_tuple_type<Print<two_matvec>::type, Detail::RangeTuple<0,two_matvec::size-1>>(d);
	std::cout << std::endl;

	typedef ExpressionMatrixTrace<two_mtx> 	two_trace;
	std::cout << "Matrix Trace: " << two_trace::get(d) << std::endl;

	typedef ExpressionMatrixDiagonal<0, two_mtx> 	two_diag;
	std::cout << "Matrix diagonal: " << std::endl;
	nested_for_each_tuple_type<Print<two_diag>::type, Detail::RangeTuple<0,two_diag::size-1>>(d);
	std::cout << std::endl;

	typedef ExpressionMatrixDiagonal<-1, two_mtx> 	two_diag_m1;
	std::cout << "Matrix diagonal: " << std::endl;
	nested_for_each_tuple_type<Print<two_diag_m1>::type, Detail::RangeTuple<0,two_diag_m1::size-1>>(d);
	std::cout << std::endl;

	typedef ExpressionMatrixDiagonal<-2, two_mtx> 	two_diag_m2;
	std::cout << "Matrix diagonal: " << std::endl;
	nested_for_each_tuple_type<Print<two_diag_m2>::type, Detail::RangeTuple<0,two_diag_m2::size-1>>(d);
	std::cout << std::endl;

	typedef ExpressionMatrixDiagonal<2, two_mtx> 	two_diag_p2;
	std::cout << "Matrix diagonal: " << std::endl;
	nested_for_each_tuple_type<Print<two_diag_p2>::type, Detail::RangeTuple<0,two_diag_p2::size-1>>(d);
	std::cout << std::endl;

	typedef ExpressionMatrixSymmetricPart<two_mtx> 	two_symm;
	std::cout << "Matrix symmetric part: " << std::endl;
	nested_for_each_tuple_type<Print<two_symm>::type, Detail::RangeTuple<0,two_symm::rows-1>, Detail::RangeTuple<0,two_symm::cols-1>>(d);
	std::cout << std::endl;

	typedef ExpressionMatrixAntisymmetricPart<two_mtx> 	two_asymm;
	std::cout << "Matrix antisymmetric part: " << std::endl;
	nested_for_each_tuple_type<Print<two_asymm>::type, Detail::RangeTuple<0,two_asymm::rows-1>, Detail::RangeTuple<0,two_asymm::cols-1>>(d);
	std::cout << std::endl;

	typedef ExpressionMatrixSum<two_symm, two_asymm> 	two_sum;
	std::cout << "Matrix sum: " << std::endl;
	nested_for_each_tuple_type<Print<two_sum>::type, Detail::RangeTuple<0,two_sum::rows-1>, Detail::RangeTuple<0,two_sum::cols-1>>(d);
	std::cout << std::endl;

	typedef ExpressionMatrixMatrixProduct<two_mtx, two_sum> 	two_prod;
	std::cout << "Matrix-Matrix Product: " << std::endl;
	nested_for_each_tuple_type<Print<two_prod>::type, Detail::RangeTuple<0,two_prod::rows-1>, Detail::RangeTuple<0,two_prod::cols-1>>(d);
	std::cout << std::endl;

	typedef ExpressionMatrixConcatVert<two_mtx, two_sum> 	two_vertcat;
	std::cout << "Matrix-Matrix Concat Vert: " << std::endl;
	nested_for_each_tuple_type<Print<two_vertcat>::type, Detail::RangeTuple<0,two_vertcat::rows-1>, Detail::RangeTuple<0,two_vertcat::cols-1>>(d);
	std::cout << std::endl;

	typedef ExpressionMatrixConcatHorz<two_mtx, two_sum> 	two_cathorz;
	std::cout << "Matrix-Matrix Concat Horz: " << std::endl;
	nested_for_each_tuple_type<Print<two_cathorz>::type, Detail::RangeTuple<0,two_cathorz::rows-1>, Detail::RangeTuple<0,two_cathorz::cols-1>>(d);
	std::cout << std::endl;

	std::cout << "Matrix: " << std::endl;
	nested_for_each_tuple_type<Print<two_mtx>::type, Detail::RangeTuple<0,two_mtx::rows-1>, Detail::RangeTuple<0,two_mtx::cols-1>>(d);
	std::cout << std::endl;
	return 0;
}