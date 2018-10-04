/** @file Macros.hpp
 *  @brief File with preprocessor Macros
 *
 *
 *  @author D. Pederson
 *  @bug No known bugs. 
 */
 
#ifndef _MACROS_H
#define _MACROS_H

namespace libra{


// concatenate
#define LIBRA_CONCATENATE(arg1, arg2)  arg1##arg2
	
// create a string from the argument
#define LIBRA_STRINGIZE(arg) #arg

// a comma macro for convenience
#define LIBRA_COMMA ,

// a semicolon macro for convenience
#define LIBRA_SEMICOLON ;

// this macro tool will defer the expansion of a macro
#define LIBRA_DEFER(...) __VA_ARGS__


// LIBRA_FOR_EACH takes a macro called "what", and applies it
// to every item passed in the variadic argument list. The argument "op"
// is a macro that is expanded in between each application of "what"
#define LIBRA_FOR_EACH_1(what, op, x, ...) what(x)
#define LIBRA_FOR_EACH_2(what, op, x, ...) what(x) op LIBRA_FOR_EACH_1(what, LIBRA_DEFER(op), __VA_ARGS__)
#define LIBRA_FOR_EACH_3(what, op, x, ...) what(x) op LIBRA_FOR_EACH_2(what, LIBRA_DEFER(op), __VA_ARGS__)
#define LIBRA_FOR_EACH_4(what, op, x, ...) what(x) op LIBRA_FOR_EACH_3(what, LIBRA_DEFER(op), __VA_ARGS__)
#define LIBRA_FOR_EACH_5(what, op, x, ...) what(x) op LIBRA_FOR_EACH_4(what, LIBRA_DEFER(op), __VA_ARGS__)
#define LIBRA_FOR_EACH_6(what, op, x, ...) what(x) op LIBRA_FOR_EACH_5(what, LIBRA_DEFER(op), __VA_ARGS__)
#define LIBRA_FOR_EACH_7(what, op, x, ...) what(x) op LIBRA_FOR_EACH_6(what, LIBRA_DEFER(op), __VA_ARGS__)
#define LIBRA_FOR_EACH_8(what, op, x, ...) what(x) op LIBRA_FOR_EACH_7(what, LIBRA_DEFER(op), __VA_ARGS__)
#define LIBRA_FOR_EACH_9(what, op, x, ...) what(x) op LIBRA_FOR_EACH_8(what, LIBRA_DEFER(op), __VA_ARGS__)
#define LIBRA_FOR_EACH_10(what, op, x, ...)	what(x) op LIBRA_FOR_EACH_9(what, LIBRA_DEFER(op), __VA_ARGS__)
#define LIBRA_FOR_EACH_11(what, op, x, ...)	what(x) op LIBRA_FOR_EACH_10(what, LIBRA_DEFER(op), __VA_ARGS__)
#define LIBRA_FOR_EACH_12(what, op, x, ...)	what(x) op LIBRA_FOR_EACH_11(what, LIBRA_DEFER(op), __VA_ARGS__)
#define LIBRA_FOR_EACH_13(what, op, x, ...)	what(x) op LIBRA_FOR_EACH_12(what, LIBRA_DEFER(op), __VA_ARGS__)
#define LIBRA_FOR_EACH_14(what, op, x, ...) what(x) op LIBRA_FOR_EACH_13(what, LIBRA_DEFER(op), __VA_ARGS__)
#define LIBRA_FOR_EACH_15(what, op, x, ...)	what(x) op LIBRA_FOR_EACH_14(what, LIBRA_DEFER(op), __VA_ARGS__)
#define LIBRA_FOR_EACH_16(what, op, x, ...)	what(x) op LIBRA_FOR_EACH_15(what, LIBRA_DEFER(op), __VA_ARGS__)
#define LIBRA_FOR_EACH_17(what, op, x, ...)	what(x) op LIBRA_FOR_EACH_16(what, LIBRA_DEFER(op), __VA_ARGS__)
#define LIBRA_FOR_EACH_18(what, op, x, ...)	what(x) op LIBRA_FOR_EACH_17(what, LIBRA_DEFER(op), __VA_ARGS__)
#define LIBRA_FOR_EACH_19(what, op, x, ...) what(x) op LIBRA_FOR_EACH_18(what, LIBRA_DEFER(op), __VA_ARGS__)
#define LIBRA_FOR_EACH_20(what, op, x, ...)	what(x) op LIBRA_FOR_EACH_19(what, LIBRA_DEFER(op), __VA_ARGS__)
#define LIBRA_FOR_EACH_21(what, op, x, ...)	what(x) op LIBRA_FOR_EACH_20(what, LIBRA_DEFER(op), __VA_ARGS__)
#define LIBRA_FOR_EACH_22(what, op, x, ...)	what(x) op LIBRA_FOR_EACH_21(what, LIBRA_DEFER(op), __VA_ARGS__)
#define LIBRA_FOR_EACH_23(what, op, x, ...)	what(x) op LIBRA_FOR_EACH_22(what, LIBRA_DEFER(op), __VA_ARGS__)
#define LIBRA_FOR_EACH_24(what, op, x, ...) what(x) op LIBRA_FOR_EACH_23(what, LIBRA_DEFER(op), __VA_ARGS__)
#define LIBRA_FOR_EACH_25(what, op, x, ...)	what(x) op LIBRA_FOR_EACH_24(what, LIBRA_DEFER(op), __VA_ARGS__)
#define LIBRA_FOR_EACH_26(what, op, x, ...)	what(x) op LIBRA_FOR_EACH_25(what, LIBRA_DEFER(op), __VA_ARGS__)
#define LIBRA_FOR_EACH_27(what, op, x, ...)	what(x) op LIBRA_FOR_EACH_26(what, LIBRA_DEFER(op), __VA_ARGS__)
#define LIBRA_FOR_EACH_28(what, op, x, ...)	what(x) op LIBRA_FOR_EACH_27(what, LIBRA_DEFER(op), __VA_ARGS__)
#define LIBRA_FOR_EACH_29(what, op, x, ...) what(x) op LIBRA_FOR_EACH_28(what, LIBRA_DEFER(op), __VA_ARGS__)
#define LIBRA_FOR_EACH_30(what, op, x, ...)	what(x) op LIBRA_FOR_EACH_29(what, LIBRA_DEFER(op), __VA_ARGS__)
#define LIBRA_FOR_EACH_31(what, op, x, ...)	what(x) op LIBRA_FOR_EACH_30(what, LIBRA_DEFER(op), __VA_ARGS__)
#define LIBRA_FOR_EACH_32(what, op, x, ...)	what(x) op LIBRA_FOR_EACH_31(what, LIBRA_DEFER(op), __VA_ARGS__)
#define LIBRA_FOR_EACH_33(what, op, x, ...)	what(x) op LIBRA_FOR_EACH_32(what, LIBRA_DEFER(op), __VA_ARGS__)
#define LIBRA_FOR_EACH_34(what, op, x, ...) what(x) op LIBRA_FOR_EACH_33(what, LIBRA_DEFER(op), __VA_ARGS__)
#define LIBRA_FOR_EACH_35(what, op, x, ...)	what(x) op LIBRA_FOR_EACH_34(what, LIBRA_DEFER(op), __VA_ARGS__)
#define LIBRA_FOR_EACH_36(what, op, x, ...)	what(x) op LIBRA_FOR_EACH_35(what, LIBRA_DEFER(op), __VA_ARGS__)
#define LIBRA_FOR_EACH_37(what, op, x, ...)	what(x) op LIBRA_FOR_EACH_36(what, LIBRA_DEFER(op), __VA_ARGS__)
#define LIBRA_FOR_EACH_38(what, op, x, ...)	what(x) op LIBRA_FOR_EACH_37(what, LIBRA_DEFER(op), __VA_ARGS__)
#define LIBRA_FOR_EACH_39(what, op, x, ...) what(x) op LIBRA_FOR_EACH_38(what, LIBRA_DEFER(op), __VA_ARGS__)
#define LIBRA_FOR_EACH_40(what, op, x, ...)	what(x) op LIBRA_FOR_EACH_39(what, LIBRA_DEFER(op), __VA_ARGS__)
#define LIBRA_FOR_EACH_41(what, op, x, ...)	what(x) op LIBRA_FOR_EACH_40(what, LIBRA_DEFER(op), __VA_ARGS__)
#define LIBRA_FOR_EACH_42(what, op, x, ...)	what(x) op LIBRA_FOR_EACH_41(what, LIBRA_DEFER(op), __VA_ARGS__)
#define LIBRA_FOR_EACH_43(what, op, x, ...)	what(x) op LIBRA_FOR_EACH_42(what, LIBRA_DEFER(op), __VA_ARGS__)
#define LIBRA_FOR_EACH_44(what, op, x, ...) what(x) op LIBRA_FOR_EACH_43(what, LIBRA_DEFER(op), __VA_ARGS__)
#define LIBRA_FOR_EACH_45(what, op, x, ...)	what(x) op LIBRA_FOR_EACH_44(what, LIBRA_DEFER(op), __VA_ARGS__)
#define LIBRA_FOR_EACH_46(what, op, x, ...)	what(x) op LIBRA_FOR_EACH_45(what, LIBRA_DEFER(op), __VA_ARGS__)
#define LIBRA_FOR_EACH_47(what, op, x, ...)	what(x) op LIBRA_FOR_EACH_46(what, LIBRA_DEFER(op), __VA_ARGS__)
#define LIBRA_FOR_EACH_48(what, op, x, ...)	what(x) op LIBRA_FOR_EACH_47(what, LIBRA_DEFER(op), __VA_ARGS__)
#define LIBRA_FOR_EACH_49(what, op, x, ...) what(x) op LIBRA_FOR_EACH_48(what, LIBRA_DEFER(op), __VA_ARGS__)
#define LIBRA_FOR_EACH_50(what, op, x, ...)	what(x) op LIBRA_FOR_EACH_49(what, LIBRA_DEFER(op), __VA_ARGS__)
#define LIBRA_FOR_EACH_51(what, op, x, ...)	what(x) op LIBRA_FOR_EACH_50(what, LIBRA_DEFER(op), __VA_ARGS__)
#define LIBRA_FOR_EACH_52(what, op, x, ...)	what(x) op LIBRA_FOR_EACH_51(what, LIBRA_DEFER(op), __VA_ARGS__)
#define LIBRA_FOR_EACH_53(what, op, x, ...)	what(x) op LIBRA_FOR_EACH_52(what, LIBRA_DEFER(op), __VA_ARGS__)
#define LIBRA_FOR_EACH_54(what, op, x, ...) what(x) op LIBRA_FOR_EACH_53(what, LIBRA_DEFER(op), __VA_ARGS__)
#define LIBRA_FOR_EACH_55(what, op, x, ...)	what(x) op LIBRA_FOR_EACH_54(what, LIBRA_DEFER(op), __VA_ARGS__)
#define LIBRA_FOR_EACH_56(what, op, x, ...)	what(x) op LIBRA_FOR_EACH_55(what, LIBRA_DEFER(op), __VA_ARGS__)
#define LIBRA_FOR_EACH_57(what, op, x, ...)	what(x) op LIBRA_FOR_EACH_56(what, LIBRA_DEFER(op), __VA_ARGS__)
#define LIBRA_FOR_EACH_58(what, op, x, ...)	what(x) op LIBRA_FOR_EACH_57(what, LIBRA_DEFER(op), __VA_ARGS__)
#define LIBRA_FOR_EACH_59(what, op, x, ...) what(x) op LIBRA_FOR_EACH_58(what, LIBRA_DEFER(op), __VA_ARGS__)

#define LIBRA_FOR_EACH_NARG(...) LIBRA_FOR_EACH_NARG_(__VA_ARGS__, LIBRA_FOR_EACH_RSEQ_N())
#define LIBRA_FOR_EACH_NARG_(...) LIBRA_FOR_EACH_ARG_N(__VA_ARGS__) 
#define LIBRA_FOR_EACH_ARG_N(_1, _2, _3, _4, _5, _6, _7, _8, _9, _10, _11, _12, _13, _14, _15, _16, _17, _18, _19, _20, _21, _22, _23, _24, _25, _26, _27, _28, _29, _30, _31, _32, _33, _34, _35, _36, _37, _38, _39, _40, _41, _42, _43, _44, _45, _46, _47, _48, _49, _50, _51, _52, _53, _54, _55, _56, _57, _58, _59, N, ...) N 
#define LIBRA_FOR_EACH_RSEQ_N() 59, 58, 57, 56, 55, 54, 53, 52, 51, 50, 49, 48, 47, 46, 45, 44, 43, 42, 41, 40, 39, 38, 37, 36, 35, 34, 33, 32, 31, 30, 29, 28, 27, 26, 25, 24, 23, 22, 21, 20, 19, 18, 17, 16, 15, 14, 13, 12, 11, 10, 9, 8, 7, 6, 5, 4, 3, 2, 1, 0

#define LIBRA_FOR_EACH_(N, what, op, x, ...) LIBRA_CONCATENATE(LIBRA_FOR_EACH_, N)(what, LIBRA_DEFER(op), x, __VA_ARGS__)
#define LIBRA_FOR_EACH(what, op, x, ...) LIBRA_FOR_EACH_(LIBRA_FOR_EACH_NARG(x, __VA_ARGS__), what, LIBRA_DEFER(op), x, __VA_ARGS__)


// for each that creates a sequence... just inserts commas between items
#define LIBRA_FOR_EACH_SEQ(what, x, ...) LIBRA_FOR_EACH(what, LIBRA_DEFER(LIBRA_COMMA), x, __VA_ARGS__)

// for each that creates seperation... just inserts semicolons between items
#define LIBRA_FOR_EACH_SEP(what, x, ...) LIBRA_FOR_EACH(what, LIBRA_DEFER(LIBRA_SEMICOLON), x, __VA_ARGS__)

// define a function macro overload on number of arguments (up to 2)
#define LIBRA_GET_MACRO(_1, _2, NAME, ...) NAME


struct DefaultNullInterfacePolicy{
	template <typename iterator_reference>
	static decltype(auto) get(iterator_reference & it) {
		return std::forward<iterator_reference>(it);
	}
};



// helper dummy type
template <typename InterfacePolicy, typename... Args>
struct LibraHasMethodHelper{};

// macro for defining the struct LibraHasMethod##MethodName
#define LIBRA_HAS_METHOD_DEF(MethodName) 						\
	template<typename T, typename _ = void> 					\
	struct LibraHasMethod##MethodName : std::false_type {};	\
																\
	template <typename T>										\
	struct LibraHasMethod##MethodName <						\
	        T,													\
	        std::conditional_t<									\
	            false,											\
	            libra::LibraHasMethodHelper<					\
	                decltype(std::declval<T>().MethodName())	\
	                >,											\
	            void											\
	            >												\
	        > : public std::true_type {};						

// macro for instantiating the struct LibraHasMethod##MethodName
#define LIBRA_HAS_METHOD(MethodName) LibraHasMethod##MethodName









// macro for name of the struct LibraFunctorFor##MethodName
#define LIBRA_FUNCTOR_FOR(FunctionName) LibraFunctorFor##FunctionName
#define LIBRA_INSTANTIATE_FUNCTOR_FOR(Name) LIBRA_FUNCTOR_FOR(Name)()
// macro for defining the struct LibraFunctorFor##FunctionName
#define LIBRA_FUNCTOR_FOR_DEF_NOINTERFACE(FunctionName) 									\
	struct LIBRA_FUNCTOR_FOR(FunctionName) { 									\
		template <typename Iterator,											\
				  typename T = 													\
				  typename std::enable_if<										\
				  		   LIBRA_HAS_METHOD(FunctionName)<decltype(*std::declval<Iterator>())>::value, 		\
				  		   Iterator 											\
				  		   >::type 												\
				  		   > 													\
		static decltype(auto) get(Iterator & it){return it->FunctionName();}; 	\
																				\
		template <typename Iterator>											\
		decltype(auto) operator()(Iterator & it){return get(it);};				\
	};





// macro for defining the struct LibraFunctorFor##FunctionName
#define LIBRA_FUNCTOR_FOR_DEF_INTERFACE(FunctionName) 									\
	struct LIBRA_FUNCTOR_FOR(FunctionName) { 							 		\
		template <typename Iterator,											\
				  typename T = 													\
				  typename std::enable_if<										\
				  		   LIBRA_HAS_METHOD(FunctionName)<decltype(InterfacePolicy::get(*std::declval<Iterator>()))>::value, 		\
				  		   Iterator 											\
				  		   >::type 												\
				  		   > 													\
		static decltype(auto) get(Iterator & it){return InterfacePolicy::get(*it).FunctionName();}; 	\
																				\
		template <typename Iterator>											\
		decltype(auto) operator()(Iterator & it){return get(it);};				\
	};


// define an overloaded version of FUNCTOR_FOR_DEF
#define LIBRA_FUNCTOR_FOR_DEF(...) LIBRA_GET_MACRO(__VA_ARGS__, LIBRA_FUNCTOR_FOR_DEF_INTERFACE, LIBRA_FUNCTOR_FOR_DEF_NOINTERFACE)(__VA_ARGS__)


					

#define LIBRA_FUNCTOR_PREPARE(FunctionName, ...)	\
	LIBRA_HAS_METHOD_DEF(FunctionName);			\
	LIBRA_FUNCTOR_FOR_DEF_NOINTERFACE(FunctionName);			


#define LIBRA_FUNCTOR_PREPARE_INTERFACE(FunctionName)	\
	LIBRA_HAS_METHOD_DEF(FunctionName);			\
	LIBRA_FUNCTOR_FOR_DEF_INTERFACE(FunctionName);	


} // end namespace libra
#endif