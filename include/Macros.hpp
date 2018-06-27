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


#define LIBRA_FOR_EACH_NARG(...) LIBRA_FOR_EACH_NARG_(__VA_ARGS__, LIBRA_FOR_EACH_RSEQ_N())
#define LIBRA_FOR_EACH_NARG_(...) LIBRA_FOR_EACH_ARG_N(__VA_ARGS__) 
#define LIBRA_FOR_EACH_ARG_N(_1, _2, _3, _4, _5, _6, _7, _8, _9, _10, _11, _12, _13, _14, _15, _16, _17, _18, N, ...) N 
#define LIBRA_FOR_EACH_RSEQ_N() 18, 17, 16, 15, 14, 13, 12, 11, 10, 9, 8, 7, 6, 5, 4, 3, 2, 1, 0

#define LIBRA_FOR_EACH_(N, what, op, x, ...) LIBRA_CONCATENATE(LIBRA_FOR_EACH_, N)(what, LIBRA_DEFER(op), x, __VA_ARGS__)
#define LIBRA_FOR_EACH(what, op, x, ...) LIBRA_FOR_EACH_(LIBRA_FOR_EACH_NARG(x, __VA_ARGS__), what, LIBRA_DEFER(op), x, __VA_ARGS__)


// for each that creates a sequence... just inserts commas between items
#define LIBRA_FOR_EACH_SEQ(what, x, ...) LIBRA_FOR_EACH(what, LIBRA_DEFER(LIBRA_COMMA), x, __VA_ARGS__)

// for each that creates seperation... just inserts semicolons between items
#define LIBRA_FOR_EACH_SEP(what, x, ...) LIBRA_FOR_EACH(what, LIBRA_DEFER(LIBRA_SEMICOLON), x, __VA_ARGS__)




// helper dummy type
template <typename... Args>
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
#define LIBRA_FUNCTOR_FOR_DEF(FunctionName) 									\
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






					

#define LIBRA_FUNCTOR_PREPARE(FunctionName, ...)	\
	LIBRA_HAS_METHOD_DEF(FunctionName);			\
	LIBRA_FUNCTOR_FOR_DEF(FunctionName);			


} // end namespace libra
#endif