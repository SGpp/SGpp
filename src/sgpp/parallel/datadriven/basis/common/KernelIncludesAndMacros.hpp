#ifndef KERNELMACROS_HPP
#define KERNELMACROS_HPP

#include "base/grid/GridStorage.hpp"

#if defined(__SSE3__) || defined(__AVX__)
#include <immintrin.h>
#endif
#if defined(__FMA4__)
#include <x86intrin.h>
#endif

#ifdef __USEAVX128__
#undef __AVX__
#endif

#define CHECK_KERNEL_CALLS
#ifdef CHECK_KERNEL_CALLS
#include "base/exception/operation_exception.hpp"

// helper macros
#define ASSERT_EQUAL(arg1, arg2) { \
	if ( (arg1) != (arg2) ) {\
		std::cerr << #arg1 << " and " << #arg2 << " are not equal: " << arg1 << " != " << arg2 \
		<< " (file:" << __FILE__ << ", line:" << __LINE__ << ")" << std::endl; \
		throw sg::base::operation_exception("values " #arg1 " and " #arg2 " are not equal");\
	}\
}
#define ASSERT_ALIGNMENT(arg, alignment) {\
	if( ( (arg)%(alignment) ) != 0 ) {\
		std::cout << #arg << "(" << arg << ") not aligned to " << #alignment << "(" << alignment << ") (file:" << __FILE__ << ", line:" << __LINE__ << ")" << std::endl; \
		throw sg::base::operation_exception("argument " #arg " not aligned!");\
	}\
}
#define ASSERT_GEQ_THAN(arg1, arg2) {\
	if( (arg1) < (arg2) ){\
		std::cout << #arg1 << "(" << arg1 << ") is smaller than " << #arg2 << "(" << arg2 << ") (file:" << __FILE__ << ", line:" << __LINE__ << ")" << std::endl; \
		throw sg::base::operation_exception(#arg1 " is smaller than " #arg2);\
	}\
}
#define ASSERT_LEQ_THAN(arg1, arg2) {\
	if( (arg1) > (arg2) ){\
		std::cout << #arg1 << "(" << arg1 << ") is greater than " << #arg2 << "(" << arg2 << ") (file:" << __FILE__ << ", line:" << __LINE__ << ")" << std::endl; \
		throw sg::base::operation_exception(#arg1 " is greater than " #arg2);\
	}\
}

#define CHECK_INDEX_ARG(arg, min, max, alignment){\
	ASSERT_GEQ_THAN(arg, min)\
	ASSERT_LEQ_THAN(arg, max)\
	ASSERT_ALIGNMENT(arg, alignment)\
}
#define CHECK_DATASET_AND_SOURCE(dataset, source) {\
	ASSERT_ALIGNMENT(dataset->getNcols(), getChunkDataPoints())\
	ASSERT_EQUAL(dataset->getNcols(), source.getSize()) \
}
#define CHECK_DATASET_AND_RESULT(dataset, result) {\
	ASSERT_ALIGNMENT(dataset->getNcols(), getChunkDataPoints())\
	ASSERT_EQUAL(dataset->getNcols(), result.getSize()) \
}

// use only the following two
#define CHECK_ARGS_MULT(level, dataset, result, s_grid, e_grid, s_data, e_data) {\
	CHECK_INDEX_ARG(s_grid, 0, level->getNrows(), 1);\
	CHECK_INDEX_ARG(e_grid, 0, level->getNrows(), 1);\
	CHECK_INDEX_ARG(s_data, 0, result.getSize(), getChunkDataPoints());\
	CHECK_INDEX_ARG(e_data, 0, result.getSize(), getChunkDataPoints());\
	CHECK_DATASET_AND_RESULT(dataset, result);\
}
#define CHECK_ARGS_MULTTRANSPOSE(level, dataset, source, s_grid, e_grid, s_data, e_data) {\
	CHECK_INDEX_ARG(s_grid, 0, level->getNrows(), 1);\
	CHECK_INDEX_ARG(e_grid, 0, level->getNrows(), 1);\
	CHECK_INDEX_ARG(s_data, 0, source.getSize(), getChunkDataPoints());\
	CHECK_INDEX_ARG(e_data, 0, source.getSize(), getChunkDataPoints());\
	CHECK_DATASET_AND_SOURCE(dataset, source);\
}
#else
#define CHECK_ARGS_MULT(level, dataset, result, s_grid, e_grid, s_data, e_data)
#define CHECK_ARGS_MULTTRANSPOSE(level, dataset, source, s_grid, e_grid, s_data, e_data)
#endif



#endif // KERNELMACROS_HPP
