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
#define ASSERT_INDEX_ARG(arg, min, max, alignment){\
	if(arg<min){\
		throw sg::base::operation_exception("argument to small!");\
	}\
	if(arg>max){\
		throw sg::base::operation_exception("argument to big!");\
	}\
	if(arg%alignment != 0){\
		std::cout << "argument " << arg << " not aligned to " << alignment << " (file:" << __FILE__ << ", line:" << __LINE__ << ")" << std::endl; \
		throw sg::base::operation_exception("argument not aligned!");\
	}\
}
#else
#define ASSERT_INDEX_ARG(arg, min, max, alignment)
#endif



#endif // KERNELMACROS_HPP
