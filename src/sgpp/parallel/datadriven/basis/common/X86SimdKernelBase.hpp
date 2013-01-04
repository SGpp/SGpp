/******************************************************************************
* Copyright (C) 2010 Technische Universitaet Muenchen                         *
* This file is part of the SG++ project. For conditions of distribution and   *
* use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
******************************************************************************/
// @author Alexander Heinecke (Alexander.Heinecke@mytum.de)
// @author Roman Karlstetter (karlstetter@mytum.de)

#ifndef X86SIMDKERNELBASE_HPP
#define X86SIMDKERNELBASE_HPP

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

namespace sg {
namespace parallel {

class X86SimdKernelBase
{
public:
	static inline size_t getChunkGridPoints(){return 12;}
	static inline size_t getChunkDataPoints(){return 24;} // must be divisible by 24
};

}
}

#endif // X86SIMDKERNELBASE_HPP
