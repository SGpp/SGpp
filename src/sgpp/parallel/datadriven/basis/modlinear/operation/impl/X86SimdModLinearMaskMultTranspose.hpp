/******************************************************************************
* Copyright (C) 2010 Technische Universitaet Muenchen                         *
* This file is part of the SG++ project. For conditions of distribution and   *
* use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
******************************************************************************/
// @author Alexander Heinecke (Alexander.Heinecke@mytum.de)
// @author Roman Karlstetter (karlstetter@mytum.de)

#ifndef X86SIMDMODLINEARMASKMULTTRANSPOSE_HPP
#define X86SIMDMODLINEARMASKMULTTRANSPOSE_HPP

#if defined(__SSE3__) || defined(__AVX__)
#include <immintrin.h>
#endif
#if defined(__FMA4__)
#include <x86intrin.h>
#endif

#ifdef __USEAVX128__
#undef __AVX__
#endif

#include "base/grid/GridStorage.hpp"
namespace sg {
namespace parallel {

class X86SimdModLinearMultTranspose
{
public:
	static inline size_t getChunkGridPoints(){return 12;}
	static inline size_t getChunkDataPoints(){return 24;}

	/// @todo implement
};

}
}

#endif // X86SIMDMODLINEARMASKMULTTRANSPOSE_HPP
