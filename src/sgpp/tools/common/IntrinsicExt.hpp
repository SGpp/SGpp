/******************************************************************************
* Copyright (C) 2011 Technische Universitaet Muenchen                         *
* This file is part of the SG++ project. For conditions of distribution and   *
* use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
******************************************************************************/
// @author Alexander Heinecke (Alexander.Heinecke@mytum.de)

#ifndef INTRINSICEXT_HPP
#define INTRINSICEXT_HPP

#include <x86intrin.h>

#ifdef __SSE3__
const __m128 _mm_abs_ps( const __m128& x);
const __m128d _mm_abs_pd( const __m128d& x);
#endif

#ifdef __AVX__
const __m256 _mm256_abs_ps( const __m256& x);
const __m256d _mm256_abs_pd( const __m256d& x);
#endif

#endif /* INTRINSICEXT_HPP */
