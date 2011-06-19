/******************************************************************************
* Copyright (C) 2011 Technische Universitaet Muenchen                         *
* This file is part of the SG++ project. For conditions of distribution and   *
* use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
******************************************************************************/
// @author Alexander Heinecke (Alexander.Heinecke@mytum.de)

#ifdef __ICC
// include SSE3 intrinsics
#include <pmmintrin.h>

union floatAbsMask
{
   const float f;
   const int i;

   floatAbsMask() : i(0x7FFFFFFF) {}
};

_MM_ALIGN16 const floatAbsMask absMask;
static const __m128 abs2Mask = _mm_load1_ps( &absMask.f );

const __m128 _mm_abs_ps( const __m128& x)
{
       return _mm_and_ps( abs2Mask, x);
}
#endif
