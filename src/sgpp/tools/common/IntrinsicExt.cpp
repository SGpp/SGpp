/******************************************************************************
* Copyright (C) 2011 Technische Universitaet Muenchen                         *
* This file is part of the SG++ project. For conditions of distribution and   *
* use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
******************************************************************************/
// @author Alexander Heinecke (Alexander.Heinecke@mytum.de)

#include "tools/common/IntrinsicExt.hpp"

#ifdef __SSE3__
union floatAbsMask
{
   const float f;
   const int i;

   floatAbsMask() : i(0x7FFFFFFF) {}
};

__attribute__((aligned(16))) const floatAbsMask absMask;
static const __m128 abs2Mask = _mm_load1_ps( &absMask.f );

const __m128 _mm_abs_ps( const __m128& x)
{
       return _mm_and_ps( abs2Mask, x);
}

union doubleAbsMask
{
   const double d;
   const long long ilong;

   doubleAbsMask() : ilong(0x7FFFFFFFFFFFFFFF) {}
};

__attribute__((aligned(16))) const doubleAbsMask absMaskd;
static const __m128d abs2Maskd = _mm_load1_pd( &absMaskd.d );

const __m128d _mm_abs_pd( const __m128d& x)
{
       return _mm_and_pd( abs2Maskd, x);
}
#endif

#ifdef __AVX__
union floatAbsMaskAVX
{
   const float fAVX;
   const int iAVX;

   floatAbsMaskAVX() : iAVX(0x7FFFFFFF) {}
};

__attribute__((aligned(32))) const floatAbsMaskAVX absMaskSPAVX;

static const __m256 abs2MaskSPAVX = _mm256_broadcast_ss( &(absMaskSPAVX.fAVX) );

const __m256 _mm256_abs_ps( const __m256& x)
{
       return _mm256_and_ps( abs2MaskSPAVX, x);
}

union doubleAbsMaskAVX
{
   const double dAVX;
   const long long ilongAVX;

   doubleAbsMaskAVX() : ilongAVX(0x7FFFFFFFFFFFFFFF) {}
};

__attribute__((aligned(32))) const doubleAbsMaskAVX absMaskAVX;

static const __m256d abs2MaskAVX = _mm256_broadcast_sd( &(absMaskAVX.dAVX) );

const __m256d _mm256_abs_pd( const __m256d& x)
{
       return _mm256_and_pd( abs2MaskAVX, x);
}
#endif
