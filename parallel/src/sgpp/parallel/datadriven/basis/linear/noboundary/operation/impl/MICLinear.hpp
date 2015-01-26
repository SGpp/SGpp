/* ****************************************************************************
* Copyright (C) 2012 - 2013 Technische Universitaet Muenchen                  *
* This file is part of the SG++ project. For conditions of distribution and   *
* use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
**************************************************************************** */
// @author Alexander Heinecke (Alexander.Heinecke@mytum.de)
// @author Roman Karlstetter (karlstetter@mytum.de)

#ifndef MICLINEAR_HPP
#define MICLINEAR_HPP
#ifdef USEMIC

#ifdef __MIC__
#include <immintrin.h>
#endif

#include <sgpp/parallel/datadriven/basis/common/mic/MICKernelBase.hpp>

#include <sgpp/globaldef.hpp>


namespace SGPP {
  namespace parallel {
    /**
     * @brief The MICLinear class
     *
     * multImpl and multTransImpl have to be run on the mic, either in offload mode or native
     */
    class MICLinear : public MICKernelBase {
      public:
        static const KernelType kernelType = Standard;
        static inline void multImpl(
          double* ptrLevel,
          double* ptrIndex,
          double* /*ptrMask*/, //unused for this specialization
          double* /*ptrOffset*/, //unused for this specialization
          double* ptrData,
          double* ptrAlpha,
          double* ptrResult,
          size_t result_size,
          size_t dims,
          const size_t start_index_grid,
          const size_t end_index_grid,
          const size_t start_index_data,
          const size_t end_index_data) {
          for (size_t i = start_index_data; i < end_index_data; i += getChunkDataPoints()) {
#ifdef __MIC__

            for (size_t j = start_index_grid; j < end_index_grid; j++) {
              _mm_prefetch((const char*) & (ptrAlpha[j + 1]), _MM_HINT_T0);
              _mm_prefetch((const char*) & (ptrLevel[((j + 1)*dims)]), _MM_HINT_T0);
              _mm_prefetch((const char*) & (ptrIndex[((j + 1)*dims)]), _MM_HINT_T0);

              __m512d support_0 = _mm512_extload_pd(&(ptrAlpha[j]), _MM_UPCONV_PD_NONE, _MM_BROADCAST_1X8, _MM_HINT_NONE);
              __m512d support_1 = _mm512_extload_pd(&(ptrAlpha[j]), _MM_UPCONV_PD_NONE, _MM_BROADCAST_1X8, _MM_HINT_NONE);
              __m512d support_2 = _mm512_extload_pd(&(ptrAlpha[j]), _MM_UPCONV_PD_NONE, _MM_BROADCAST_1X8, _MM_HINT_NONE);
#if  ( MIC_UNROLLING_WIDTH >  24 )
              __m512d support_3 = _mm512_extload_pd(&(ptrAlpha[j]), _MM_UPCONV_PD_NONE, _MM_BROADCAST_1X8, _MM_HINT_NONE);
#endif
#if  ( MIC_UNROLLING_WIDTH >  32 )
              __m512d support_4 = _mm512_extload_pd(&(ptrAlpha[j]), _MM_UPCONV_PD_NONE, _MM_BROADCAST_1X8, _MM_HINT_NONE);
              __m512d support_5 = _mm512_extload_pd(&(ptrAlpha[j]), _MM_UPCONV_PD_NONE, _MM_BROADCAST_1X8, _MM_HINT_NONE);
#endif
#if  ( MIC_UNROLLING_WIDTH >  48 )
              __m512d support_6 = _mm512_extload_pd(&(ptrAlpha[j]), _MM_UPCONV_PD_NONE, _MM_BROADCAST_1X8, _MM_HINT_NONE);
              __m512d support_7 = _mm512_extload_pd(&(ptrAlpha[j]), _MM_UPCONV_PD_NONE, _MM_BROADCAST_1X8, _MM_HINT_NONE);
#endif
#if  ( MIC_UNROLLING_WIDTH >  64 )
              __m512d support_8 = _mm512_extload_pd(&(ptrAlpha[j]), _MM_UPCONV_PD_NONE, _MM_BROADCAST_1X8, _MM_HINT_NONE);
#endif
#if  ( MIC_UNROLLING_WIDTH >  72 )
              __m512d support_9 = _mm512_extload_pd(&(ptrAlpha[j]), _MM_UPCONV_PD_NONE, _MM_BROADCAST_1X8, _MM_HINT_NONE);
#endif
#if  ( MIC_UNROLLING_WIDTH >  80 )
              __m512d support_10 = _mm512_extload_pd(&(ptrAlpha[j]), _MM_UPCONV_PD_NONE, _MM_BROADCAST_1X8, _MM_HINT_NONE);
              __m512d support_11 = _mm512_extload_pd(&(ptrAlpha[j]), _MM_UPCONV_PD_NONE, _MM_BROADCAST_1X8, _MM_HINT_NONE);
#endif

              for (size_t d = 0; d < dims - 1; d++) {
                __m512d eval_0 = _mm512_load_pd(&(ptrData[(d * result_size) + i + 0]));
                __m512d eval_1 = _mm512_load_pd(&(ptrData[(d * result_size) + i + 8]));
                __m512d eval_2 = _mm512_load_pd(&(ptrData[(d * result_size) + i + 16]));
#if  ( MIC_UNROLLING_WIDTH >  24 )
                __m512d eval_3 = _mm512_load_pd(&(ptrData[(d * result_size) + i + 24]));
#endif
#if  ( MIC_UNROLLING_WIDTH >  32 )
                __m512d eval_4 = _mm512_load_pd(&(ptrData[(d * result_size) + i + 32]));
                __m512d eval_5 = _mm512_load_pd(&(ptrData[(d * result_size) + i + 40]));
#endif
#if  ( MIC_UNROLLING_WIDTH >  48 )
                __m512d eval_6 = _mm512_load_pd(&(ptrData[(d * result_size) + i + 48]));
                __m512d eval_7 = _mm512_load_pd(&(ptrData[(d * result_size) + i + 56]));
#endif
#if  ( MIC_UNROLLING_WIDTH >  64 )
                __m512d eval_8 = _mm512_load_pd(&(ptrData[(d * result_size) + i + 64]));
#endif
#if  ( MIC_UNROLLING_WIDTH >  72 )
                __m512d eval_9 = _mm512_load_pd(&(ptrData[(d * result_size) + i + 72]));
#endif
#if  ( MIC_UNROLLING_WIDTH >  80 )
                __m512d eval_10 = _mm512_load_pd(&(ptrData[(d * result_size) + i + 80]));
                __m512d eval_11 = _mm512_load_pd(&(ptrData[(d * result_size) + i + 88]));
#endif
                __m512d level = _mm512_extload_pd(&(ptrLevel[(j * dims) + d]), _MM_UPCONV_PD_NONE, _MM_BROADCAST_1X8, _MM_HINT_NONE);
                __m512d index = _mm512_extload_pd(&(ptrIndex[(j * dims) + d]), _MM_UPCONV_PD_NONE, _MM_BROADCAST_1X8, _MM_HINT_NONE);

                eval_0 = _mm512_fmsub_pd(eval_0, level, index);
                eval_1 = _mm512_fmsub_pd(eval_1, level, index);
                eval_2 = _mm512_fmsub_pd(eval_2, level, index);
#if  ( MIC_UNROLLING_WIDTH >  24 )
                eval_3 = _mm512_fmsub_pd(eval_3, level, index);
#endif
#if  ( MIC_UNROLLING_WIDTH >  32 )
                eval_4 = _mm512_fmsub_pd(eval_4, level, index);
                eval_5 = _mm512_fmsub_pd(eval_5, level, index);
#endif
#if  ( MIC_UNROLLING_WIDTH >  48 )
                eval_6 = _mm512_fmsub_pd(eval_6, level, index);
                eval_7 = _mm512_fmsub_pd(eval_7, level, index);
#endif
#if  ( MIC_UNROLLING_WIDTH >  64 )
                eval_8 = _mm512_fmsub_pd(eval_8, level, index);
#endif
#if  ( MIC_UNROLLING_WIDTH >  72 )
                eval_9 = _mm512_fmsub_pd(eval_9, level, index);
#endif
#if  ( MIC_UNROLLING_WIDTH >  80 )
                eval_10 = _mm512_fmsub_pd(eval_10, level, index);
                eval_11 = _mm512_fmsub_pd(eval_11, level, index);
#endif
                __m512i abs2MaskLRBni = _mm512_set_1to8_epi64(0x7FFFFFFFFFFFFFFF);
                eval_0 = _mm512_castsi512_pd(_mm512_and_epi64( abs2MaskLRBni, _mm512_castpd_si512(eval_0)));
                eval_1 = _mm512_castsi512_pd(_mm512_and_epi64( abs2MaskLRBni, _mm512_castpd_si512(eval_1)));
                eval_2 = _mm512_castsi512_pd(_mm512_and_epi64( abs2MaskLRBni, _mm512_castpd_si512(eval_2)));
#if  ( MIC_UNROLLING_WIDTH >  24 )
                eval_3 = _mm512_castsi512_pd(_mm512_and_epi64( abs2MaskLRBni, _mm512_castpd_si512(eval_3)));
#endif
#if  ( MIC_UNROLLING_WIDTH >  32 )
                eval_4 = _mm512_castsi512_pd(_mm512_and_epi64( abs2MaskLRBni, _mm512_castpd_si512(eval_4)));
                eval_5 = _mm512_castsi512_pd(_mm512_and_epi64( abs2MaskLRBni, _mm512_castpd_si512(eval_5)));
#endif
#if  ( MIC_UNROLLING_WIDTH >  48 )
                eval_6 = _mm512_castsi512_pd(_mm512_and_epi64( abs2MaskLRBni, _mm512_castpd_si512(eval_6)));
                eval_7 = _mm512_castsi512_pd(_mm512_and_epi64( abs2MaskLRBni, _mm512_castpd_si512(eval_7)));
#endif
#if  ( MIC_UNROLLING_WIDTH >  64 )
                eval_8 = _mm512_castsi512_pd(_mm512_and_epi64( abs2MaskLRBni, _mm512_castpd_si512(eval_8)));
#endif
#if  ( MIC_UNROLLING_WIDTH >  72 )
                eval_9 = _mm512_castsi512_pd(_mm512_and_epi64( abs2MaskLRBni, _mm512_castpd_si512(eval_9)));
#endif
#if  ( MIC_UNROLLING_WIDTH >  80 )
                eval_10 = _mm512_castsi512_pd(_mm512_and_epi64( abs2MaskLRBni, _mm512_castpd_si512(eval_10)));
                eval_11 = _mm512_castsi512_pd(_mm512_and_epi64( abs2MaskLRBni, _mm512_castpd_si512(eval_11)));
#endif

                __m512d one = _mm512_set_1to8_pd(1.0);


                eval_0 = _mm512_sub_pd(one, eval_0);
                eval_1 = _mm512_sub_pd(one, eval_1);
                eval_2 = _mm512_sub_pd(one, eval_2);
#if  ( MIC_UNROLLING_WIDTH >  24 )
                eval_3 = _mm512_sub_pd(one, eval_3);
#endif
#if  ( MIC_UNROLLING_WIDTH >  32 )
                eval_4 = _mm512_sub_pd(one, eval_4);
                eval_5 = _mm512_sub_pd(one, eval_5);
#endif
#if  ( MIC_UNROLLING_WIDTH >  48 )
                eval_6 = _mm512_sub_pd(one, eval_6);
                eval_7 = _mm512_sub_pd(one, eval_7);
#endif
#if  ( MIC_UNROLLING_WIDTH >  64 )
                eval_8 = _mm512_sub_pd(one, eval_8);
#endif
#if  ( MIC_UNROLLING_WIDTH >  72 )
                eval_9 = _mm512_sub_pd(one, eval_9);
#endif
#if  ( MIC_UNROLLING_WIDTH >  80 )
                eval_10 = _mm512_sub_pd(one, eval_10);
                eval_11 = _mm512_sub_pd(one, eval_11);
#endif
                __m512d zero = _mm512_set_1to8_pd(0.0);

                eval_0 = _mm512_gmax_pd(zero, eval_0);
                eval_1 = _mm512_gmax_pd(zero, eval_1);
                eval_2 = _mm512_gmax_pd(zero, eval_2);
#if  ( MIC_UNROLLING_WIDTH >  24 )
                eval_3 = _mm512_gmax_pd(zero, eval_3);
#endif
#if  ( MIC_UNROLLING_WIDTH >  32 )
                eval_4 = _mm512_gmax_pd(zero, eval_4);
                eval_5 = _mm512_gmax_pd(zero, eval_5);
#endif
#if  ( MIC_UNROLLING_WIDTH >  48 )
                eval_6 = _mm512_gmax_pd(zero, eval_6);
                eval_7 = _mm512_gmax_pd(zero, eval_7);
#endif
#if  ( MIC_UNROLLING_WIDTH >  64 )
                eval_8 = _mm512_gmax_pd(zero, eval_8);
#endif
#if  ( MIC_UNROLLING_WIDTH >  72 )
                eval_9 = _mm512_gmax_pd(zero, eval_9);
#endif
#if  ( MIC_UNROLLING_WIDTH >  80 )
                eval_10 = _mm512_gmax_pd(zero, eval_10);
                eval_11 = _mm512_gmax_pd(zero, eval_11);
#endif

                support_0 = _mm512_mul_pd(support_0, eval_0);
                support_1 = _mm512_mul_pd(support_1, eval_1);
                support_2 = _mm512_mul_pd(support_2, eval_2);
#if  ( MIC_UNROLLING_WIDTH >  24 )
                support_3 = _mm512_mul_pd(support_3, eval_3);
#endif
#if  ( MIC_UNROLLING_WIDTH >  32 )
                support_4 = _mm512_mul_pd(support_4, eval_4);
                support_5 = _mm512_mul_pd(support_5, eval_5);
#endif
#if  ( MIC_UNROLLING_WIDTH >  48 )
                support_6 = _mm512_mul_pd(support_6, eval_6);
                support_7 = _mm512_mul_pd(support_7, eval_7);
#endif
#if  ( MIC_UNROLLING_WIDTH >  64 )
                support_8 = _mm512_mul_pd(support_8, eval_8);
#endif
#if  ( MIC_UNROLLING_WIDTH >  72 )
                support_9 = _mm512_mul_pd(support_9, eval_9);
#endif
#if  ( MIC_UNROLLING_WIDTH >  80 )
                support_10 = _mm512_mul_pd(support_10, eval_10);
                support_11 = _mm512_mul_pd(support_11, eval_11);
#endif
              }

              size_t d = dims - 1;

              __m512d eval_0 = _mm512_load_pd(&(ptrData[(d * result_size) + i + 0]));
              __m512d eval_1 = _mm512_load_pd(&(ptrData[(d * result_size) + i + 8]));
              __m512d eval_2 = _mm512_load_pd(&(ptrData[(d * result_size) + i + 16]));
#if  ( MIC_UNROLLING_WIDTH >  24 )
              __m512d eval_3 = _mm512_load_pd(&(ptrData[(d * result_size) + i + 24]));
#endif
#if  ( MIC_UNROLLING_WIDTH >  32 )
              __m512d eval_4 = _mm512_load_pd(&(ptrData[(d * result_size) + i + 32]));
              __m512d eval_5 = _mm512_load_pd(&(ptrData[(d * result_size) + i + 40]));
#endif
#if  ( MIC_UNROLLING_WIDTH >  48 )
              __m512d eval_6 = _mm512_load_pd(&(ptrData[(d * result_size) + i + 48]));
              __m512d eval_7 = _mm512_load_pd(&(ptrData[(d * result_size) + i + 56]));
#endif
#if  ( MIC_UNROLLING_WIDTH >  64 )
              __m512d eval_8 = _mm512_load_pd(&(ptrData[(d * result_size) + i + 64]));
#endif
#if  ( MIC_UNROLLING_WIDTH >  72 )
              __m512d eval_9 = _mm512_load_pd(&(ptrData[(d * result_size) + i + 72]));
#endif
#if  ( MIC_UNROLLING_WIDTH >  80 )
              __m512d eval_10 = _mm512_load_pd(&(ptrData[(d * result_size) + i + 80]));
              __m512d eval_11 = _mm512_load_pd(&(ptrData[(d * result_size) + i + 88]));
#endif
              __m512d level = _mm512_extload_pd(&(ptrLevel[(j * dims) + d]), _MM_UPCONV_PD_NONE, _MM_BROADCAST_1X8, _MM_HINT_NONE);
              __m512d index = _mm512_extload_pd(&(ptrIndex[(j * dims) + d]), _MM_UPCONV_PD_NONE, _MM_BROADCAST_1X8, _MM_HINT_NONE);

              eval_0 = _mm512_fmsub_pd(eval_0, level, index);
              eval_1 = _mm512_fmsub_pd(eval_1, level, index);
              eval_2 = _mm512_fmsub_pd(eval_2, level, index);
#if  ( MIC_UNROLLING_WIDTH >  24 )
              eval_3 = _mm512_fmsub_pd(eval_3, level, index);
#endif
#if  ( MIC_UNROLLING_WIDTH >  32 )
              eval_4 = _mm512_fmsub_pd(eval_4, level, index);
              eval_5 = _mm512_fmsub_pd(eval_5, level, index);
#endif
#if  ( MIC_UNROLLING_WIDTH >  48 )
              eval_6 = _mm512_fmsub_pd(eval_6, level, index);
              eval_7 = _mm512_fmsub_pd(eval_7, level, index);
#endif
#if  ( MIC_UNROLLING_WIDTH >  64 )
              eval_8 = _mm512_fmsub_pd(eval_8, level, index);
#endif
#if  ( MIC_UNROLLING_WIDTH >  72 )
              eval_9 = _mm512_fmsub_pd(eval_9, level, index);
#endif
#if  ( MIC_UNROLLING_WIDTH >  80 )
              eval_10 = _mm512_fmsub_pd(eval_10, level, index);
              eval_11 = _mm512_fmsub_pd(eval_11, level, index);
#endif
              __m512i abs2MaskLRBni = _mm512_set_1to8_epi64(0x7FFFFFFFFFFFFFFF);
              eval_0 = _mm512_castsi512_pd(_mm512_and_epi64( abs2MaskLRBni, _mm512_castpd_si512(eval_0)));
              eval_1 = _mm512_castsi512_pd(_mm512_and_epi64( abs2MaskLRBni, _mm512_castpd_si512(eval_1)));
              eval_2 = _mm512_castsi512_pd(_mm512_and_epi64( abs2MaskLRBni, _mm512_castpd_si512(eval_2)));
#if  ( MIC_UNROLLING_WIDTH >  24 )
              eval_3 = _mm512_castsi512_pd(_mm512_and_epi64( abs2MaskLRBni, _mm512_castpd_si512(eval_3)));
#endif
#if  ( MIC_UNROLLING_WIDTH >  32 )
              eval_4 = _mm512_castsi512_pd(_mm512_and_epi64( abs2MaskLRBni, _mm512_castpd_si512(eval_4)));
              eval_5 = _mm512_castsi512_pd(_mm512_and_epi64( abs2MaskLRBni, _mm512_castpd_si512(eval_5)));
#endif
#if  ( MIC_UNROLLING_WIDTH >  48 )
              eval_6 = _mm512_castsi512_pd(_mm512_and_epi64( abs2MaskLRBni, _mm512_castpd_si512(eval_6)));
              eval_7 = _mm512_castsi512_pd(_mm512_and_epi64( abs2MaskLRBni, _mm512_castpd_si512(eval_7)));
#endif
#if  ( MIC_UNROLLING_WIDTH >  64 )
              eval_8 = _mm512_castsi512_pd(_mm512_and_epi64( abs2MaskLRBni, _mm512_castpd_si512(eval_8)));
#endif
#if  ( MIC_UNROLLING_WIDTH >  72 )
              eval_9 = _mm512_castsi512_pd(_mm512_and_epi64( abs2MaskLRBni, _mm512_castpd_si512(eval_9)));
#endif
#if  ( MIC_UNROLLING_WIDTH >  80 )
              eval_10 = _mm512_castsi512_pd(_mm512_and_epi64( abs2MaskLRBni, _mm512_castpd_si512(eval_10)));
              eval_11 = _mm512_castsi512_pd(_mm512_and_epi64( abs2MaskLRBni, _mm512_castpd_si512(eval_11)));
#endif
              __m512d one = _mm512_set_1to8_pd(1.0);

              eval_0 = _mm512_sub_pd(one, eval_0);
              eval_1 = _mm512_sub_pd(one, eval_1);
              eval_2 = _mm512_sub_pd(one, eval_2);
#if  ( MIC_UNROLLING_WIDTH >  24 )
              eval_3 = _mm512_sub_pd(one, eval_3);
#endif
#if  ( MIC_UNROLLING_WIDTH >  32 )
              eval_4 = _mm512_sub_pd(one, eval_4);
              eval_5 = _mm512_sub_pd(one, eval_5);
#endif
#if  ( MIC_UNROLLING_WIDTH >  48 )
              eval_6 = _mm512_sub_pd(one, eval_6);
              eval_7 = _mm512_sub_pd(one, eval_7);
#endif
#if  ( MIC_UNROLLING_WIDTH >  64 )
              eval_8 = _mm512_sub_pd(one, eval_8);
#endif
#if  ( MIC_UNROLLING_WIDTH >  72 )
              eval_9 = _mm512_sub_pd(one, eval_9);
#endif
#if  ( MIC_UNROLLING_WIDTH >  80 )
              eval_10 = _mm512_sub_pd(one, eval_10);
              eval_11 = _mm512_sub_pd(one, eval_11);
#endif
              __m512d zero = _mm512_set_1to8_pd(0.0);

              eval_0 = _mm512_gmax_pd(zero, eval_0);
              eval_1 = _mm512_gmax_pd(zero, eval_1);
              eval_2 = _mm512_gmax_pd(zero, eval_2);
#if  ( MIC_UNROLLING_WIDTH >  24 )
              eval_3 = _mm512_gmax_pd(zero, eval_3);
#endif
#if  ( MIC_UNROLLING_WIDTH >  32 )
              eval_4 = _mm512_gmax_pd(zero, eval_4);
              eval_5 = _mm512_gmax_pd(zero, eval_5);
#endif
#if  ( MIC_UNROLLING_WIDTH >  48 )
              eval_6 = _mm512_gmax_pd(zero, eval_6);
              eval_7 = _mm512_gmax_pd(zero, eval_7);
#endif
#if  ( MIC_UNROLLING_WIDTH >  64 )
              eval_8 = _mm512_gmax_pd(zero, eval_8);
#endif
#if  ( MIC_UNROLLING_WIDTH >  72 )
              eval_9 = _mm512_gmax_pd(zero, eval_9);
#endif
#if  ( MIC_UNROLLING_WIDTH >  80 )
              eval_10 = _mm512_gmax_pd(zero, eval_10);
              eval_11 = _mm512_gmax_pd(zero, eval_11);
#endif
              __m512d res_0 = _mm512_load_pd(&(ptrResult[i + 0]));
              __m512d res_1 = _mm512_load_pd(&(ptrResult[i + 8]));
              __m512d res_2 = _mm512_load_pd(&(ptrResult[i + 16]));
#if  ( MIC_UNROLLING_WIDTH >  24 )
              __m512d res_3 = _mm512_load_pd(&(ptrResult[i + 24]));
#endif
#if  ( MIC_UNROLLING_WIDTH >  32 )
              __m512d res_4 = _mm512_load_pd(&(ptrResult[i + 32]));
              __m512d res_5 = _mm512_load_pd(&(ptrResult[i + 40]));
#endif

              res_0 = _mm512_add_pd(res_0, _mm512_mul_pd(support_0, eval_0));
              res_1 = _mm512_add_pd(res_1, _mm512_mul_pd(support_1, eval_1));
              res_2 = _mm512_add_pd(res_2, _mm512_mul_pd(support_2, eval_2));
#if  ( MIC_UNROLLING_WIDTH >  24 )
              res_3 = _mm512_add_pd(res_3, _mm512_mul_pd(support_3, eval_3));
#endif
#if  ( MIC_UNROLLING_WIDTH >  32 )
              res_4 = _mm512_add_pd(res_4, _mm512_mul_pd(support_4, eval_4));
              res_5 = _mm512_add_pd(res_5, _mm512_mul_pd(support_5, eval_5));
#endif

              _mm512_store_pd(&(ptrResult[i + 0]), res_0);
              _mm512_store_pd(&(ptrResult[i + 8]), res_1);
              _mm512_store_pd(&(ptrResult[i + 16]), res_2);
#if  ( MIC_UNROLLING_WIDTH >  24 )
              _mm512_store_pd(&(ptrResult[i + 24]), res_3);
#endif
#if  ( MIC_UNROLLING_WIDTH >  32 )
              _mm512_store_pd(&(ptrResult[i + 32]), res_4);
              _mm512_store_pd(&(ptrResult[i + 40]), res_5);
#endif

#if  ( MIC_UNROLLING_WIDTH >  48 )
              __m512d res_6 = _mm512_load_pd(&(ptrResult[i + 48]));
              __m512d res_7 = _mm512_load_pd(&(ptrResult[i + 56]));
#endif
#if  ( MIC_UNROLLING_WIDTH >  64 )
              __m512d res_8 = _mm512_load_pd(&(ptrResult[i + 64]));
#endif
#if  ( MIC_UNROLLING_WIDTH >  72 )
              __m512d res_9 = _mm512_load_pd(&(ptrResult[i + 72]));
#endif
#if  ( MIC_UNROLLING_WIDTH >  80 )
              __m512d res_10 = _mm512_load_pd(&(ptrResult[i + 80]));
              __m512d res_11 = _mm512_load_pd(&(ptrResult[i + 88]));
#endif

#if  ( MIC_UNROLLING_WIDTH >  48 )
              res_6 = _mm512_add_pd(res_6, _mm512_mul_pd(support_6, eval_6));
              res_7 = _mm512_add_pd(res_7, _mm512_mul_pd(support_7, eval_7));
#endif
#if  ( MIC_UNROLLING_WIDTH >  64 )
              res_8 = _mm512_add_pd(res_8, _mm512_mul_pd(support_8, eval_8));
#endif
#if  ( MIC_UNROLLING_WIDTH >  72 )
              res_9 = _mm512_add_pd(res_9, _mm512_mul_pd(support_9, eval_9));
#endif
#if  ( MIC_UNROLLING_WIDTH >  80 )
              res_10 = _mm512_add_pd(res_10, _mm512_mul_pd(support_10, eval_10));
              res_11 = _mm512_add_pd(res_11, _mm512_mul_pd(support_11, eval_11));
#endif

#if  ( MIC_UNROLLING_WIDTH >  48 )
              _mm512_store_pd(&(ptrResult[i + 48]), res_6);
              _mm512_store_pd(&(ptrResult[i + 56]), res_7);
#endif
#if  ( MIC_UNROLLING_WIDTH >  64 )
              _mm512_store_pd(&(ptrResult[i + 64]), res_8);
#endif
#if  ( MIC_UNROLLING_WIDTH >  72 )
              _mm512_store_pd(&(ptrResult[i + 72]), res_9);
#endif
#if  ( MIC_UNROLLING_WIDTH >  80 )
              _mm512_store_pd(&(ptrResult[i + 80]), res_10);
              _mm512_store_pd(&(ptrResult[i + 88]), res_11);
#endif
            }

#else

            for (size_t ii = i; ii < i + getChunkDataPoints(); ii++) {
              for (size_t j = start_index_grid; j < end_index_grid; j++) {
                double curSupport = ptrAlpha[j];

                for (size_t d = 0; d < dims; d++) {
                  double eval = ((ptrLevel[(j * dims) + d]) * (ptrData[(d * result_size) + ii]));
                  double index_calc = eval - (ptrIndex[(j * dims) + d]);
                  double abs = fabs(index_calc);
                  double last = 1.0 - abs;
                  double localSupport = std::max<double>(last, 0.0);
                  curSupport *= localSupport;
                }

                ptrResult[ii] += curSupport;
              }
            }

#endif
          }

        }

        static inline void multTransposeImpl(
          double* ptrLevel,
          double* ptrIndex,
          double* /*ptrMask*/, //unused for this specialization
          double* /*ptrOffset*/, //unused for this specialization
          double* ptrData,
          double* ptrSource,
          double* ptrResult,
          size_t sourceSize,
          size_t dims,
          const size_t start_index_grid,
          const size_t end_index_grid,
          const size_t start_index_data,
          const size_t end_index_data) {
          if (start_index_grid >= end_index_grid) {
            // special handling for grid index (prefetch for last item)
            // so we have to check it here
            return;
          }

#ifdef __MIC__

          for (size_t i = start_index_data; i < end_index_data; i += getChunkDataPoints()) {
            for (size_t j = start_index_grid; j < end_index_grid - 1; j++) {
              __m512d support_0 = _mm512_load_pd(&(ptrSource[i + 0]));
              __m512d support_1 = _mm512_load_pd(&(ptrSource[i + 8]));
              __m512d support_2 = _mm512_load_pd(&(ptrSource[i + 16]));
#if  ( MIC_UNROLLING_WIDTH >  24 )
              __m512d support_3 = _mm512_load_pd(&(ptrSource[i + 24]));
#endif
#if  ( MIC_UNROLLING_WIDTH >  32 )
              __m512d support_4 = _mm512_load_pd(&(ptrSource[i + 32]));
              __m512d support_5 = _mm512_load_pd(&(ptrSource[i + 40]));
#endif
#if  ( MIC_UNROLLING_WIDTH >  48 )
              __m512d support_6 = _mm512_load_pd(&(ptrSource[i + 48]));
              __m512d support_7 = _mm512_load_pd(&(ptrSource[i + 56]));
#endif
#if  ( MIC_UNROLLING_WIDTH >  64 )
              __m512d support_8 = _mm512_load_pd(&(ptrSource[i + 64]));
#endif
#if  ( MIC_UNROLLING_WIDTH >  72 )
              __m512d support_9 = _mm512_load_pd(&(ptrSource[i + 72]));
#endif
#if  ( MIC_UNROLLING_WIDTH >  80 )
              __m512d support_10 = _mm512_load_pd(&(ptrSource[i + 80]));
              __m512d support_11 = _mm512_load_pd(&(ptrSource[i + 88]));
#endif

              _mm_prefetch((const char*) & (ptrLevel[((j + 1)*dims)]), _MM_HINT_T0);
              _mm_prefetch((const char*) & (ptrIndex[((j + 1)*dims)]), _MM_HINT_T0);
              // Prefetch exclusive -> non sync is needed
              _mm_prefetch((const char*) & (ptrResult[j]), _MM_HINT_ET0);

              for (size_t d = 0; d < dims; d++) {
                __m512d eval_0 = _mm512_load_pd(&(ptrData[(d * sourceSize) + i + 0]));
                __m512d eval_1 = _mm512_load_pd(&(ptrData[(d * sourceSize) + i + 8]));
                __m512d eval_2 = _mm512_load_pd(&(ptrData[(d * sourceSize) + i + 16]));
#if  ( MIC_UNROLLING_WIDTH >  24 )
                __m512d eval_3 = _mm512_load_pd(&(ptrData[(d * sourceSize) + i + 24]));
#endif
#if  ( MIC_UNROLLING_WIDTH >  32 )
                __m512d eval_4 = _mm512_load_pd(&(ptrData[(d * sourceSize) + i + 32]));
                __m512d eval_5 = _mm512_load_pd(&(ptrData[(d * sourceSize) + i + 40]));
#endif
#if  ( MIC_UNROLLING_WIDTH >  48 )
                __m512d eval_6 = _mm512_load_pd(&(ptrData[(d * sourceSize) + i + 48]));
                __m512d eval_7 = _mm512_load_pd(&(ptrData[(d * sourceSize) + i + 56]));
#endif
#if  ( MIC_UNROLLING_WIDTH >  64 )
                __m512d eval_8 = _mm512_load_pd(&(ptrData[(d * sourceSize) + i + 64]));
#endif
#if  ( MIC_UNROLLING_WIDTH >  72 )
                __m512d eval_9 = _mm512_load_pd(&(ptrData[(d * sourceSize) + i + 72]));
#endif
#if  ( MIC_UNROLLING_WIDTH >  80 )
                __m512d eval_10 = _mm512_load_pd(&(ptrData[(d * sourceSize) + i + 80]));
                __m512d eval_11 = _mm512_load_pd(&(ptrData[(d * sourceSize) + i + 88]));
#endif

                __m512d level = _mm512_extload_pd(&(ptrLevel[(j * dims) + d]),  _MM_UPCONV_PD_NONE, _MM_BROADCAST_1X8, _MM_HINT_NONE);
                __m512d index = _mm512_extload_pd(&(ptrIndex[(j * dims) + d]),  _MM_UPCONV_PD_NONE, _MM_BROADCAST_1X8, _MM_HINT_NONE);

                eval_0 = _mm512_fmsub_pd(eval_0, level, index);
                eval_1 = _mm512_fmsub_pd(eval_1, level, index);
                eval_2 = _mm512_fmsub_pd(eval_2, level, index);
#if  ( MIC_UNROLLING_WIDTH >  24 )
                eval_3 = _mm512_fmsub_pd(eval_3, level, index);
#endif
#if  ( MIC_UNROLLING_WIDTH >  32 )
                eval_4 = _mm512_fmsub_pd(eval_4, level, index);
                eval_5 = _mm512_fmsub_pd(eval_5, level, index);
#endif
#if  ( MIC_UNROLLING_WIDTH >  48 )
                eval_6 = _mm512_fmsub_pd(eval_6, level, index);
                eval_7 = _mm512_fmsub_pd(eval_7, level, index);
#endif
#if  ( MIC_UNROLLING_WIDTH >  64 )
                eval_8 = _mm512_fmsub_pd(eval_8, level, index);
#endif
#if  ( MIC_UNROLLING_WIDTH >  72 )
                eval_9 = _mm512_fmsub_pd(eval_9, level, index);
#endif
#if  ( MIC_UNROLLING_WIDTH >  80 )
                eval_10 = _mm512_fmsub_pd(eval_10, level, index);
                eval_11 = _mm512_fmsub_pd(eval_11, level, index);
#endif

                __m512d one = _mm512_set_1to8_pd(1.0);
                __m512d zero = _mm512_set_1to8_pd(0.0);
                __m512i abs2MaskLRBni = _mm512_set_1to8_epi64(0x7FFFFFFFFFFFFFFF);

                eval_0 = _mm512_castsi512_pd(_mm512_and_epi64( abs2MaskLRBni, _mm512_castpd_si512(eval_0)));
                eval_1 = _mm512_castsi512_pd(_mm512_and_epi64( abs2MaskLRBni, _mm512_castpd_si512(eval_1)));
                eval_2 = _mm512_castsi512_pd(_mm512_and_epi64( abs2MaskLRBni, _mm512_castpd_si512(eval_2)));
#if  ( MIC_UNROLLING_WIDTH >  24 )
                eval_3 = _mm512_castsi512_pd(_mm512_and_epi64( abs2MaskLRBni, _mm512_castpd_si512(eval_3)));
#endif
#if  ( MIC_UNROLLING_WIDTH >  32 )
                eval_4 = _mm512_castsi512_pd(_mm512_and_epi64( abs2MaskLRBni, _mm512_castpd_si512(eval_4)));
                eval_5 = _mm512_castsi512_pd(_mm512_and_epi64( abs2MaskLRBni, _mm512_castpd_si512(eval_5)));
#endif
#if  ( MIC_UNROLLING_WIDTH >  48 )
                eval_6 = _mm512_castsi512_pd(_mm512_and_epi64( abs2MaskLRBni, _mm512_castpd_si512(eval_6)));
                eval_7 = _mm512_castsi512_pd(_mm512_and_epi64( abs2MaskLRBni, _mm512_castpd_si512(eval_7)));
#endif
#if  ( MIC_UNROLLING_WIDTH >  64 )
                eval_8 = _mm512_castsi512_pd(_mm512_and_epi64( abs2MaskLRBni, _mm512_castpd_si512(eval_8)));
#endif
#if  ( MIC_UNROLLING_WIDTH >  72 )
                eval_9 = _mm512_castsi512_pd(_mm512_and_epi64( abs2MaskLRBni, _mm512_castpd_si512(eval_9)));
#endif
#if  ( MIC_UNROLLING_WIDTH >  80 )
                eval_10 = _mm512_castsi512_pd(_mm512_and_epi64( abs2MaskLRBni, _mm512_castpd_si512(eval_10)));
                eval_11 = _mm512_castsi512_pd(_mm512_and_epi64( abs2MaskLRBni, _mm512_castpd_si512(eval_11)));
#endif

                eval_0 = _mm512_sub_pd(one, eval_0);
                eval_1 = _mm512_sub_pd(one, eval_1);
                eval_2 = _mm512_sub_pd(one, eval_2);
#if  ( MIC_UNROLLING_WIDTH >  24 )
                eval_3 = _mm512_sub_pd(one, eval_3);
#endif
#if  ( MIC_UNROLLING_WIDTH >  32 )
                eval_4 = _mm512_sub_pd(one, eval_4);
                eval_5 = _mm512_sub_pd(one, eval_5);
#endif
#if  ( MIC_UNROLLING_WIDTH >  48 )
                eval_6 = _mm512_sub_pd(one, eval_6);
                eval_7 = _mm512_sub_pd(one, eval_7);
#endif
#if  ( MIC_UNROLLING_WIDTH >  64 )
                eval_8 = _mm512_sub_pd(one, eval_8);
#endif
#if  ( MIC_UNROLLING_WIDTH >  72 )
                eval_9 = _mm512_sub_pd(one, eval_9);
#endif
#if  ( MIC_UNROLLING_WIDTH >  80 )
                eval_10 = _mm512_sub_pd(one, eval_10);
                eval_11 = _mm512_sub_pd(one, eval_11);
#endif
                eval_0 = _mm512_gmax_pd(zero, eval_0);
                eval_1 = _mm512_gmax_pd(zero, eval_1);
                eval_2 = _mm512_gmax_pd(zero, eval_2);
#if  ( MIC_UNROLLING_WIDTH >  24 )
                eval_3 = _mm512_gmax_pd(zero, eval_3);
#endif
#if  ( MIC_UNROLLING_WIDTH >  32 )
                eval_4 = _mm512_gmax_pd(zero, eval_4);
                eval_5 = _mm512_gmax_pd(zero, eval_5);
#endif
#if  ( MIC_UNROLLING_WIDTH >  48 )
                eval_6 = _mm512_gmax_pd(zero, eval_6);
                eval_7 = _mm512_gmax_pd(zero, eval_7);
#endif
#if  ( MIC_UNROLLING_WIDTH >  64 )
                eval_8 = _mm512_gmax_pd(zero, eval_8);
#endif
#if  ( MIC_UNROLLING_WIDTH >  72 )
                eval_9 = _mm512_gmax_pd(zero, eval_9);
#endif
#if  ( MIC_UNROLLING_WIDTH >  80 )
                eval_10 = _mm512_gmax_pd(zero, eval_10);
                eval_11 = _mm512_gmax_pd(zero, eval_11);
#endif
                support_0 = _mm512_mul_pd(support_0, eval_0);
                support_1 = _mm512_mul_pd(support_1, eval_1);
                support_2 = _mm512_mul_pd(support_2, eval_2);
#if  ( MIC_UNROLLING_WIDTH >  24 )
                support_3 = _mm512_mul_pd(support_3, eval_3);
#endif
#if  ( MIC_UNROLLING_WIDTH >  32 )
                support_4 = _mm512_mul_pd(support_4, eval_4);
                support_5 = _mm512_mul_pd(support_5, eval_5);
#endif
#if  ( MIC_UNROLLING_WIDTH >  48 )
                support_6 = _mm512_mul_pd(support_6, eval_6);
                support_7 = _mm512_mul_pd(support_7, eval_7);
#endif
#if  ( MIC_UNROLLING_WIDTH >  64 )
                support_8 = _mm512_mul_pd(support_8, eval_8);
#endif
#if  ( MIC_UNROLLING_WIDTH >  72 )
                support_9 = _mm512_mul_pd(support_9, eval_9);
#endif
#if  ( MIC_UNROLLING_WIDTH >  80 )
                support_10 = _mm512_mul_pd(support_10, eval_10);
                support_11 = _mm512_mul_pd(support_11, eval_11);
#endif
              }

              support_0 = _mm512_add_pd(support_0, support_1);
#if  ( MIC_UNROLLING_WIDTH >  24 )
              support_2 = _mm512_add_pd(support_2, support_3);
#endif
#if  ( MIC_UNROLLING_WIDTH >  32 )
              support_4 = _mm512_add_pd(support_4, support_5);
#endif
#if  ( MIC_UNROLLING_WIDTH >  48 )
              support_6 = _mm512_add_pd(support_6, support_7);
#endif
#if  ( MIC_UNROLLING_WIDTH >  72 )
              support_8 = _mm512_add_pd(support_8, support_9);
#endif
#if  ( MIC_UNROLLING_WIDTH >  80 )
              support_10 = _mm512_add_pd(support_10, support_11);
#endif

              support_0 = _mm512_add_pd(support_0, support_2);
#if  ( MIC_UNROLLING_WIDTH >  48 )
              support_4 = _mm512_add_pd(support_4, support_6);
#endif
#if  ( MIC_UNROLLING_WIDTH >  80 )
              support_8 = _mm512_add_pd(support_8, support_10);
#endif

#if  ( MIC_UNROLLING_WIDTH >  32 )
              support_0 = _mm512_add_pd(support_0, support_4);
#endif
#if  ( MIC_UNROLLING_WIDTH >  64 )
              support_0 = _mm512_add_pd(support_0, support_8);
#endif
              ptrResult[j] += _mm512_reduce_add_pd(support_0);
            }

            _mm_prefetch((const char*) & (ptrSource[i + MIC_UNROLLING_WIDTH + 0]), _MM_HINT_T0);
            _mm_prefetch((const char*) & (ptrSource[i + MIC_UNROLLING_WIDTH + 8]), _MM_HINT_T0);
            _mm_prefetch((const char*) & (ptrSource[i + MIC_UNROLLING_WIDTH + 16]), _MM_HINT_T0);
#if  ( MIC_UNROLLING_WIDTH >  24 )
            _mm_prefetch((const char*) & (ptrSource[i + MIC_UNROLLING_WIDTH + 24]), _MM_HINT_T0);
#endif
#if  ( MIC_UNROLLING_WIDTH >  32 )
            _mm_prefetch((const char*) & (ptrSource[i + MIC_UNROLLING_WIDTH + 32]), _MM_HINT_T0);
            _mm_prefetch((const char*) & (ptrSource[i + MIC_UNROLLING_WIDTH + 40]), _MM_HINT_T0);
#endif
#if  ( MIC_UNROLLING_WIDTH >  48 )
            _mm_prefetch((const char*) & (ptrSource[i + MIC_UNROLLING_WIDTH + 48]), _MM_HINT_T0);
            _mm_prefetch((const char*) & (ptrSource[i + MIC_UNROLLING_WIDTH + 56]), _MM_HINT_T0);
#endif
#if  ( MIC_UNROLLING_WIDTH >  64 )
            _mm_prefetch((const char*) & (ptrSource[i + MIC_UNROLLING_WIDTH + 64]), _MM_HINT_T0);
#endif
#if  ( MIC_UNROLLING_WIDTH >  72 )
            _mm_prefetch((const char*) & (ptrSource[i + MIC_UNROLLING_WIDTH + 72]), _MM_HINT_T0);
#endif
#if  ( MIC_UNROLLING_WIDTH >  80 )
            _mm_prefetch((const char*) & (ptrSource[i + MIC_UNROLLING_WIDTH + 80]), _MM_HINT_T0);
            _mm_prefetch((const char*) & (ptrSource[i + MIC_UNROLLING_WIDTH + 88]), _MM_HINT_T0);
#endif

            __m512d support_0 = _mm512_load_pd(&(ptrSource[i + 0]));
            __m512d support_1 = _mm512_load_pd(&(ptrSource[i + 8]));
            __m512d support_2 = _mm512_load_pd(&(ptrSource[i + 16]));
#if  ( MIC_UNROLLING_WIDTH >  24 )
            __m512d support_3 = _mm512_load_pd(&(ptrSource[i + 24]));
#endif
#if  ( MIC_UNROLLING_WIDTH >  32 )
            __m512d support_4 = _mm512_load_pd(&(ptrSource[i + 32]));
            __m512d support_5 = _mm512_load_pd(&(ptrSource[i + 40]));
#endif
#if  ( MIC_UNROLLING_WIDTH >  48 )
            __m512d support_6 = _mm512_load_pd(&(ptrSource[i + 48]));
            __m512d support_7 = _mm512_load_pd(&(ptrSource[i + 56]));
#endif
#if  ( MIC_UNROLLING_WIDTH >  64 )
            __m512d support_8 = _mm512_load_pd(&(ptrSource[i + 64]));
#endif
#if  ( MIC_UNROLLING_WIDTH >  72 )
            __m512d support_9 = _mm512_load_pd(&(ptrSource[i + 72]));
#endif
#if  ( MIC_UNROLLING_WIDTH >  80 )
            __m512d support_10 = _mm512_load_pd(&(ptrSource[i + 80]));
            __m512d support_11 = _mm512_load_pd(&(ptrSource[i + 88]));
#endif

            size_t j = end_index_grid - 1;

            for (size_t d = 0; d < dims; d++) {
              __m512d eval_0 = _mm512_load_pd(&(ptrData[(d * sourceSize) + i + 0]));
              __m512d eval_1 = _mm512_load_pd(&(ptrData[(d * sourceSize) + i + 8]));
              __m512d eval_2 = _mm512_load_pd(&(ptrData[(d * sourceSize) + i + 16]));
#if  ( MIC_UNROLLING_WIDTH >  24 )
              __m512d eval_3 = _mm512_load_pd(&(ptrData[(d * sourceSize) + i + 24]));
#endif
#if  ( MIC_UNROLLING_WIDTH >  32 )
              __m512d eval_4 = _mm512_load_pd(&(ptrData[(d * sourceSize) + i + 32]));
              __m512d eval_5 = _mm512_load_pd(&(ptrData[(d * sourceSize) + i + 40]));
#endif
#if  ( MIC_UNROLLING_WIDTH >  48 )
              __m512d eval_6 = _mm512_load_pd(&(ptrData[(d * sourceSize) + i + 48]));
              __m512d eval_7 = _mm512_load_pd(&(ptrData[(d * sourceSize) + i + 56]));
#endif
#if  ( MIC_UNROLLING_WIDTH >  64 )
              __m512d eval_8 = _mm512_load_pd(&(ptrData[(d * sourceSize) + i + 64]));
#endif
#if  ( MIC_UNROLLING_WIDTH >  72 )
              __m512d eval_9 = _mm512_load_pd(&(ptrData[(d * sourceSize) + i + 72]));
#endif
#if  ( MIC_UNROLLING_WIDTH >  80 )
              __m512d eval_10 = _mm512_load_pd(&(ptrData[(d * sourceSize) + i + 80]));
              __m512d eval_11 = _mm512_load_pd(&(ptrData[(d * sourceSize) + i + 88]));
#endif

              _mm_prefetch((const char*) & (ptrData[(d * sourceSize) + i + MIC_UNROLLING_WIDTH + 0]), _MM_HINT_T0);
              _mm_prefetch((const char*) & (ptrData[(d * sourceSize) + i + MIC_UNROLLING_WIDTH + 8]), _MM_HINT_T0);
              _mm_prefetch((const char*) & (ptrData[(d * sourceSize) + i + MIC_UNROLLING_WIDTH + 16]), _MM_HINT_T0);

              __m512d level = _mm512_extload_pd(&(ptrLevel[(j * dims) + d]),  _MM_UPCONV_PD_NONE, _MM_BROADCAST_1X8, _MM_HINT_NONE);
              __m512d index = _mm512_extload_pd(&(ptrIndex[(j * dims) + d]),  _MM_UPCONV_PD_NONE, _MM_BROADCAST_1X8, _MM_HINT_NONE);

              eval_0 = _mm512_fmsub_pd(eval_0, level, index);
              eval_1 = _mm512_fmsub_pd(eval_1, level, index);
              eval_2 = _mm512_fmsub_pd(eval_2, level, index);
#if  ( MIC_UNROLLING_WIDTH >  24 )
              eval_3 = _mm512_fmsub_pd(eval_3, level, index);
#endif
#if  ( MIC_UNROLLING_WIDTH >  32 )
              eval_4 = _mm512_fmsub_pd(eval_4, level, index);
              eval_5 = _mm512_fmsub_pd(eval_5, level, index);
#endif
#if  ( MIC_UNROLLING_WIDTH >  48 )
              eval_6 = _mm512_fmsub_pd(eval_6, level, index);
              eval_7 = _mm512_fmsub_pd(eval_7, level, index);
#endif
#if  ( MIC_UNROLLING_WIDTH >  64 )
              eval_8 = _mm512_fmsub_pd(eval_8, level, index);
#endif
#if  ( MIC_UNROLLING_WIDTH >  72 )
              eval_9 = _mm512_fmsub_pd(eval_9, level, index);
#endif
#if  ( MIC_UNROLLING_WIDTH >  80 )
              eval_10 = _mm512_fmsub_pd(eval_10, level, index);
              eval_11 = _mm512_fmsub_pd(eval_11, level, index);
#endif

#if  ( MIC_UNROLLING_WIDTH >  24 )
              _mm_prefetch((const char*) & (ptrData[(d * sourceSize) + i + MIC_UNROLLING_WIDTH + 24]), _MM_HINT_T0);
#endif
#if  ( MIC_UNROLLING_WIDTH >  32 )
              _mm_prefetch((const char*) & (ptrData[(d * sourceSize) + i + MIC_UNROLLING_WIDTH + 32]), _MM_HINT_T0);
              _mm_prefetch((const char*) & (ptrData[(d * sourceSize) + i + MIC_UNROLLING_WIDTH + 40]), _MM_HINT_T0);
#endif
              __m512d one = _mm512_set_1to8_pd(1.0);
              __m512d zero = _mm512_set_1to8_pd(0.0);

              __m512i abs2MaskLRBni = _mm512_set_1to8_epi64(0x7FFFFFFFFFFFFFFF);
              eval_0 = _mm512_castsi512_pd(_mm512_and_epi64( abs2MaskLRBni, _mm512_castpd_si512(eval_0)));
              eval_1 = _mm512_castsi512_pd(_mm512_and_epi64( abs2MaskLRBni, _mm512_castpd_si512(eval_1)));
              eval_2 = _mm512_castsi512_pd(_mm512_and_epi64( abs2MaskLRBni, _mm512_castpd_si512(eval_2)));
#if  ( MIC_UNROLLING_WIDTH >  24 )
              eval_3 = _mm512_castsi512_pd(_mm512_and_epi64( abs2MaskLRBni, _mm512_castpd_si512(eval_3)));
#endif
#if  ( MIC_UNROLLING_WIDTH >  32 )
              eval_4 = _mm512_castsi512_pd(_mm512_and_epi64( abs2MaskLRBni, _mm512_castpd_si512(eval_4)));
              eval_5 = _mm512_castsi512_pd(_mm512_and_epi64( abs2MaskLRBni, _mm512_castpd_si512(eval_5)));
#endif
#if  ( MIC_UNROLLING_WIDTH >  48 )
              eval_6 = _mm512_castsi512_pd(_mm512_and_epi64( abs2MaskLRBni, _mm512_castpd_si512(eval_6)));
              eval_7 = _mm512_castsi512_pd(_mm512_and_epi64( abs2MaskLRBni, _mm512_castpd_si512(eval_7)));
#endif
#if  ( MIC_UNROLLING_WIDTH >  64 )
              eval_8 = _mm512_castsi512_pd(_mm512_and_epi64( abs2MaskLRBni, _mm512_castpd_si512(eval_8)));
#endif
#if  ( MIC_UNROLLING_WIDTH >  72 )
              eval_9 = _mm512_castsi512_pd(_mm512_and_epi64( abs2MaskLRBni, _mm512_castpd_si512(eval_9)));
#endif
#if  ( MIC_UNROLLING_WIDTH >  80 )
              eval_10 = _mm512_castsi512_pd(_mm512_and_epi64( abs2MaskLRBni, _mm512_castpd_si512(eval_10)));
              eval_11 = _mm512_castsi512_pd(_mm512_and_epi64( abs2MaskLRBni, _mm512_castpd_si512(eval_11)));
#endif
#if  ( MIC_UNROLLING_WIDTH >  48 )
              _mm_prefetch((const char*) & (ptrData[(d * sourceSize) + i + MIC_UNROLLING_WIDTH + 48]), _MM_HINT_T0);
              _mm_prefetch((const char*) & (ptrData[(d * sourceSize) + i + MIC_UNROLLING_WIDTH + 56]), _MM_HINT_T0);
#endif
#if  ( MIC_UNROLLING_WIDTH >  64 )
              _mm_prefetch((const char*) & (ptrData[(d * sourceSize) + i + MIC_UNROLLING_WIDTH + 64]), _MM_HINT_T0);
#endif

              eval_0 = _mm512_sub_pd(one, eval_0);
              eval_1 = _mm512_sub_pd(one, eval_1);
              eval_2 = _mm512_sub_pd(one, eval_2);
#if  ( MIC_UNROLLING_WIDTH >  24 )
              eval_3 = _mm512_sub_pd(one, eval_3);
#endif
#if  ( MIC_UNROLLING_WIDTH >  32 )
              eval_4 = _mm512_sub_pd(one, eval_4);
              eval_5 = _mm512_sub_pd(one, eval_5);
#endif
#if  ( MIC_UNROLLING_WIDTH >  48 )
              eval_6 = _mm512_sub_pd(one, eval_6);
              eval_7 = _mm512_sub_pd(one, eval_7);
#endif
#if  ( MIC_UNROLLING_WIDTH >  64 )
              eval_8 = _mm512_sub_pd(one, eval_8);
#endif
#if  ( MIC_UNROLLING_WIDTH >  72 )
              eval_9 = _mm512_sub_pd(one, eval_9);
#endif
#if  ( MIC_UNROLLING_WIDTH >  80 )
              eval_10 = _mm512_sub_pd(one, eval_10);
              eval_11 = _mm512_sub_pd(one, eval_11);
#endif

#if  ( MIC_UNROLLING_WIDTH >  72 )
              _mm_prefetch((const char*) & (ptrData[(d * sourceSize) + i + MIC_UNROLLING_WIDTH + 72]), _MM_HINT_T0);
#endif
#if  ( MIC_UNROLLING_WIDTH >  80 )
              _mm_prefetch((const char*) & (ptrData[(d * sourceSize) + i + MIC_UNROLLING_WIDTH + 80]), _MM_HINT_T0);
              _mm_prefetch((const char*) & (ptrData[(d * sourceSize) + i + MIC_UNROLLING_WIDTH + 88]), _MM_HINT_T0);

#endif
              eval_0 = _mm512_gmax_pd(zero, eval_0);
              eval_1 = _mm512_gmax_pd(zero, eval_1);
              eval_2 = _mm512_gmax_pd(zero, eval_2);
#if  ( MIC_UNROLLING_WIDTH >  24 )
              eval_3 = _mm512_gmax_pd(zero, eval_3);
#endif
#if  ( MIC_UNROLLING_WIDTH >  32 )
              eval_4 = _mm512_gmax_pd(zero, eval_4);
              eval_5 = _mm512_gmax_pd(zero, eval_5);
#endif
#if  ( MIC_UNROLLING_WIDTH >  48 )
              eval_6 = _mm512_gmax_pd(zero, eval_6);
              eval_7 = _mm512_gmax_pd(zero, eval_7);
#endif
#if  ( MIC_UNROLLING_WIDTH >  64 )
              eval_8 = _mm512_gmax_pd(zero, eval_8);
#endif
#if  ( MIC_UNROLLING_WIDTH >  72 )
              eval_9 = _mm512_gmax_pd(zero, eval_9);
#endif
#if  ( MIC_UNROLLING_WIDTH >  80 )
              eval_10 = _mm512_gmax_pd(zero, eval_10);
              eval_11 = _mm512_gmax_pd(zero, eval_11);
#endif

              support_0 = _mm512_mul_pd(support_0, eval_0);
              support_1 = _mm512_mul_pd(support_1, eval_1);
              support_2 = _mm512_mul_pd(support_2, eval_2);
#if  ( MIC_UNROLLING_WIDTH >  24 )
              support_3 = _mm512_mul_pd(support_3, eval_3);
#endif
#if  ( MIC_UNROLLING_WIDTH >  32 )
              support_4 = _mm512_mul_pd(support_4, eval_4);
              support_5 = _mm512_mul_pd(support_5, eval_5);
#endif
#if  ( MIC_UNROLLING_WIDTH >  48 )
              support_6 = _mm512_mul_pd(support_6, eval_6);
              support_7 = _mm512_mul_pd(support_7, eval_7);
#endif
#if  ( MIC_UNROLLING_WIDTH >  64 )
              support_8 = _mm512_mul_pd(support_8, eval_8);
#endif
#if  ( MIC_UNROLLING_WIDTH >  72 )
              support_9 = _mm512_mul_pd(support_9, eval_9);
#endif
#if  ( MIC_UNROLLING_WIDTH >  80 )
              support_10 = _mm512_mul_pd(support_10, eval_10);
              support_11 = _mm512_mul_pd(support_11, eval_11);
#endif
            }

            support_0 = _mm512_add_pd(support_0, support_1);
#if  ( MIC_UNROLLING_WIDTH >  24 )
            support_2 = _mm512_add_pd(support_2, support_3);
#endif
#if  ( MIC_UNROLLING_WIDTH >  32 )
            support_4 = _mm512_add_pd(support_4, support_5);
#endif
#if  ( MIC_UNROLLING_WIDTH >  48 )
            support_6 = _mm512_add_pd(support_6, support_7);
#endif
#if  ( MIC_UNROLLING_WIDTH >  72 )
            support_8 = _mm512_add_pd(support_8, support_9);
#endif
#if  ( MIC_UNROLLING_WIDTH >  80 )
            support_10 = _mm512_add_pd(support_10, support_11);
#endif

            support_0 = _mm512_add_pd(support_0, support_2);
#if  ( MIC_UNROLLING_WIDTH >  48 )
            support_4 = _mm512_add_pd(support_4, support_6);
#endif
#if  ( MIC_UNROLLING_WIDTH >  80 )
            support_8 = _mm512_add_pd(support_8, support_10);
#endif

#if  ( MIC_UNROLLING_WIDTH >  32 )
            support_0 = _mm512_add_pd(support_0, support_4);
#endif
#if  ( MIC_UNROLLING_WIDTH >  64 )
            support_0 = _mm512_add_pd(support_0, support_8);
#endif

            ptrResult[j] += _mm512_reduce_add_pd(support_0);
          }

#else

          for (size_t i = start_index_data; i < end_index_data; i++) {
            for (size_t j = start_index_grid; j < end_index_grid; j++) {
              double curSupport = ptrSource[i];

              for (size_t d = 0; d < dims; d++) {
                double eval = ((ptrLevel[(j * dims) + d]) * (ptrData[(d * sourceSize) + i]));
                double index_calc = eval - (ptrIndex[(j * dims) + d]);
                double abs = fabs(index_calc);
                double last = 1.0 - abs;
                double localSupport = std::max<double>(last, 0.0);
                curSupport *= localSupport;
              }

              ptrResult[j] += curSupport;
            }
          }

#endif
        }
    };
  }
}
#endif // USEMIC

#endif // MICLINEAR_HPP
