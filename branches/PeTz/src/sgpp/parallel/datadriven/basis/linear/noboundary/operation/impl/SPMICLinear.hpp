/* ****************************************************************************
* Copyright (C) 2012 - 2013 Technische Universitaet Muenchen                  *
* This file is part of the SG++ project. For conditions of distribution and   *
* use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
**************************************************************************** */
// @author Alexander Heinecke (Alexander.Heinecke@mytum.de)
// @author Roman Karlstetter (karlstetter@mytum.de)

#ifndef SPMICLINEAR_HPP
#define SPMICLINEAR_HPP
#ifdef USEMIC

#ifdef __MIC__
#include <immintrin.h>
#endif

#include "parallel/datadriven/basis/common/mic/MICKernelBase.hpp"

namespace sg {
  namespace parallel {
    /**
     * @brief The MICLinear class
     *
     * multImpl and multTransImpl have to be run on the mic, either in offload mode or native
     */
    class SPMICLinear : public SPMICKernelBase {
      public:
        static const KernelType kernelType = Standard;
        static inline void multImpl(
          float* ptrLevel,
          float* ptrIndex,
          float* /*ptrMask*/, //unused for this specialization
          float* /*ptrOffset*/, //unused for this specialization
          float* ptrData,
          float* ptrAlpha,
          float* ptrResult,
          size_t result_size,
          size_t dims,
          const size_t start_index_grid,
          const size_t end_index_grid,
          const size_t start_index_data,
          const size_t end_index_data) {
          for (size_t i = start_index_data; i < end_index_data; i += getChunkDataPoints()) {
#ifdef __MIC__

            for (size_t j = start_index_grid; j < end_index_grid; j++) {
              __m512 support_0 = _mm512_extload_ps(&(ptrAlpha[j]), _MM_UPCONV_PS_NONE, _MM_BROADCAST_1X16, _MM_HINT_NONE);
              __m512 support_1 = _mm512_extload_ps(&(ptrAlpha[j]), _MM_UPCONV_PS_NONE, _MM_BROADCAST_1X16, _MM_HINT_NONE);
              __m512 support_2 = _mm512_extload_ps(&(ptrAlpha[j]), _MM_UPCONV_PS_NONE, _MM_BROADCAST_1X16, _MM_HINT_NONE);
#if  ( MIC_UNROLLING_WIDTH_SP >  48 )
              __m512 support_3 = _mm512_extload_ps(&(ptrAlpha[j]), _MM_UPCONV_PS_NONE, _MM_BROADCAST_1X16, _MM_HINT_NONE);
#endif
#if  ( MIC_UNROLLING_WIDTH_SP >  64 )
              __m512 support_4 = _mm512_extload_ps(&(ptrAlpha[j]), _MM_UPCONV_PS_NONE, _MM_BROADCAST_1X16, _MM_HINT_NONE);
              __m512 support_5 = _mm512_extload_ps(&(ptrAlpha[j]), _MM_UPCONV_PS_NONE, _MM_BROADCAST_1X16, _MM_HINT_NONE);
#endif
#if  ( MIC_UNROLLING_WIDTH_SP >  96 )
              __m512 support_6 = _mm512_extload_ps(&(ptrAlpha[j]), _MM_UPCONV_PS_NONE, _MM_BROADCAST_1X16, _MM_HINT_NONE);
              __m512 support_7 = _mm512_extload_ps(&(ptrAlpha[j]), _MM_UPCONV_PS_NONE, _MM_BROADCAST_1X16, _MM_HINT_NONE);
#endif
#if  ( MIC_UNROLLING_WIDTH_SP >  128 )
              __m512 support_8 = _mm512_extload_ps(&(ptrAlpha[j]), _MM_UPCONV_PS_NONE, _MM_BROADCAST_1X16, _MM_HINT_NONE);
#endif
#if  ( MIC_UNROLLING_WIDTH_SP >  144 )
              __m512 support_9 = _mm512_extload_ps(&(ptrAlpha[j]), _MM_UPCONV_PS_NONE, _MM_BROADCAST_1X16, _MM_HINT_NONE);
#endif
#if  ( MIC_UNROLLING_WIDTH_SP >  160 )
              __m512 support_10 = _mm512_extload_ps(&(ptrAlpha[j]), _MM_UPCONV_PS_NONE, _MM_BROADCAST_1X16, _MM_HINT_NONE);
              __m512 support_11 = _mm512_extload_ps(&(ptrAlpha[j]), _MM_UPCONV_PS_NONE, _MM_BROADCAST_1X16, _MM_HINT_NONE);
#endif

              for (size_t d = 0; d < dims - 1; d++) {
                __m512 eval_0 = _mm512_load_ps(&(ptrData[(d * result_size) + i + 0]));
                __m512 eval_1 = _mm512_load_ps(&(ptrData[(d * result_size) + i + 16]));
                __m512 eval_2 = _mm512_load_ps(&(ptrData[(d * result_size) + i + 32]));
#if  ( MIC_UNROLLING_WIDTH_SP >  48 )
                __m512 eval_3 = _mm512_load_ps(&(ptrData[(d * result_size) + i + 48]));
#endif
#if  ( MIC_UNROLLING_WIDTH_SP >  64 )
                __m512 eval_4 = _mm512_load_ps(&(ptrData[(d * result_size) + i + 64]));
                __m512 eval_5 = _mm512_load_ps(&(ptrData[(d * result_size) + i + 80]));
#endif
#if  ( MIC_UNROLLING_WIDTH_SP >  96 )
                __m512 eval_6 = _mm512_load_ps(&(ptrData[(d * result_size) + i + 96]));
                __m512 eval_7 = _mm512_load_ps(&(ptrData[(d * result_size) + i + 112]));
#endif
#if  ( MIC_UNROLLING_WIDTH_SP >  128 )
                __m512 eval_8 = _mm512_load_ps(&(ptrData[(d * result_size) + i + 128]));
#endif
#if  ( MIC_UNROLLING_WIDTH_SP >  144 )
                __m512 eval_9 = _mm512_load_ps(&(ptrData[(d * result_size) + i + 144]));
#endif
#if  ( MIC_UNROLLING_WIDTH_SP >  160 )
                __m512 eval_10 = _mm512_load_ps(&(ptrData[(d * result_size) + i + 160]));
                __m512 eval_11 = _mm512_load_ps(&(ptrData[(d * result_size) + i + 176]));
#endif
                __m512 level = _mm512_extload_ps(&(ptrLevel[(j * dims) + d]), _MM_UPCONV_PS_NONE, _MM_BROADCAST_1X16, _MM_HINT_NONE);
                __m512 index = _mm512_extload_ps(&(ptrIndex[(j * dims) + d]), _MM_UPCONV_PS_NONE, _MM_BROADCAST_1X16, _MM_HINT_NONE);

                eval_0 = _mm512_fmsub_ps(eval_0, level, index);
                eval_1 = _mm512_fmsub_ps(eval_1, level, index);
                eval_2 = _mm512_fmsub_ps(eval_2, level, index);
#if  ( MIC_UNROLLING_WIDTH_SP >  48 )
                eval_3 = _mm512_fmsub_ps(eval_3, level, index);
#endif
#if  ( MIC_UNROLLING_WIDTH_SP >  64 )
                eval_4 = _mm512_fmsub_ps(eval_4, level, index);
                eval_5 = _mm512_fmsub_ps(eval_5, level, index);
#endif
#if  ( MIC_UNROLLING_WIDTH_SP >  96 )
                eval_6 = _mm512_fmsub_ps(eval_6, level, index);
                eval_7 = _mm512_fmsub_ps(eval_7, level, index);
#endif
#if  ( MIC_UNROLLING_WIDTH_SP >  128 )
                eval_8 = _mm512_fmsub_ps(eval_8, level, index);
#endif
#if  ( MIC_UNROLLING_WIDTH_SP >  144 )
                eval_9 = _mm512_fmsub_ps(eval_9, level, index);
#endif
#if  ( MIC_UNROLLING_WIDTH_SP >  160 )
                eval_10 = _mm512_fmsub_ps(eval_10, level, index);
                eval_11 = _mm512_fmsub_ps(eval_11, level, index);
#endif

                eval_0 = _mm512_gmaxabs_ps( eval_0, eval_0);
                eval_1 = _mm512_gmaxabs_ps( eval_1, eval_1);
                eval_2 = _mm512_gmaxabs_ps( eval_2, eval_2);
#if  ( MIC_UNROLLING_WIDTH_SP >  48 )
                eval_3 = _mm512_gmaxabs_ps( eval_3, eval_3);
#endif
#if  ( MIC_UNROLLING_WIDTH_SP >  64 )
                eval_4 = _mm512_gmaxabs_ps( eval_4, eval_4);
                eval_5 = _mm512_gmaxabs_ps( eval_5, eval_5);
#endif
#if  ( MIC_UNROLLING_WIDTH_SP >  96 )
                eval_6 = _mm512_gmaxabs_ps( eval_6, eval_6);
                eval_7 = _mm512_gmaxabs_ps( eval_7, eval_7);
#endif
#if  ( MIC_UNROLLING_WIDTH_SP >  128 )
                eval_8 = _mm512_gmaxabs_ps( eval_8, eval_8);
#endif
#if  ( MIC_UNROLLING_WIDTH_SP >  144 )
                eval_9 = _mm512_gmaxabs_ps( eval_9, eval_9);
#endif
#if  ( MIC_UNROLLING_WIDTH_SP >  160 )
                eval_10 = _mm512_gmaxabs_ps( eval_10, eval_10);
                eval_11 = _mm512_gmaxabs_ps( eval_11, eval_11);
#endif

                __m512 one = _mm512_set_1to16_ps(1.0f);

                eval_0 = _mm512_sub_ps(one, eval_0);
                eval_1 = _mm512_sub_ps(one, eval_1);
                eval_2 = _mm512_sub_ps(one, eval_2);
#if  ( MIC_UNROLLING_WIDTH_SP >  48 )
                eval_3 = _mm512_sub_ps(one, eval_3);
#endif
#if  ( MIC_UNROLLING_WIDTH_SP >  64 )
                eval_4 = _mm512_sub_ps(one, eval_4);
                eval_5 = _mm512_sub_ps(one, eval_5);
#endif
#if  ( MIC_UNROLLING_WIDTH_SP >  96 )
                eval_6 = _mm512_sub_ps(one, eval_6);
                eval_7 = _mm512_sub_ps(one, eval_7);
#endif
#if  ( MIC_UNROLLING_WIDTH_SP >  128 )
                eval_8 = _mm512_sub_ps(one, eval_8);
#endif
#if  ( MIC_UNROLLING_WIDTH_SP >  144 )
                eval_9 = _mm512_sub_ps(one, eval_9);
#endif
#if  ( MIC_UNROLLING_WIDTH_SP >  160 )
                eval_10 = _mm512_sub_ps(one, eval_10);
                eval_11 = _mm512_sub_ps(one, eval_11);
#endif

                __m512 zero = _mm512_set_1to16_ps(0.0f);

                eval_0 = _mm512_gmax_ps(zero, eval_0);
                eval_1 = _mm512_gmax_ps(zero, eval_1);
                eval_2 = _mm512_gmax_ps(zero, eval_2);
#if  ( MIC_UNROLLING_WIDTH_SP >  48 )
                eval_3 = _mm512_gmax_ps(zero, eval_3);
#endif
#if  ( MIC_UNROLLING_WIDTH_SP >  64 )
                eval_4 = _mm512_gmax_ps(zero, eval_4);
                eval_5 = _mm512_gmax_ps(zero, eval_5);
#endif
#if  ( MIC_UNROLLING_WIDTH_SP >  96 )
                eval_6 = _mm512_gmax_ps(zero, eval_6);
                eval_7 = _mm512_gmax_ps(zero, eval_7);
#endif
#if  ( MIC_UNROLLING_WIDTH_SP >  128 )
                eval_8 = _mm512_gmax_ps(zero, eval_8);
#endif
#if  ( MIC_UNROLLING_WIDTH_SP >  144 )
                eval_9 = _mm512_gmax_ps(zero, eval_9);
#endif
#if  ( MIC_UNROLLING_WIDTH_SP >  160 )
                eval_10 = _mm512_gmax_ps(zero, eval_10);
                eval_11 = _mm512_gmax_ps(zero, eval_11);
#endif

                support_0 = _mm512_mul_ps(support_0, eval_0);
                support_1 = _mm512_mul_ps(support_1, eval_1);
                support_2 = _mm512_mul_ps(support_2, eval_2);
#if  ( MIC_UNROLLING_WIDTH_SP >  48 )
                support_3 = _mm512_mul_ps(support_3, eval_3);
#endif
#if  ( MIC_UNROLLING_WIDTH_SP >  64 )
                support_4 = _mm512_mul_ps(support_4, eval_4);
                support_5 = _mm512_mul_ps(support_5, eval_5);
#endif
#if  ( MIC_UNROLLING_WIDTH_SP >  96 )
                support_6 = _mm512_mul_ps(support_6, eval_6);
                support_7 = _mm512_mul_ps(support_7, eval_7);
#endif
#if  ( MIC_UNROLLING_WIDTH_SP >  128 )
                support_8 = _mm512_mul_ps(support_8, eval_8);
#endif
#if  ( MIC_UNROLLING_WIDTH_SP >  144 )
                support_9 = _mm512_mul_ps(support_9, eval_9);
#endif
#if  ( MIC_UNROLLING_WIDTH_SP >  160 )
                support_10 = _mm512_mul_ps(support_10, eval_10);
                support_11 = _mm512_mul_ps(support_11, eval_11);
#endif
              }

              size_t d = dims - 1;

              __m512 eval_0 = _mm512_load_ps(&(ptrData[(d * result_size) + i + 0]));
              __m512 eval_1 = _mm512_load_ps(&(ptrData[(d * result_size) + i + 16]));
              __m512 eval_2 = _mm512_load_ps(&(ptrData[(d * result_size) + i + 32]));
#if  ( MIC_UNROLLING_WIDTH_SP >  48 )
              __m512 eval_3 = _mm512_load_ps(&(ptrData[(d * result_size) + i + 48]));
#endif
#if  ( MIC_UNROLLING_WIDTH_SP >  64 )
              __m512 eval_4 = _mm512_load_ps(&(ptrData[(d * result_size) + i + 64]));
              __m512 eval_5 = _mm512_load_ps(&(ptrData[(d * result_size) + i + 80]));
#endif
#if  ( MIC_UNROLLING_WIDTH_SP >  96 )
              __m512 eval_6 = _mm512_load_ps(&(ptrData[(d * result_size) + i + 96]));
              __m512 eval_7 = _mm512_load_ps(&(ptrData[(d * result_size) + i + 112]));
#endif
#if  ( MIC_UNROLLING_WIDTH_SP >  128 )
              __m512 eval_8 = _mm512_load_ps(&(ptrData[(d * result_size) + i + 128]));
#endif
#if  ( MIC_UNROLLING_WIDTH_SP >  144 )
              __m512 eval_9 = _mm512_load_ps(&(ptrData[(d * result_size) + i + 144]));
#endif
#if  ( MIC_UNROLLING_WIDTH_SP >  160 )
              __m512 eval_10 = _mm512_load_ps(&(ptrData[(d * result_size) + i + 160]));
              __m512 eval_11 = _mm512_load_ps(&(ptrData[(d * result_size) + i + 176]));
#endif

              __m512 level = _mm512_extload_ps(&(ptrLevel[(j * dims) + d]), _MM_UPCONV_PS_NONE, _MM_BROADCAST_1X16, _MM_HINT_NONE);
              __m512 index = _mm512_extload_ps(&(ptrIndex[(j * dims) + d]), _MM_UPCONV_PS_NONE, _MM_BROADCAST_1X16, _MM_HINT_NONE);

              eval_0 = _mm512_fmsub_ps(eval_0, level, index);
              eval_1 = _mm512_fmsub_ps(eval_1, level, index);
              eval_2 = _mm512_fmsub_ps(eval_2, level, index);
#if  ( MIC_UNROLLING_WIDTH_SP >  48 )
              eval_3 = _mm512_fmsub_ps(eval_3, level, index);
#endif
#if  ( MIC_UNROLLING_WIDTH_SP >  64 )
              eval_4 = _mm512_fmsub_ps(eval_4, level, index);
              eval_5 = _mm512_fmsub_ps(eval_5, level, index);
#endif
#if  ( MIC_UNROLLING_WIDTH_SP >  96 )
              eval_6 = _mm512_fmsub_ps(eval_6, level, index);
              eval_7 = _mm512_fmsub_ps(eval_7, level, index);
#endif
#if  ( MIC_UNROLLING_WIDTH_SP >  128 )
              eval_8 = _mm512_fmsub_ps(eval_8, level, index);
#endif
#if  ( MIC_UNROLLING_WIDTH_SP >  144 )
              eval_9 = _mm512_fmsub_ps(eval_9, level, index);
#endif
#if  ( MIC_UNROLLING_WIDTH_SP >  160 )
              eval_10 = _mm512_fmsub_ps(eval_10, level, index);
              eval_11 = _mm512_fmsub_ps(eval_11, level, index);
#endif

              eval_0 = _mm512_gmaxabs_ps( eval_0, eval_0);
              eval_1 = _mm512_gmaxabs_ps( eval_1, eval_1);
              eval_2 = _mm512_gmaxabs_ps( eval_2, eval_2);
#if  ( MIC_UNROLLING_WIDTH_SP >  48 )
              eval_3 = _mm512_gmaxabs_ps( eval_3, eval_3);
#endif
#if  ( MIC_UNROLLING_WIDTH_SP >  64 )
              eval_4 = _mm512_gmaxabs_ps( eval_4, eval_4);
              eval_5 = _mm512_gmaxabs_ps( eval_5, eval_5);
#endif
#if  ( MIC_UNROLLING_WIDTH_SP >  96 )
              eval_6 = _mm512_gmaxabs_ps( eval_6, eval_6);
              eval_7 = _mm512_gmaxabs_ps( eval_7, eval_7);
#endif
#if  ( MIC_UNROLLING_WIDTH_SP >  128 )
              eval_8 = _mm512_gmaxabs_ps( eval_8, eval_8);
#endif
#if  ( MIC_UNROLLING_WIDTH_SP >  144 )
              eval_9 = _mm512_gmaxabs_ps( eval_9, eval_9);
#endif
#if  ( MIC_UNROLLING_WIDTH_SP >  160 )
              eval_10 = _mm512_gmaxabs_ps( eval_10, eval_10);
              eval_11 = _mm512_gmaxabs_ps( eval_11, eval_11);
#endif

              __m512 one = _mm512_set_1to16_ps(1.0f);

              eval_0 = _mm512_sub_ps(one, eval_0);
              eval_1 = _mm512_sub_ps(one, eval_1);
              eval_2 = _mm512_sub_ps(one, eval_2);
#if  ( MIC_UNROLLING_WIDTH_SP >  48 )
              eval_3 = _mm512_sub_ps(one, eval_3);
#endif
#if  ( MIC_UNROLLING_WIDTH_SP >  64 )
              eval_4 = _mm512_sub_ps(one, eval_4);
              eval_5 = _mm512_sub_ps(one, eval_5);
#endif
#if  ( MIC_UNROLLING_WIDTH_SP >  96 )
              eval_6 = _mm512_sub_ps(one, eval_6);
              eval_7 = _mm512_sub_ps(one, eval_7);
#endif
#if  ( MIC_UNROLLING_WIDTH_SP >  128 )
              eval_8 = _mm512_sub_ps(one, eval_8);
#endif
#if  ( MIC_UNROLLING_WIDTH_SP >  144 )
              eval_9 = _mm512_sub_ps(one, eval_9);
#endif
#if  ( MIC_UNROLLING_WIDTH_SP >  160 )
              eval_10 = _mm512_sub_ps(one, eval_10);
              eval_11 = _mm512_sub_ps(one, eval_11);
#endif

              __m512 zero = _mm512_set_1to16_ps(0.0f);

              eval_0 = _mm512_gmax_ps(zero, eval_0);
              eval_1 = _mm512_gmax_ps(zero, eval_1);
              eval_2 = _mm512_gmax_ps(zero, eval_2);
#if  ( MIC_UNROLLING_WIDTH_SP >  48 )
              eval_3 = _mm512_gmax_ps(zero, eval_3);
#endif
#if  ( MIC_UNROLLING_WIDTH_SP >  64 )
              eval_4 = _mm512_gmax_ps(zero, eval_4);
              eval_5 = _mm512_gmax_ps(zero, eval_5);
#endif
#if  ( MIC_UNROLLING_WIDTH_SP >  96 )
              eval_6 = _mm512_gmax_ps(zero, eval_6);
              eval_7 = _mm512_gmax_ps(zero, eval_7);
#endif
#if  ( MIC_UNROLLING_WIDTH_SP >  128 )
              eval_8 = _mm512_gmax_ps(zero, eval_8);
#endif
#if  ( MIC_UNROLLING_WIDTH_SP >  144 )
              eval_9 = _mm512_gmax_ps(zero, eval_9);
#endif
#if  ( MIC_UNROLLING_WIDTH_SP >  160 )
              eval_10 = _mm512_gmax_ps(zero, eval_10);
              eval_11 = _mm512_gmax_ps(zero, eval_11);
#endif

              __m512 res_0 = _mm512_load_ps(&(ptrResult[i + 0]));
              __m512 res_1 = _mm512_load_ps(&(ptrResult[i + 16]));
              __m512 res_2 = _mm512_load_ps(&(ptrResult[i + 32]));
#if  ( MIC_UNROLLING_WIDTH_SP >  48 )
              __m512 res_3 = _mm512_load_ps(&(ptrResult[i + 48]));
#endif
#if  ( MIC_UNROLLING_WIDTH_SP >  64 )
              __m512 res_4 = _mm512_load_ps(&(ptrResult[i + 64]));
              __m512 res_5 = _mm512_load_ps(&(ptrResult[i + 80]));
#endif

              res_0 = _mm512_add_ps(res_0, _mm512_mul_ps(support_0, eval_0));
              res_1 = _mm512_add_ps(res_1, _mm512_mul_ps(support_1, eval_1));
              res_2 = _mm512_add_ps(res_2, _mm512_mul_ps(support_2, eval_2));
#if  ( MIC_UNROLLING_WIDTH_SP >  48 )
              res_3 = _mm512_add_ps(res_3, _mm512_mul_ps(support_3, eval_3));
#endif
#if  ( MIC_UNROLLING_WIDTH_SP >  64 )
              res_4 = _mm512_add_ps(res_4, _mm512_mul_ps(support_4, eval_4));
              res_5 = _mm512_add_ps(res_5, _mm512_mul_ps(support_5, eval_5));
#endif

              _mm512_store_ps(&(ptrResult[i + 0]), res_0);
              _mm512_store_ps(&(ptrResult[i + 16]), res_1);
              _mm512_store_ps(&(ptrResult[i + 32]), res_2);
#if  ( MIC_UNROLLING_WIDTH_SP >  48 )
              _mm512_store_ps(&(ptrResult[i + 48]), res_3);
#endif
#if  ( MIC_UNROLLING_WIDTH_SP >  64 )
              _mm512_store_ps(&(ptrResult[i + 64]), res_4);
              _mm512_store_ps(&(ptrResult[i + 80]), res_5);
#endif

#if  ( MIC_UNROLLING_WIDTH_SP >  96 )
              __m512 res_6 = _mm512_load_ps(&(ptrResult[i + 96]));
              __m512 res_7 = _mm512_load_ps(&(ptrResult[i + 112]));
#endif
#if  ( MIC_UNROLLING_WIDTH_SP >  128 )
              __m512 res_8 = _mm512_load_ps(&(ptrResult[i + 128]));
#endif
#if  ( MIC_UNROLLING_WIDTH_SP >  144 )
              __m512 res_9 = _mm512_load_ps(&(ptrResult[i + 144]));
#endif
#if  ( MIC_UNROLLING_WIDTH_SP >  160 )
              __m512 res_10 = _mm512_load_ps(&(ptrResult[i + 160]));
              __m512 res_11 = _mm512_load_ps(&(ptrResult[i + 176]));
#endif

#if  ( MIC_UNROLLING_WIDTH_SP >  96 )
              res_6 = _mm512_add_ps(res_6, _mm512_mul_ps(support_6, eval_6));
              res_7 = _mm512_add_ps(res_7, _mm512_mul_ps(support_7, eval_7));
#endif
#if  ( MIC_UNROLLING_WIDTH_SP >  128 )
              res_8 = _mm512_add_ps(res_8, _mm512_mul_ps(support_8, eval_8));
#endif
#if  ( MIC_UNROLLING_WIDTH_SP >  144 )
              res_9 = _mm512_add_ps(res_9, _mm512_mul_ps(support_9, eval_9));
#endif
#if  ( MIC_UNROLLING_WIDTH_SP >  160 )
              res_10 = _mm512_add_ps(res_10, _mm512_mul_ps(support_10, eval_10));
              res_11 = _mm512_add_ps(res_11, _mm512_mul_ps(support_11, eval_11));
#endif

#if  ( MIC_UNROLLING_WIDTH_SP >  96 )
              _mm512_store_ps(&(ptrResult[i + 96]), res_6);
              _mm512_store_ps(&(ptrResult[i + 112]), res_7);
#endif
#if  ( MIC_UNROLLING_WIDTH_SP >  128 )
              _mm512_store_ps(&(ptrResult[i + 128]), res_8);
#endif
#if  ( MIC_UNROLLING_WIDTH_SP >  144 )
              _mm512_store_ps(&(ptrResult[i + 144]), res_9);
#endif
#if  ( MIC_UNROLLING_WIDTH_SP >  160 )
              _mm512_store_ps(&(ptrResult[i + 160]), res_10);
              _mm512_store_ps(&(ptrResult[i + 176]), res_11);
#endif
            }

#else

            for (size_t ii = i; ii < i + getChunkDataPoints(); ii++) {
              for (size_t j = start_index_grid; j < end_index_grid; j++) {
                float curSupport = ptrAlpha[j];

                for (size_t d = 0; d < dims; d++) {
                  float eval = ((ptrLevel[(j * dims) + d]) * (ptrData[(d * result_size) + ii]));
                  float index_calc = eval - (ptrIndex[(j * dims) + d]);
                  float abs = static_cast<float>(fabs(index_calc));
                  float last = 1.0f - abs;
                  float localSupport = std::max<float>(last, 0.0f);
                  curSupport *= localSupport;
                }

                ptrResult[ii] += curSupport;
              }
            }

#endif
          }
        }

        static inline void multTransposeImpl(
          float* ptrLevel,
          float* ptrIndex,
          float* /*ptrMask*/, //unused for this specialization
          float* /*ptrOffset*/, //unused for this specialization
          float* ptrData,
          float* ptrSource,
          float* ptrResult,
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
              __m512 support_0 = _mm512_load_ps(&(ptrSource[i + 0]));
              __m512 support_1 = _mm512_load_ps(&(ptrSource[i + 16]));
              __m512 support_2 = _mm512_load_ps(&(ptrSource[i + 32]));
#if  ( MIC_UNROLLING_WIDTH_SP >  48 )
              __m512 support_3 = _mm512_load_ps(&(ptrSource[i + 48]));
#endif
#if  ( MIC_UNROLLING_WIDTH_SP >  64 )
              __m512 support_4 = _mm512_load_ps(&(ptrSource[i + 64]));
              __m512 support_5 = _mm512_load_ps(&(ptrSource[i + 80]));
#endif
#if  ( MIC_UNROLLING_WIDTH_SP >  96 )
              __m512 support_6 = _mm512_load_ps(&(ptrSource[i + 96]));
              __m512 support_7 = _mm512_load_ps(&(ptrSource[i + 112]));
#endif
#if  ( MIC_UNROLLING_WIDTH_SP >  128 )
              __m512 support_8 = _mm512_load_ps(&(ptrSource[i + 128]));
#endif
#if  ( MIC_UNROLLING_WIDTH_SP >  144 )
              __m512 support_9 = _mm512_load_ps(&(ptrSource[i + 144]));
#endif
#if  ( MIC_UNROLLING_WIDTH_SP >  160 )
              __m512 support_10 = _mm512_load_ps(&(ptrSource[i + 160]));
              __m512 support_11 = _mm512_load_ps(&(ptrSource[i + 176]));
#endif

              for (size_t d = 0; d < dims; d++) {
                __m512 eval_0 = _mm512_load_ps(&(ptrData[(d * sourceSize) + i + 0]));
                __m512 eval_1 = _mm512_load_ps(&(ptrData[(d * sourceSize) + i + 16]));
                __m512 eval_2 = _mm512_load_ps(&(ptrData[(d * sourceSize) + i + 32]));
#if  ( MIC_UNROLLING_WIDTH_SP >  48 )
                __m512 eval_3 = _mm512_load_ps(&(ptrData[(d * sourceSize) + i + 48]));
#endif
#if  ( MIC_UNROLLING_WIDTH_SP >  64 )
                __m512 eval_4 = _mm512_load_ps(&(ptrData[(d * sourceSize) + i + 64]));
                __m512 eval_5 = _mm512_load_ps(&(ptrData[(d * sourceSize) + i + 80]));
#endif
#if  ( MIC_UNROLLING_WIDTH_SP >  96 )
                __m512 eval_6 = _mm512_load_ps(&(ptrData[(d * sourceSize) + i + 96]));
                __m512 eval_7 = _mm512_load_ps(&(ptrData[(d * sourceSize) + i + 112]));
#endif
#if  ( MIC_UNROLLING_WIDTH_SP >  128 )
                __m512 eval_8 = _mm512_load_ps(&(ptrData[(d * sourceSize) + i + 128]));
#endif
#if  ( MIC_UNROLLING_WIDTH_SP >  144 )
                __m512 eval_9 = _mm512_load_ps(&(ptrData[(d * sourceSize) + i + 144]));
#endif
#if  ( MIC_UNROLLING_WIDTH_SP >  160 )
                __m512 eval_10 = _mm512_load_ps(&(ptrData[(d * sourceSize) + i + 160]));
                __m512 eval_11 = _mm512_load_ps(&(ptrData[(d * sourceSize) + i + 176]));
#endif

                __m512 level = _mm512_extload_ps(&(ptrLevel[(j * dims) + d]), _MM_UPCONV_PS_NONE, _MM_BROADCAST_1X16, _MM_HINT_NONE);
                __m512 index = _mm512_extload_ps(&(ptrIndex[(j * dims) + d]), _MM_UPCONV_PS_NONE, _MM_BROADCAST_1X16, _MM_HINT_NONE);

                eval_0 = _mm512_fmsub_ps(eval_0, level, index);
                eval_1 = _mm512_fmsub_ps(eval_1, level, index);
                eval_2 = _mm512_fmsub_ps(eval_2, level, index);
#if  ( MIC_UNROLLING_WIDTH_SP >  48 )
                eval_3 = _mm512_fmsub_ps(eval_3, level, index);
#endif
#if  ( MIC_UNROLLING_WIDTH_SP >  64 )
                eval_4 = _mm512_fmsub_ps(eval_4, level, index);
                eval_5 = _mm512_fmsub_ps(eval_5, level, index);
#endif
#if  ( MIC_UNROLLING_WIDTH_SP >  96 )
                eval_6 = _mm512_fmsub_ps(eval_6, level, index);
                eval_7 = _mm512_fmsub_ps(eval_7, level, index);
#endif
#if  ( MIC_UNROLLING_WIDTH_SP >  128 )
                eval_8 = _mm512_fmsub_ps(eval_8, level, index);
#endif
#if  ( MIC_UNROLLING_WIDTH_SP >  144 )
                eval_9 = _mm512_fmsub_ps(eval_9, level, index);
#endif
#if  ( MIC_UNROLLING_WIDTH_SP >  160 )
                eval_10 = _mm512_fmsub_ps(eval_10, level, index);
                eval_11 = _mm512_fmsub_ps(eval_11, level, index);
#endif

                __m512 one = _mm512_set_1to16_ps(1.0f);
                __m512 zero = _mm512_set_1to16_ps(0.0f);

                eval_0 = _mm512_gmaxabs_ps( eval_0, eval_0);
                eval_1 = _mm512_gmaxabs_ps( eval_1, eval_1);
                eval_2 = _mm512_gmaxabs_ps( eval_2, eval_2);
#if  ( MIC_UNROLLING_WIDTH_SP >  48 )
                eval_3 = _mm512_gmaxabs_ps( eval_3, eval_3);
#endif
#if  ( MIC_UNROLLING_WIDTH_SP >  64 )
                eval_4 = _mm512_gmaxabs_ps( eval_4, eval_4);
                eval_5 = _mm512_gmaxabs_ps( eval_5, eval_5);
#endif
#if  ( MIC_UNROLLING_WIDTH_SP >  96 )
                eval_6 = _mm512_gmaxabs_ps( eval_6, eval_6);
                eval_7 = _mm512_gmaxabs_ps( eval_7, eval_7);
#endif
#if  ( MIC_UNROLLING_WIDTH_SP >  128 )
                eval_8 = _mm512_gmaxabs_ps( eval_8, eval_8);
#endif
#if  ( MIC_UNROLLING_WIDTH_SP >  144 )
                eval_9 = _mm512_gmaxabs_ps( eval_9, eval_9);
#endif
#if  ( MIC_UNROLLING_WIDTH_SP >  160 )
                eval_10 = _mm512_gmaxabs_ps( eval_10, eval_10);
                eval_11 = _mm512_gmaxabs_ps( eval_11, eval_11);
#endif

                eval_0 = _mm512_sub_ps(one, eval_0);
                eval_1 = _mm512_sub_ps(one, eval_1);
                eval_2 = _mm512_sub_ps(one, eval_2);
#if  ( MIC_UNROLLING_WIDTH_SP >  48 )
                eval_3 = _mm512_sub_ps(one, eval_3);
#endif
#if  ( MIC_UNROLLING_WIDTH_SP >  64 )
                eval_4 = _mm512_sub_ps(one, eval_4);
                eval_5 = _mm512_sub_ps(one, eval_5);
#endif
#if  ( MIC_UNROLLING_WIDTH_SP >  96 )
                eval_6 = _mm512_sub_ps(one, eval_6);
                eval_7 = _mm512_sub_ps(one, eval_7);
#endif
#if  ( MIC_UNROLLING_WIDTH_SP >  128 )
                eval_8 = _mm512_sub_ps(one, eval_8);
#endif
#if  ( MIC_UNROLLING_WIDTH_SP >  144 )
                eval_9 = _mm512_sub_ps(one, eval_9);
#endif
#if  ( MIC_UNROLLING_WIDTH_SP >  160 )
                eval_10 = _mm512_sub_ps(one, eval_10);
                eval_11 = _mm512_sub_ps(one, eval_11);
#endif

                eval_0 = _mm512_gmax_ps(zero, eval_0);
                eval_1 = _mm512_gmax_ps(zero, eval_1);
                eval_2 = _mm512_gmax_ps(zero, eval_2);
#if  ( MIC_UNROLLING_WIDTH_SP >  48 )
                eval_3 = _mm512_gmax_ps(zero, eval_3);
#endif
#if  ( MIC_UNROLLING_WIDTH_SP >  64 )
                eval_4 = _mm512_gmax_ps(zero, eval_4);
                eval_5 = _mm512_gmax_ps(zero, eval_5);
#endif
#if  ( MIC_UNROLLING_WIDTH_SP >  96 )
                eval_6 = _mm512_gmax_ps(zero, eval_6);
                eval_7 = _mm512_gmax_ps(zero, eval_7);
#endif
#if  ( MIC_UNROLLING_WIDTH_SP >  128 )
                eval_8 = _mm512_gmax_ps(zero, eval_8);
#endif
#if  ( MIC_UNROLLING_WIDTH_SP >  144 )
                eval_9 = _mm512_gmax_ps(zero, eval_9);
#endif
#if  ( MIC_UNROLLING_WIDTH_SP >  160 )
                eval_10 = _mm512_gmax_ps(zero, eval_10);
                eval_11 = _mm512_gmax_ps(zero, eval_11);
#endif

                support_0 = _mm512_mul_ps(support_0, eval_0);
                support_1 = _mm512_mul_ps(support_1, eval_1);
                support_2 = _mm512_mul_ps(support_2, eval_2);
#if  ( MIC_UNROLLING_WIDTH_SP >  48 )
                support_3 = _mm512_mul_ps(support_3, eval_3);
#endif
#if  ( MIC_UNROLLING_WIDTH_SP >  64 )
                support_4 = _mm512_mul_ps(support_4, eval_4);
                support_5 = _mm512_mul_ps(support_5, eval_5);
#endif
#if  ( MIC_UNROLLING_WIDTH_SP >  96 )
                support_6 = _mm512_mul_ps(support_6, eval_6);
                support_7 = _mm512_mul_ps(support_7, eval_7);
#endif
#if  ( MIC_UNROLLING_WIDTH_SP >  128 )
                support_8 = _mm512_mul_ps(support_8, eval_8);
#endif
#if  ( MIC_UNROLLING_WIDTH_SP >  144 )
                support_9 = _mm512_mul_ps(support_9, eval_9);
#endif
#if  ( MIC_UNROLLING_WIDTH_SP >  160 )
                support_10 = _mm512_mul_ps(support_10, eval_10);
                support_11 = _mm512_mul_ps(support_11, eval_11);
#endif
              }

              _mm_prefetch((const char*) & (ptrLevel[((j + 1)*dims)]), _MM_HINT_T0);
              _mm_prefetch((const char*) & (ptrIndex[((j + 1)*dims)]), _MM_HINT_T0);
              // Prefetch exclusive -> non sync is needed
              _mm_prefetch((const char*) & (ptrResult[j + 1]), _MM_HINT_ET0);

              support_0 = _mm512_add_ps(support_0, support_1);
#if  ( MIC_UNROLLING_WIDTH_SP >  48 )
              support_2 = _mm512_add_ps(support_2, support_3);
#endif
#if  ( MIC_UNROLLING_WIDTH_SP >  64 )
              support_4 = _mm512_add_ps(support_4, support_5);
#endif
#if  ( MIC_UNROLLING_WIDTH_SP >  96 )
              support_6 = _mm512_add_ps(support_6, support_7);
#endif
#if  ( MIC_UNROLLING_WIDTH_SP >  144 )
              support_8 = _mm512_add_ps(support_8, support_9);
#endif
#if  ( MIC_UNROLLING_WIDTH_SP >  160 )
              support_10 = _mm512_add_ps(support_10, support_11);
#endif

              support_0 = _mm512_add_ps(support_0, support_2);
#if  ( MIC_UNROLLING_WIDTH_SP >  96 )
              support_4 = _mm512_add_ps(support_4, support_6);
#endif
#if  ( MIC_UNROLLING_WIDTH_SP >  160 )
              support_8 = _mm512_add_ps(support_8, support_10);
#endif

#if  ( MIC_UNROLLING_WIDTH_SP >  64 )
              support_0 = _mm512_add_ps(support_0, support_4);
#endif
#if  ( MIC_UNROLLING_WIDTH_SP >  128 )
              support_0 = _mm512_add_ps(support_0, support_8);
#endif

              ptrResult[j] += _mm512_reduce_add_ps(support_0);
            }

            _mm_prefetch((const char*) & (ptrSource[i + MIC_UNROLLING_WIDTH_SP + 0]), _MM_HINT_T0);
            _mm_prefetch((const char*) & (ptrSource[i + MIC_UNROLLING_WIDTH_SP + 16]), _MM_HINT_T0);
            _mm_prefetch((const char*) & (ptrSource[i + MIC_UNROLLING_WIDTH_SP + 32]), _MM_HINT_T0);
#if  ( MIC_UNROLLING_WIDTH_SP >  48 )
            _mm_prefetch((const char*) & (ptrSource[i + MIC_UNROLLING_WIDTH_SP + 48]), _MM_HINT_T0);
#endif
#if  ( MIC_UNROLLING_WIDTH_SP >  64 )
            _mm_prefetch((const char*) & (ptrSource[i + MIC_UNROLLING_WIDTH_SP + 64]), _MM_HINT_T0);
            _mm_prefetch((const char*) & (ptrSource[i + MIC_UNROLLING_WIDTH_SP + 80]), _MM_HINT_T0);
#endif
#if  ( MIC_UNROLLING_WIDTH_SP >  96 )
            _mm_prefetch((const char*) & (ptrSource[i + MIC_UNROLLING_WIDTH_SP + 96]), _MM_HINT_T0);
            _mm_prefetch((const char*) & (ptrSource[i + MIC_UNROLLING_WIDTH_SP + 112]), _MM_HINT_T0);
#endif
#if  ( MIC_UNROLLING_WIDTH_SP >  128 )
            _mm_prefetch((const char*) & (ptrSource[i + MIC_UNROLLING_WIDTH_SP + 128]), _MM_HINT_T0);
#endif
#if  ( MIC_UNROLLING_WIDTH_SP >  144 )
            _mm_prefetch((const char*) & (ptrSource[i + MIC_UNROLLING_WIDTH_SP + 144]), _MM_HINT_T0);
#endif
#if  ( MIC_UNROLLING_WIDTH_SP >  160 )
            _mm_prefetch((const char*) & (ptrSource[i + MIC_UNROLLING_WIDTH_SP + 160]), _MM_HINT_T0);
            _mm_prefetch((const char*) & (ptrSource[i + MIC_UNROLLING_WIDTH_SP + 176]), _MM_HINT_T0);
#endif

            size_t j = end_index_grid - 1;

            __m512 support_0 = _mm512_load_ps(&(ptrSource[i + 0]));
            __m512 support_1 = _mm512_load_ps(&(ptrSource[i + 16]));
            __m512 support_2 = _mm512_load_ps(&(ptrSource[i + 32]));
#if  ( MIC_UNROLLING_WIDTH_SP >  48 )
            __m512 support_3 = _mm512_load_ps(&(ptrSource[i + 48]));
#endif
#if  ( MIC_UNROLLING_WIDTH_SP >  64 )
            __m512 support_4 = _mm512_load_ps(&(ptrSource[i + 64]));
            __m512 support_5 = _mm512_load_ps(&(ptrSource[i + 80]));
#endif
#if  ( MIC_UNROLLING_WIDTH_SP >  96 )
            __m512 support_6 = _mm512_load_ps(&(ptrSource[i + 96]));
            __m512 support_7 = _mm512_load_ps(&(ptrSource[i + 112]));
#endif
#if  ( MIC_UNROLLING_WIDTH_SP >  128 )
            __m512 support_8 = _mm512_load_ps(&(ptrSource[i + 128]));
#endif
#if  ( MIC_UNROLLING_WIDTH_SP >  144 )
            __m512 support_9 = _mm512_load_ps(&(ptrSource[i + 144]));
#endif
#if  ( MIC_UNROLLING_WIDTH_SP >  160 )
            __m512 support_10 = _mm512_load_ps(&(ptrSource[i + 160]));
            __m512 support_11 = _mm512_load_ps(&(ptrSource[i + 176]));
#endif

            for (size_t d = 0; d < dims; d++) {
              __m512 eval_0 = _mm512_load_ps(&(ptrData[(d * sourceSize) + i + 0]));
              __m512 eval_1 = _mm512_load_ps(&(ptrData[(d * sourceSize) + i + 16]));
              __m512 eval_2 = _mm512_load_ps(&(ptrData[(d * sourceSize) + i + 32]));
#if  ( MIC_UNROLLING_WIDTH_SP >  48 )
              __m512 eval_3 = _mm512_load_ps(&(ptrData[(d * sourceSize) + i + 48]));
#endif
#if  ( MIC_UNROLLING_WIDTH_SP >  64 )
              __m512 eval_4 = _mm512_load_ps(&(ptrData[(d * sourceSize) + i + 64]));
              __m512 eval_5 = _mm512_load_ps(&(ptrData[(d * sourceSize) + i + 80]));
#endif
#if  ( MIC_UNROLLING_WIDTH_SP >  96 )
              __m512 eval_6 = _mm512_load_ps(&(ptrData[(d * sourceSize) + i + 96]));
              __m512 eval_7 = _mm512_load_ps(&(ptrData[(d * sourceSize) + i + 112]));
#endif
#if  ( MIC_UNROLLING_WIDTH_SP >  128 )
              __m512 eval_8 = _mm512_load_ps(&(ptrData[(d * sourceSize) + i + 128]));
#endif
#if  ( MIC_UNROLLING_WIDTH_SP >  144 )
              __m512 eval_9 = _mm512_load_ps(&(ptrData[(d * sourceSize) + i + 144]));
#endif
#if  ( MIC_UNROLLING_WIDTH_SP >  160 )
              __m512 eval_10 = _mm512_load_ps(&(ptrData[(d * sourceSize) + i + 160]));
              __m512 eval_11 = _mm512_load_ps(&(ptrData[(d * sourceSize) + i + 176]));
#endif

              _mm_prefetch((const char*) & (ptrData[(d * sourceSize) + i + MIC_UNROLLING_WIDTH_SP + 0]), _MM_HINT_T0);
              _mm_prefetch((const char*) & (ptrData[(d * sourceSize) + i + MIC_UNROLLING_WIDTH_SP + 16]), _MM_HINT_T0);
              _mm_prefetch((const char*) & (ptrData[(d * sourceSize) + i + MIC_UNROLLING_WIDTH_SP + 32]), _MM_HINT_T0);

              __m512 level = _mm512_extload_ps(&(ptrLevel[(j * dims) + d]), _MM_UPCONV_PS_NONE, _MM_BROADCAST_1X16, _MM_HINT_NONE);
              __m512 index = _mm512_extload_ps(&(ptrIndex[(j * dims) + d]), _MM_UPCONV_PS_NONE, _MM_BROADCAST_1X16, _MM_HINT_NONE);

              __m512 one = _mm512_set_1to16_ps(1.0f);
              __m512 zero = _mm512_set_1to16_ps(0.0f);

              eval_0 = _mm512_fmsub_ps(eval_0, level, index);
              eval_1 = _mm512_fmsub_ps(eval_1, level, index);
              eval_2 = _mm512_fmsub_ps(eval_2, level, index);
#if  ( MIC_UNROLLING_WIDTH_SP >  48 )
              eval_3 = _mm512_fmsub_ps(eval_3, level, index);
#endif
#if  ( MIC_UNROLLING_WIDTH_SP >  64 )
              eval_4 = _mm512_fmsub_ps(eval_4, level, index);
              eval_5 = _mm512_fmsub_ps(eval_5, level, index);
#endif
#if  ( MIC_UNROLLING_WIDTH_SP >  96 )
              eval_6 = _mm512_fmsub_ps(eval_6, level, index);
              eval_7 = _mm512_fmsub_ps(eval_7, level, index);
#endif
#if  ( MIC_UNROLLING_WIDTH_SP >  128 )
              eval_8 = _mm512_fmsub_ps(eval_8, level, index);
#endif
#if  ( MIC_UNROLLING_WIDTH_SP >  144 )
              eval_9 = _mm512_fmsub_ps(eval_9, level, index);
#endif
#if  ( MIC_UNROLLING_WIDTH_SP >  160 )
              eval_10 = _mm512_fmsub_ps(eval_10, level, index);
              eval_11 = _mm512_fmsub_ps(eval_11, level, index);
#endif

#if  ( MIC_UNROLLING_WIDTH_SP >  48 )
              _mm_prefetch((const char*) & (ptrData[(d * sourceSize) + i + MIC_UNROLLING_WIDTH_SP + 48]), _MM_HINT_T0);
#endif
#if  ( MIC_UNROLLING_WIDTH_SP >  64 )
              _mm_prefetch((const char*) & (ptrData[(d * sourceSize) + i + MIC_UNROLLING_WIDTH_SP + 64]), _MM_HINT_T0);
              _mm_prefetch((const char*) & (ptrData[(d * sourceSize) + i + MIC_UNROLLING_WIDTH_SP + 80]), _MM_HINT_T0);
#endif

              eval_0 = _mm512_gmaxabs_ps( eval_0, eval_0);
              eval_1 = _mm512_gmaxabs_ps( eval_1, eval_1);
              eval_2 = _mm512_gmaxabs_ps( eval_2, eval_2);
#if  ( MIC_UNROLLING_WIDTH_SP >  48 )
              eval_3 = _mm512_gmaxabs_ps( eval_3, eval_3);
#endif
#if  ( MIC_UNROLLING_WIDTH_SP >  64 )
              eval_4 = _mm512_gmaxabs_ps( eval_4, eval_4);
              eval_5 = _mm512_gmaxabs_ps( eval_5, eval_5);
#endif
#if  ( MIC_UNROLLING_WIDTH_SP >  96 )
              eval_6 = _mm512_gmaxabs_ps( eval_6, eval_6);
              eval_7 = _mm512_gmaxabs_ps( eval_7, eval_7);
#endif
#if  ( MIC_UNROLLING_WIDTH_SP >  128 )
              eval_8 = _mm512_gmaxabs_ps( eval_8, eval_8);
#endif
#if  ( MIC_UNROLLING_WIDTH_SP >  144 )
              eval_9 = _mm512_gmaxabs_ps( eval_9, eval_9);
#endif
#if  ( MIC_UNROLLING_WIDTH_SP >  160 )
              eval_10 = _mm512_gmaxabs_ps( eval_10, eval_10);
              eval_11 = _mm512_gmaxabs_ps( eval_11, eval_11);
#endif

              eval_0 = _mm512_sub_ps(one, eval_0);
              eval_1 = _mm512_sub_ps(one, eval_1);
              eval_2 = _mm512_sub_ps(one, eval_2);
#if  ( MIC_UNROLLING_WIDTH_SP >  48 )
              eval_3 = _mm512_sub_ps(one, eval_3);
#endif
#if  ( MIC_UNROLLING_WIDTH_SP >  64 )
              eval_4 = _mm512_sub_ps(one, eval_4);
              eval_5 = _mm512_sub_ps(one, eval_5);
#endif
#if  ( MIC_UNROLLING_WIDTH_SP >  96 )
              eval_6 = _mm512_sub_ps(one, eval_6);
              eval_7 = _mm512_sub_ps(one, eval_7);
#endif
#if  ( MIC_UNROLLING_WIDTH_SP >  128 )
              eval_8 = _mm512_sub_ps(one, eval_8);
#endif
#if  ( MIC_UNROLLING_WIDTH_SP >  144 )
              eval_9 = _mm512_sub_ps(one, eval_9);
#endif
#if  ( MIC_UNROLLING_WIDTH_SP >  160 )
              eval_10 = _mm512_sub_ps(one, eval_10);
              eval_11 = _mm512_sub_ps(one, eval_11);
#endif

#if  ( MIC_UNROLLING_WIDTH_SP >  96 )
              _mm_prefetch((const char*) & (ptrData[(d * sourceSize) + i + MIC_UNROLLING_WIDTH_SP + 96]), _MM_HINT_T0);
              _mm_prefetch((const char*) & (ptrData[(d * sourceSize) + i + MIC_UNROLLING_WIDTH_SP + 112]), _MM_HINT_T0);
#endif
#if  ( MIC_UNROLLING_WIDTH_SP >  128 )
              _mm_prefetch((const char*) & (ptrData[(d * sourceSize) + i + MIC_UNROLLING_WIDTH_SP + 128]), _MM_HINT_T0);
#endif

              eval_0 = _mm512_gmax_ps(zero, eval_0);
              eval_1 = _mm512_gmax_ps(zero, eval_1);
              eval_2 = _mm512_gmax_ps(zero, eval_2);
#if  ( MIC_UNROLLING_WIDTH_SP >  48 )
              eval_3 = _mm512_gmax_ps(zero, eval_3);
#endif
#if  ( MIC_UNROLLING_WIDTH_SP >  64 )
              eval_4 = _mm512_gmax_ps(zero, eval_4);
              eval_5 = _mm512_gmax_ps(zero, eval_5);
#endif
#if  ( MIC_UNROLLING_WIDTH_SP >  96 )
              eval_6 = _mm512_gmax_ps(zero, eval_6);
              eval_7 = _mm512_gmax_ps(zero, eval_7);
#endif
#if  ( MIC_UNROLLING_WIDTH_SP >  128 )
              eval_8 = _mm512_gmax_ps(zero, eval_8);
#endif
#if  ( MIC_UNROLLING_WIDTH_SP >  144 )
              eval_9 = _mm512_gmax_ps(zero, eval_9);
#endif
#if  ( MIC_UNROLLING_WIDTH_SP >  160 )
              eval_10 = _mm512_gmax_ps(zero, eval_10);
              eval_11 = _mm512_gmax_ps(zero, eval_11);
#endif

#if  ( MIC_UNROLLING_WIDTH_SP >  128 )
              _mm_prefetch((const char*) & (ptrData[(d * sourceSize) + i + MIC_UNROLLING_WIDTH_SP + 144]), _MM_HINT_T0);
#endif
#if  ( MIC_UNROLLING_WIDTH_SP >  160 )
              _mm_prefetch((const char*) & (ptrData[(d * sourceSize) + i + MIC_UNROLLING_WIDTH_SP + 160]), _MM_HINT_T0);
              _mm_prefetch((const char*) & (ptrData[(d * sourceSize) + i + MIC_UNROLLING_WIDTH_SP + 172]), _MM_HINT_T0);
#endif

              support_0 = _mm512_mul_ps(support_0, eval_0);
              support_1 = _mm512_mul_ps(support_1, eval_1);
              support_2 = _mm512_mul_ps(support_2, eval_2);
#if  ( MIC_UNROLLING_WIDTH_SP >  48 )
              support_3 = _mm512_mul_ps(support_3, eval_3);
#endif
#if  ( MIC_UNROLLING_WIDTH_SP >  64 )
              support_4 = _mm512_mul_ps(support_4, eval_4);
              support_5 = _mm512_mul_ps(support_5, eval_5);
#endif
#if  ( MIC_UNROLLING_WIDTH_SP >  96 )
              support_6 = _mm512_mul_ps(support_6, eval_6);
              support_7 = _mm512_mul_ps(support_7, eval_7);
#endif
#if  ( MIC_UNROLLING_WIDTH_SP >  128 )
              support_8 = _mm512_mul_ps(support_8, eval_8);
#endif
#if  ( MIC_UNROLLING_WIDTH_SP >  144 )
              support_9 = _mm512_mul_ps(support_9, eval_9);
#endif
#if  ( MIC_UNROLLING_WIDTH_SP >  160 )
              support_10 = _mm512_mul_ps(support_10, eval_10);
              support_11 = _mm512_mul_ps(support_11, eval_11);
#endif
            }

            support_0 = _mm512_add_ps(support_0, support_1);
#if  ( MIC_UNROLLING_WIDTH_SP >  48 )
            support_2 = _mm512_add_ps(support_2, support_3);
#endif
#if  ( MIC_UNROLLING_WIDTH_SP >  64 )
            support_4 = _mm512_add_ps(support_4, support_5);
#endif
#if  ( MIC_UNROLLING_WIDTH_SP >  96 )
            support_6 = _mm512_add_ps(support_6, support_7);
#endif
#if  ( MIC_UNROLLING_WIDTH_SP >  144 )
            support_8 = _mm512_add_ps(support_8, support_9);
#endif
#if  ( MIC_UNROLLING_WIDTH_SP >  160 )
            support_10 = _mm512_add_ps(support_10, support_11);
#endif

            support_0 = _mm512_add_ps(support_0, support_2);
#if  ( MIC_UNROLLING_WIDTH_SP >  96 )
            support_4 = _mm512_add_ps(support_4, support_6);
#endif
#if  ( MIC_UNROLLING_WIDTH_SP >  160 )
            support_8 = _mm512_add_ps(support_8, support_10);
#endif

#if  ( MIC_UNROLLING_WIDTH_SP >  64 )
            support_0 = _mm512_add_ps(support_0, support_4);
#endif
#if  ( MIC_UNROLLING_WIDTH_SP >  128 )
            support_0 = _mm512_add_ps(support_0, support_8);
#endif
            ptrResult[j] += _mm512_reduce_add_ps(support_0);
          }

#else

          for (size_t i = start_index_data; i < end_index_data; i++) {
            for (size_t j = start_index_grid; j < end_index_grid; j++) {
              float curSupport = ptrSource[i];

              for (size_t d = 0; d < dims; d++) {
                float eval = ((ptrLevel[(j * dims) + d]) * (ptrData[(d * sourceSize) + i]));
                float index_calc = eval - (ptrIndex[(j * dims) + d]);
                float abs = static_cast<float>(fabs(index_calc));
                float last = 1.0f - abs;
                float localSupport = std::max<float>(last, 0.0f);
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
#endif // SPMICLINEAR_HPP
