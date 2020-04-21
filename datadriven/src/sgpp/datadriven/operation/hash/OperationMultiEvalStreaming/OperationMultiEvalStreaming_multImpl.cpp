// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <algorithm>
#include <cmath>

#if defined(__SSE3__) && !defined(__AVX__)
#include <pmmintrin.h>
#endif

#if defined(__SSE3__) && defined(__AVX__)
#include <immintrin.h>
#endif

#if defined(__MIC__)
#include <immintrin.h>  // NOLINT(build/include)
#endif

#include <sgpp/datadriven/operation/hash/OperationMultiEvalStreaming/OperationMultiEvalStreaming.hpp>
#include <sgpp/globaldef.hpp>

namespace sgpp {
namespace datadriven {

#if defined(__SSE3__) && !defined(__AVX__) && !defined(__AVX512F__)
void OperationMultiEvalStreaming::multImpl(
    sgpp::base::DataMatrix* level, sgpp::base::DataMatrix* index, sgpp::base::DataMatrix* dataset,
    sgpp::base::DataVector& alpha, sgpp::base::DataVector& result, const size_t start_index_grid,
    const size_t end_index_grid, const size_t start_index_data, const size_t end_index_data) {
  double* ptrLevel = level->getPointer();
  double* ptrIndex = index->getPointer();
  double* ptrAlpha = alpha.getPointer();
  double* ptrData = dataset->getPointer();
  double* ptrResult = result.getPointer();
  size_t result_size = result.getSize();
  size_t dims = dataset->getNrows();

  for (size_t c = start_index_data; c < end_index_data;
       c += std::min<size_t>(getChunkDataPoints(), (end_index_data - c))) {
#ifdef __ICC
#pragma ivdep
#pragma vector aligned
#endif

    for (size_t m = start_index_grid; m < end_index_grid;
         m += std::min<size_t>(static_cast<size_t>(getChunkGridPoints()),
                               (end_index_grid - m))) {
      size_t grid_inc = std::min<size_t>(
          static_cast<size_t>(getChunkGridPoints()), (end_index_grid - m));

      uint64_t imask = 0x7FFFFFFFFFFFFFFF;
      double* fmask = reinterpret_cast<double*>(&imask);

      for (size_t i = c; i < c + getChunkDataPoints(); i += 12) {
        for (size_t j = m; j < m + grid_inc; j++) {
          __m128d support_0 = _mm_loaddup_pd(&(ptrAlpha[j]));
          __m128d support_1 = _mm_loaddup_pd(&(ptrAlpha[j]));
          __m128d support_2 = _mm_loaddup_pd(&(ptrAlpha[j]));
          __m128d support_3 = _mm_loaddup_pd(&(ptrAlpha[j]));
          __m128d support_4 = _mm_loaddup_pd(&(ptrAlpha[j]));
          __m128d support_5 = _mm_loaddup_pd(&(ptrAlpha[j]));

          __m128d mask = _mm_set1_pd(*fmask);
          __m128d one = _mm_set1_pd(1.0);
          __m128d zero = _mm_set1_pd(0.0);

          for (size_t d = 0; d < dims; d++) {
            __m128d eval_0 = _mm_load_pd(&(ptrData[(d * result_size) + i]));
            __m128d eval_1 = _mm_load_pd(&(ptrData[(d * result_size) + i + 2]));
            __m128d eval_2 = _mm_load_pd(&(ptrData[(d * result_size) + i + 4]));
            __m128d eval_3 = _mm_load_pd(&(ptrData[(d * result_size) + i + 6]));
            __m128d eval_4 = _mm_load_pd(&(ptrData[(d * result_size) + i + 8]));
            __m128d eval_5 = _mm_load_pd(&(ptrData[(d * result_size) + i + 10]));

            __m128d level = _mm_loaddup_pd(&(ptrLevel[(j * dims) + d]));
            __m128d index = _mm_loaddup_pd(&(ptrIndex[(j * dims) + d]));
#ifdef __FMA4__
            eval_0 = _mm_msub_pd(eval_0, level, index);
            eval_1 = _mm_msub_pd(eval_1, level, index);
            eval_2 = _mm_msub_pd(eval_2, level, index);
            eval_3 = _mm_msub_pd(eval_3, level, index);
            eval_4 = _mm_msub_pd(eval_4, level, index);
            eval_5 = _mm_msub_pd(eval_5, level, index);
#else
            eval_0 = _mm_sub_pd(_mm_mul_pd(eval_0, level), index);
            eval_1 = _mm_sub_pd(_mm_mul_pd(eval_1, level), index);
            eval_2 = _mm_sub_pd(_mm_mul_pd(eval_2, level), index);
            eval_3 = _mm_sub_pd(_mm_mul_pd(eval_3, level), index);
            eval_4 = _mm_sub_pd(_mm_mul_pd(eval_4, level), index);
            eval_5 = _mm_sub_pd(_mm_mul_pd(eval_5, level), index);
#endif
            eval_0 = _mm_and_pd(mask, eval_0);
            eval_1 = _mm_and_pd(mask, eval_1);
            eval_2 = _mm_and_pd(mask, eval_2);
            eval_3 = _mm_and_pd(mask, eval_3);
            eval_4 = _mm_and_pd(mask, eval_4);
            eval_5 = _mm_and_pd(mask, eval_5);

            eval_0 = _mm_sub_pd(one, eval_0);
            eval_1 = _mm_sub_pd(one, eval_1);
            eval_2 = _mm_sub_pd(one, eval_2);
            eval_3 = _mm_sub_pd(one, eval_3);
            eval_4 = _mm_sub_pd(one, eval_4);
            eval_5 = _mm_sub_pd(one, eval_5);

            eval_0 = _mm_max_pd(zero, eval_0);
            eval_1 = _mm_max_pd(zero, eval_1);
            eval_2 = _mm_max_pd(zero, eval_2);
            eval_3 = _mm_max_pd(zero, eval_3);
            eval_4 = _mm_max_pd(zero, eval_4);
            eval_5 = _mm_max_pd(zero, eval_5);

            support_0 = _mm_mul_pd(support_0, eval_0);
            support_1 = _mm_mul_pd(support_1, eval_1);
            support_2 = _mm_mul_pd(support_2, eval_2);
            support_3 = _mm_mul_pd(support_3, eval_3);
            support_4 = _mm_mul_pd(support_4, eval_4);
            support_5 = _mm_mul_pd(support_5, eval_5);
          }

          __m128d res_0 = _mm_load_pd(&(ptrResult[i]));
          __m128d res_1 = _mm_load_pd(&(ptrResult[i + 2]));
          __m128d res_2 = _mm_load_pd(&(ptrResult[i + 4]));
          __m128d res_3 = _mm_load_pd(&(ptrResult[i + 6]));
          __m128d res_4 = _mm_load_pd(&(ptrResult[i + 8]));
          __m128d res_5 = _mm_load_pd(&(ptrResult[i + 10]));

          res_0 = _mm_add_pd(res_0, support_0);
          res_1 = _mm_add_pd(res_1, support_1);
          res_2 = _mm_add_pd(res_2, support_2);
          res_3 = _mm_add_pd(res_3, support_3);
          res_4 = _mm_add_pd(res_4, support_4);
          res_5 = _mm_add_pd(res_5, support_5);

          _mm_store_pd(&(ptrResult[i]), res_0);
          _mm_store_pd(&(ptrResult[i + 2]), res_1);
          _mm_store_pd(&(ptrResult[i + 4]), res_2);
          _mm_store_pd(&(ptrResult[i + 6]), res_3);
          _mm_store_pd(&(ptrResult[i + 8]), res_4);
          _mm_store_pd(&(ptrResult[i + 10]), res_5);
        }
      }
    }
  }
}
#endif

#if defined(__SSE3__) && defined(__AVX__) && !defined(__AVX512F__)
void OperationMultiEvalStreaming::multImpl(
    sgpp::base::DataMatrix* level, sgpp::base::DataMatrix* index, sgpp::base::DataMatrix* dataset,
    sgpp::base::DataVector& alpha, sgpp::base::DataVector& result, const size_t start_index_grid,
    const size_t end_index_grid, const size_t start_index_data, const size_t end_index_data) {
  double* ptrLevel = level->getPointer();
  double* ptrIndex = index->getPointer();
  double* ptrAlpha = alpha.getPointer();
  double* ptrData = dataset->getPointer();
  double* ptrResult = result.getPointer();
  size_t result_size = result.getSize();
  size_t dims = dataset->getNrows();

  for (size_t c = start_index_data; c < end_index_data;
       c += std::min<size_t>(getChunkDataPoints(), (end_index_data - c))) {
#ifdef __ICC
#pragma ivdep
#pragma vector aligned
#endif

    for (size_t m = start_index_grid; m < end_index_grid;
         m += std::min<size_t>((size_t)getChunkGridPoints(), (end_index_grid - m))) {
      size_t grid_inc = std::min<size_t>((size_t)getChunkGridPoints(), (end_index_grid - m));

      int64_t imask = 0x7FFFFFFFFFFFFFFF;
      double* fmask = reinterpret_cast<double*>(&imask);

      for (size_t i = c; i < c + getChunkDataPoints(); i += 24) {
        for (size_t j = m; j < m + grid_inc; j++) {
          __m256d support_0 = _mm256_broadcast_sd(&(ptrAlpha[j]));
          __m256d support_1 = _mm256_broadcast_sd(&(ptrAlpha[j]));
          __m256d support_2 = _mm256_broadcast_sd(&(ptrAlpha[j]));
          __m256d support_3 = _mm256_broadcast_sd(&(ptrAlpha[j]));
          __m256d support_4 = _mm256_broadcast_sd(&(ptrAlpha[j]));
          __m256d support_5 = _mm256_broadcast_sd(&(ptrAlpha[j]));

          __m256d mask = _mm256_broadcast_sd(fmask);
          __m256d one = _mm256_set1_pd(1.0);
          __m256d zero = _mm256_set1_pd(0.0);

          for (size_t d = 0; d < dims; d++) {
            __m256d eval_0 = _mm256_load_pd(&(ptrData[(d * result_size) + i]));
            __m256d eval_1 = _mm256_load_pd(&(ptrData[(d * result_size) + i + 4]));
            __m256d eval_2 = _mm256_load_pd(&(ptrData[(d * result_size) + i + 8]));
            __m256d eval_3 = _mm256_load_pd(&(ptrData[(d * result_size) + i + 12]));
            __m256d eval_4 = _mm256_load_pd(&(ptrData[(d * result_size) + i + 16]));
            __m256d eval_5 = _mm256_load_pd(&(ptrData[(d * result_size) + i + 20]));

            __m256d level = _mm256_broadcast_sd(&(ptrLevel[(j * dims) + d]));
            __m256d index = _mm256_broadcast_sd(&(ptrIndex[(j * dims) + d]));
#ifdef __FMA4__
            eval_0 = _mm256_msub_pd(eval_0, level, index);
            eval_1 = _mm256_msub_pd(eval_1, level, index);
            eval_2 = _mm256_msub_pd(eval_2, level, index);
            eval_3 = _mm256_msub_pd(eval_3, level, index);
            eval_4 = _mm256_msub_pd(eval_4, level, index);
            eval_5 = _mm256_msub_pd(eval_5, level, index);
#else
#ifdef __AVX2__
            eval_0 = _mm256_fmsub_pd(eval_0, level, index);
            eval_1 = _mm256_fmsub_pd(eval_1, level, index);
            eval_2 = _mm256_fmsub_pd(eval_2, level, index);
            eval_3 = _mm256_fmsub_pd(eval_3, level, index);
            eval_4 = _mm256_fmsub_pd(eval_4, level, index);
            eval_5 = _mm256_fmsub_pd(eval_5, level, index);
#else
            eval_0 = _mm256_sub_pd(_mm256_mul_pd(eval_0, level), index);
            eval_1 = _mm256_sub_pd(_mm256_mul_pd(eval_1, level), index);
            eval_2 = _mm256_sub_pd(_mm256_mul_pd(eval_2, level), index);
            eval_3 = _mm256_sub_pd(_mm256_mul_pd(eval_3, level), index);
            eval_4 = _mm256_sub_pd(_mm256_mul_pd(eval_4, level), index);
            eval_5 = _mm256_sub_pd(_mm256_mul_pd(eval_5, level), index);
#endif
#endif
            eval_0 = _mm256_and_pd(mask, eval_0);
            eval_1 = _mm256_and_pd(mask, eval_1);
            eval_2 = _mm256_and_pd(mask, eval_2);
            eval_3 = _mm256_and_pd(mask, eval_3);
            eval_4 = _mm256_and_pd(mask, eval_4);
            eval_5 = _mm256_and_pd(mask, eval_5);

            eval_0 = _mm256_sub_pd(one, eval_0);
            eval_1 = _mm256_sub_pd(one, eval_1);
            eval_2 = _mm256_sub_pd(one, eval_2);
            eval_3 = _mm256_sub_pd(one, eval_3);
            eval_4 = _mm256_sub_pd(one, eval_4);
            eval_5 = _mm256_sub_pd(one, eval_5);

            eval_0 = _mm256_max_pd(zero, eval_0);
            eval_1 = _mm256_max_pd(zero, eval_1);
            eval_2 = _mm256_max_pd(zero, eval_2);
            eval_3 = _mm256_max_pd(zero, eval_3);
            eval_4 = _mm256_max_pd(zero, eval_4);
            eval_5 = _mm256_max_pd(zero, eval_5);

            support_0 = _mm256_mul_pd(support_0, eval_0);
            support_1 = _mm256_mul_pd(support_1, eval_1);
            support_2 = _mm256_mul_pd(support_2, eval_2);
            support_3 = _mm256_mul_pd(support_3, eval_3);
            support_4 = _mm256_mul_pd(support_4, eval_4);
            support_5 = _mm256_mul_pd(support_5, eval_5);
          }

          __m256d res_0 = _mm256_load_pd(&(ptrResult[i]));
          __m256d res_1 = _mm256_load_pd(&(ptrResult[i + 4]));
          __m256d res_2 = _mm256_load_pd(&(ptrResult[i + 8]));
          __m256d res_3 = _mm256_load_pd(&(ptrResult[i + 12]));
          __m256d res_4 = _mm256_load_pd(&(ptrResult[i + 16]));
          __m256d res_5 = _mm256_load_pd(&(ptrResult[i + 20]));

          res_0 = _mm256_add_pd(res_0, support_0);
          res_1 = _mm256_add_pd(res_1, support_1);
          res_2 = _mm256_add_pd(res_2, support_2);
          res_3 = _mm256_add_pd(res_3, support_3);
          res_4 = _mm256_add_pd(res_4, support_4);
          res_5 = _mm256_add_pd(res_5, support_5);

          _mm256_store_pd(&(ptrResult[i]), res_0);
          _mm256_store_pd(&(ptrResult[i + 4]), res_1);
          _mm256_store_pd(&(ptrResult[i + 8]), res_2);
          _mm256_store_pd(&(ptrResult[i + 12]), res_3);
          _mm256_store_pd(&(ptrResult[i + 16]), res_4);
          _mm256_store_pd(&(ptrResult[i + 20]), res_5);
        }
      }
    }
  }
}
#endif

#if defined(__MIC__) || defined(__AVX512F__)
void OperationMultiEvalStreaming::multImpl(
    sgpp::base::DataMatrix* level, sgpp::base::DataMatrix* index, sgpp::base::DataMatrix* dataset,
    sgpp::base::DataVector& alpha, sgpp::base::DataVector& result, const size_t start_index_grid,
    const size_t end_index_grid, const size_t start_index_data, const size_t end_index_data) {
  double* ptrLevel = level->getPointer();
  double* ptrIndex = index->getPointer();
  double* ptrAlpha = alpha.getPointer();
  double* ptrData = dataset->getPointer();
  double* ptrResult = result.getPointer();
  size_t result_size = result.getSize();
  size_t dims = dataset->getNrows();

#if defined(__MIC__)
#define _mm512_broadcast_sd(A) \
  _mm512_extload_pd(A, _MM_UPCONV_PD_NONE, _MM_BROADCAST_1X8, _MM_HINT_NONE)
#define _mm512_max_pd(A, B) _mm512_gmax_pd(A, B)
#define _mm512_set1_epi64(A) _mm512_set_1to8_epi64(A)
#define _mm512_set1_pd(A) _mm512_set_1to8_pd(A)
#endif
#if defined(__AVX512F__)
#define _mm512_broadcast_sd(A) _mm512_broadcastsd_pd(_mm_load_sd(A))
#endif

  for (size_t i = start_index_data; i < end_index_data; i += getChunkDataPoints()) {
    for (size_t j = start_index_grid; j < end_index_grid; j++) {
      _mm_prefetch((const char*)&(ptrAlpha[j + 1]), _MM_HINT_T0);
      _mm_prefetch((const char*)&(ptrLevel[((j + 1) * dims)]), _MM_HINT_T0);
      _mm_prefetch((const char*)&(ptrIndex[((j + 1) * dims)]), _MM_HINT_T0);

      __m512d support_0 = _mm512_broadcast_sd(&(ptrAlpha[j]));
      __m512d support_1 = _mm512_broadcast_sd(&(ptrAlpha[j]));
      __m512d support_2 = _mm512_broadcast_sd(&(ptrAlpha[j]));
#if (STREAMING_LINEAR_MIC_AVX512_UNROLLING_WIDTH > 24)
      __m512d support_3 = _mm512_broadcast_sd(&(ptrAlpha[j]));
#endif
#if (STREAMING_LINEAR_MIC_AVX512_UNROLLING_WIDTH > 32)
      __m512d support_4 = _mm512_broadcast_sd(&(ptrAlpha[j]));
      __m512d support_5 = _mm512_broadcast_sd(&(ptrAlpha[j]));
#endif
#if (STREAMING_LINEAR_MIC_AVX512_UNROLLING_WIDTH > 48)
      __m512d support_6 = _mm512_broadcast_sd(&(ptrAlpha[j]));
      __m512d support_7 = _mm512_broadcast_sd(&(ptrAlpha[j]));
#endif
#if (STREAMING_LINEAR_MIC_AVX512_UNROLLING_WIDTH > 64)
      __m512d support_8 = _mm512_broadcast_sd(&(ptrAlpha[j]));
#endif
#if (STREAMING_LINEAR_MIC_AVX512_UNROLLING_WIDTH > 72)
      __m512d support_9 = _mm512_broadcast_sd(&(ptrAlpha[j]));
#endif
#if (STREAMING_LINEAR_MIC_AVX512_UNROLLING_WIDTH > 80)
      __m512d support_10 = _mm512_broadcast_sd(&(ptrAlpha[j]));
      __m512d support_11 = _mm512_broadcast_sd(&(ptrAlpha[j]));
#endif

      for (size_t d = 0; d < dims; d++) {
        __m512d eval_0 = _mm512_load_pd(&(ptrData[(d * result_size) + i + 0]));
        __m512d eval_1 = _mm512_load_pd(&(ptrData[(d * result_size) + i + 8]));
        __m512d eval_2 = _mm512_load_pd(&(ptrData[(d * result_size) + i + 16]));
#if (STREAMING_LINEAR_MIC_AVX512_UNROLLING_WIDTH > 24)
        __m512d eval_3 = _mm512_load_pd(&(ptrData[(d * result_size) + i + 24]));
#endif
#if (STREAMING_LINEAR_MIC_AVX512_UNROLLING_WIDTH > 32)
        __m512d eval_4 = _mm512_load_pd(&(ptrData[(d * result_size) + i + 32]));
        __m512d eval_5 = _mm512_load_pd(&(ptrData[(d * result_size) + i + 40]));
#endif
#if (STREAMING_LINEAR_MIC_AVX512_UNROLLING_WIDTH > 48)
        __m512d eval_6 = _mm512_load_pd(&(ptrData[(d * result_size) + i + 48]));
        __m512d eval_7 = _mm512_load_pd(&(ptrData[(d * result_size) + i + 56]));
#endif
#if (STREAMING_LINEAR_MIC_AVX512_UNROLLING_WIDTH > 64)
        __m512d eval_8 = _mm512_load_pd(&(ptrData[(d * result_size) + i + 64]));
#endif
#if (STREAMING_LINEAR_MIC_AVX512_UNROLLING_WIDTH > 72)
        __m512d eval_9 = _mm512_load_pd(&(ptrData[(d * result_size) + i + 72]));
#endif
#if (STREAMING_LINEAR_MIC_AVX512_UNROLLING_WIDTH > 80)
        __m512d eval_10 = _mm512_load_pd(&(ptrData[(d * result_size) + i + 80]));
        __m512d eval_11 = _mm512_load_pd(&(ptrData[(d * result_size) + i + 88]));
#endif
        __m512d level = _mm512_broadcast_sd(&(ptrLevel[(j * dims) + d]));
        __m512d index = _mm512_broadcast_sd(&(ptrIndex[(j * dims) + d]));

        eval_0 = _mm512_fmsub_pd(eval_0, level, index);
        eval_1 = _mm512_fmsub_pd(eval_1, level, index);
        eval_2 = _mm512_fmsub_pd(eval_2, level, index);
#if (STREAMING_LINEAR_MIC_AVX512_UNROLLING_WIDTH > 24)
        eval_3 = _mm512_fmsub_pd(eval_3, level, index);
#endif
#if (STREAMING_LINEAR_MIC_AVX512_UNROLLING_WIDTH > 32)
        eval_4 = _mm512_fmsub_pd(eval_4, level, index);
        eval_5 = _mm512_fmsub_pd(eval_5, level, index);
#endif
#if (STREAMING_LINEAR_MIC_AVX512_UNROLLING_WIDTH > 48)
        eval_6 = _mm512_fmsub_pd(eval_6, level, index);
        eval_7 = _mm512_fmsub_pd(eval_7, level, index);
#endif
#if (STREAMING_LINEAR_MIC_AVX512_UNROLLING_WIDTH > 64)
        eval_8 = _mm512_fmsub_pd(eval_8, level, index);
#endif
#if (STREAMING_LINEAR_MIC_AVX512_UNROLLING_WIDTH > 72)
        eval_9 = _mm512_fmsub_pd(eval_9, level, index);
#endif
#if (STREAMING_LINEAR_MIC_AVX512_UNROLLING_WIDTH > 80)
        eval_10 = _mm512_fmsub_pd(eval_10, level, index);
        eval_11 = _mm512_fmsub_pd(eval_11, level, index);
#endif
        __m512i abs2Mask = _mm512_set1_epi64(0x7FFFFFFFFFFFFFFF);
        eval_0 = _mm512_castsi512_pd(_mm512_and_epi64(abs2Mask, _mm512_castpd_si512(eval_0)));
        eval_1 = _mm512_castsi512_pd(_mm512_and_epi64(abs2Mask, _mm512_castpd_si512(eval_1)));
        eval_2 = _mm512_castsi512_pd(_mm512_and_epi64(abs2Mask, _mm512_castpd_si512(eval_2)));
#if (STREAMING_LINEAR_MIC_AVX512_UNROLLING_WIDTH > 24)
        eval_3 = _mm512_castsi512_pd(_mm512_and_epi64(abs2Mask, _mm512_castpd_si512(eval_3)));
#endif
#if (STREAMING_LINEAR_MIC_AVX512_UNROLLING_WIDTH > 32)
        eval_4 = _mm512_castsi512_pd(_mm512_and_epi64(abs2Mask, _mm512_castpd_si512(eval_4)));
        eval_5 = _mm512_castsi512_pd(_mm512_and_epi64(abs2Mask, _mm512_castpd_si512(eval_5)));
#endif
#if (STREAMING_LINEAR_MIC_AVX512_UNROLLING_WIDTH > 48)
        eval_6 = _mm512_castsi512_pd(_mm512_and_epi64(abs2Mask, _mm512_castpd_si512(eval_6)));
        eval_7 = _mm512_castsi512_pd(_mm512_and_epi64(abs2Mask, _mm512_castpd_si512(eval_7)));
#endif
#if (STREAMING_LINEAR_MIC_AVX512_UNROLLING_WIDTH > 64)
        eval_8 = _mm512_castsi512_pd(_mm512_and_epi64(abs2Mask, _mm512_castpd_si512(eval_8)));
#endif
#if (STREAMING_LINEAR_MIC_AVX512_UNROLLING_WIDTH > 72)
        eval_9 = _mm512_castsi512_pd(_mm512_and_epi64(abs2Mask, _mm512_castpd_si512(eval_9)));
#endif
#if (STREAMING_LINEAR_MIC_AVX512_UNROLLING_WIDTH > 80)
        eval_10 = _mm512_castsi512_pd(_mm512_and_epi64(abs2Mask, _mm512_castpd_si512(eval_10)));
        eval_11 = _mm512_castsi512_pd(_mm512_and_epi64(abs2Mask, _mm512_castpd_si512(eval_11)));
#endif

        __m512d one = _mm512_set1_pd(1.0);

        eval_0 = _mm512_sub_pd(one, eval_0);
        eval_1 = _mm512_sub_pd(one, eval_1);
        eval_2 = _mm512_sub_pd(one, eval_2);
#if (STREAMING_LINEAR_MIC_AVX512_UNROLLING_WIDTH > 24)
        eval_3 = _mm512_sub_pd(one, eval_3);
#endif
#if (STREAMING_LINEAR_MIC_AVX512_UNROLLING_WIDTH > 32)
        eval_4 = _mm512_sub_pd(one, eval_4);
        eval_5 = _mm512_sub_pd(one, eval_5);
#endif
#if (STREAMING_LINEAR_MIC_AVX512_UNROLLING_WIDTH > 48)
        eval_6 = _mm512_sub_pd(one, eval_6);
        eval_7 = _mm512_sub_pd(one, eval_7);
#endif
#if (STREAMING_LINEAR_MIC_AVX512_UNROLLING_WIDTH > 64)
        eval_8 = _mm512_sub_pd(one, eval_8);
#endif
#if (STREAMING_LINEAR_MIC_AVX512_UNROLLING_WIDTH > 72)
        eval_9 = _mm512_sub_pd(one, eval_9);
#endif
#if (STREAMING_LINEAR_MIC_AVX512_UNROLLING_WIDTH > 80)
        eval_10 = _mm512_sub_pd(one, eval_10);
        eval_11 = _mm512_sub_pd(one, eval_11);
#endif
        __m512d zero = _mm512_setzero_pd();

        eval_0 = _mm512_max_pd(zero, eval_0);
        eval_1 = _mm512_max_pd(zero, eval_1);
        eval_2 = _mm512_max_pd(zero, eval_2);
#if (STREAMING_LINEAR_MIC_AVX512_UNROLLING_WIDTH > 24)
        eval_3 = _mm512_max_pd(zero, eval_3);
#endif
#if (STREAMING_LINEAR_MIC_AVX512_UNROLLING_WIDTH > 32)
        eval_4 = _mm512_max_pd(zero, eval_4);
        eval_5 = _mm512_max_pd(zero, eval_5);
#endif
#if (STREAMING_LINEAR_MIC_AVX512_UNROLLING_WIDTH > 48)
        eval_6 = _mm512_max_pd(zero, eval_6);
        eval_7 = _mm512_max_pd(zero, eval_7);
#endif
#if (STREAMING_LINEAR_MIC_AVX512_UNROLLING_WIDTH > 64)
        eval_8 = _mm512_max_pd(zero, eval_8);
#endif
#if (STREAMING_LINEAR_MIC_AVX512_UNROLLING_WIDTH > 72)
        eval_9 = _mm512_max_pd(zero, eval_9);
#endif
#if (STREAMING_LINEAR_MIC_AVX512_UNROLLING_WIDTH > 80)
        eval_10 = _mm512_max_pd(zero, eval_10);
        eval_11 = _mm512_max_pd(zero, eval_11);
#endif

        support_0 = _mm512_mul_pd(support_0, eval_0);
        support_1 = _mm512_mul_pd(support_1, eval_1);
        support_2 = _mm512_mul_pd(support_2, eval_2);
#if (STREAMING_LINEAR_MIC_AVX512_UNROLLING_WIDTH > 24)
        support_3 = _mm512_mul_pd(support_3, eval_3);
#endif
#if (STREAMING_LINEAR_MIC_AVX512_UNROLLING_WIDTH > 32)
        support_4 = _mm512_mul_pd(support_4, eval_4);
        support_5 = _mm512_mul_pd(support_5, eval_5);
#endif
#if (STREAMING_LINEAR_MIC_AVX512_UNROLLING_WIDTH > 48)
        support_6 = _mm512_mul_pd(support_6, eval_6);
        support_7 = _mm512_mul_pd(support_7, eval_7);
#endif
#if (STREAMING_LINEAR_MIC_AVX512_UNROLLING_WIDTH > 64)
        support_8 = _mm512_mul_pd(support_8, eval_8);
#endif
#if (STREAMING_LINEAR_MIC_AVX512_UNROLLING_WIDTH > 72)
        support_9 = _mm512_mul_pd(support_9, eval_9);
#endif
#if (STREAMING_LINEAR_MIC_AVX512_UNROLLING_WIDTH > 80)
        support_10 = _mm512_mul_pd(support_10, eval_10);
        support_11 = _mm512_mul_pd(support_11, eval_11);
#endif
      }

      __m512d res_0 = _mm512_load_pd(&(ptrResult[i + 0]));
      __m512d res_1 = _mm512_load_pd(&(ptrResult[i + 8]));
      __m512d res_2 = _mm512_load_pd(&(ptrResult[i + 16]));
#if (STREAMING_LINEAR_MIC_AVX512_UNROLLING_WIDTH > 24)
      __m512d res_3 = _mm512_load_pd(&(ptrResult[i + 24]));
#endif
#if (STREAMING_LINEAR_MIC_AVX512_UNROLLING_WIDTH > 32)
      __m512d res_4 = _mm512_load_pd(&(ptrResult[i + 32]));
      __m512d res_5 = _mm512_load_pd(&(ptrResult[i + 40]));
#endif

      res_0 = _mm512_add_pd(res_0, support_0);
      res_1 = _mm512_add_pd(res_1, support_1);
      res_2 = _mm512_add_pd(res_2, support_2);
#if (STREAMING_LINEAR_MIC_AVX512_UNROLLING_WIDTH > 24)
      res_3 = _mm512_add_pd(res_3, support_3);
#endif
#if (STREAMING_LINEAR_MIC_AVX512_UNROLLING_WIDTH > 32)
      res_4 = _mm512_add_pd(res_4, support_4);
      res_5 = _mm512_add_pd(res_5, support_5);
#endif

      _mm512_store_pd(&(ptrResult[i + 0]), res_0);
      _mm512_store_pd(&(ptrResult[i + 8]), res_1);
      _mm512_store_pd(&(ptrResult[i + 16]), res_2);
#if (STREAMING_LINEAR_MIC_AVX512_UNROLLING_WIDTH > 24)
      _mm512_store_pd(&(ptrResult[i + 24]), res_3);
#endif
#if (STREAMING_LINEAR_MIC_AVX512_UNROLLING_WIDTH > 32)
      _mm512_store_pd(&(ptrResult[i + 32]), res_4);
      _mm512_store_pd(&(ptrResult[i + 40]), res_5);
#endif

#if (STREAMING_LINEAR_MIC_AVX512_UNROLLING_WIDTH > 48)
      __m512d res_6 = _mm512_load_pd(&(ptrResult[i + 48]));
      __m512d res_7 = _mm512_load_pd(&(ptrResult[i + 56]));
#endif
#if (STREAMING_LINEAR_MIC_AVX512_UNROLLING_WIDTH > 64)
      __m512d res_8 = _mm512_load_pd(&(ptrResult[i + 64]));
#endif
#if (STREAMING_LINEAR_MIC_AVX512_UNROLLING_WIDTH > 72)
      __m512d res_9 = _mm512_load_pd(&(ptrResult[i + 72]));
#endif
#if (STREAMING_LINEAR_MIC_AVX512_UNROLLING_WIDTH > 80)
      __m512d res_10 = _mm512_load_pd(&(ptrResult[i + 80]));
      __m512d res_11 = _mm512_load_pd(&(ptrResult[i + 88]));
#endif

#if (STREAMING_LINEAR_MIC_AVX512_UNROLLING_WIDTH > 48)
      res_6 = _mm512_add_pd(res_6, support_6);
      res_7 = _mm512_add_pd(res_7, support_7);
#endif
#if (STREAMING_LINEAR_MIC_AVX512_UNROLLING_WIDTH > 64)
      res_8 = _mm512_add_pd(res_8, support_8);
#endif
#if (STREAMING_LINEAR_MIC_AVX512_UNROLLING_WIDTH > 72)
      res_9 = _mm512_add_pd(res_9, support_9);
#endif
#if (STREAMING_LINEAR_MIC_AVX512_UNROLLING_WIDTH > 80)
      res_10 = _mm512_add_pd(res_10, support_10);
      res_11 = _mm512_add_pd(res_11, support_11);
#endif

#if (STREAMING_LINEAR_MIC_AVX512_UNROLLING_WIDTH > 48)
      _mm512_store_pd(&(ptrResult[i + 48]), res_6);
      _mm512_store_pd(&(ptrResult[i + 56]), res_7);
#endif
#if (STREAMING_LINEAR_MIC_AVX512_UNROLLING_WIDTH > 64)
      _mm512_store_pd(&(ptrResult[i + 64]), res_8);
#endif
#if (STREAMING_LINEAR_MIC_AVX512_UNROLLING_WIDTH > 72)
      _mm512_store_pd(&(ptrResult[i + 72]), res_9);
#endif
#if (STREAMING_LINEAR_MIC_AVX512_UNROLLING_WIDTH > 80)
      _mm512_store_pd(&(ptrResult[i + 80]), res_10);
      _mm512_store_pd(&(ptrResult[i + 88]), res_11);
#endif
    }
  }
  //}
  //}
}
#endif

#if !defined(__SSE3__) && !defined(__AVX__) && !defined(__MIC__) && !defined(__AVX512F__)
void OperationMultiEvalStreaming::multImpl(
    sgpp::base::DataMatrix* level, sgpp::base::DataMatrix* index, sgpp::base::DataMatrix* dataset,
    sgpp::base::DataVector& alpha, sgpp::base::DataVector& result, const size_t start_index_grid,
    const size_t end_index_grid, const size_t start_index_data, const size_t end_index_data) {
  double* ptrLevel = level->getPointer();
  double* ptrIndex = index->getPointer();
  double* ptrAlpha = alpha.getPointer();
  double* ptrData = dataset->getPointer();
  double* ptrResult = result.getPointer();
  size_t result_size = result.getSize();
  size_t dims = dataset->getNrows();

#warning "warning: using fallback implementation for OperationMultiEvalStreaming mult kernel"

  for (size_t c = start_index_data; c < end_index_data;
       c += std::min<size_t>(getChunkDataPoints(), (end_index_data - c))) {
    size_t data_end = std::min<size_t>((size_t)getChunkDataPoints() + c, end_index_data);

#ifdef __ICC
#pragma ivdep
#pragma vector aligned
#endif

    for (size_t m = start_index_grid; m < end_index_grid;
         m += std::min<size_t>((size_t)getChunkGridPoints(), (end_index_grid - m))) {
      size_t grid_end = std::min<size_t>((size_t)getChunkGridPoints() + m, end_index_grid);

      for (size_t i = c; i < data_end; i++) {
        for (size_t j = m; j < grid_end; j++) {
          double curSupport = ptrAlpha[j];

          for (size_t d = 0; d < dims; d++) {
            double eval = ((ptrLevel[(j * dims) + d]) * (ptrData[(d * result_size) + i]));
            double index_calc = eval - (ptrIndex[(j * dims) + d]);
            double abs = std::fabs(index_calc);
            double last = 1.0 - abs;
            double localSupport = std::max<double>(last, 0.0);
            curSupport *= localSupport;
          }

          ptrResult[i] += curSupport;
        }
      }
    }
  }
}
#endif

}  // namespace datadriven
}  // namespace sgpp
