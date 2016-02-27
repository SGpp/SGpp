// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#ifndef X86SIMDLINEAR_HPP
#define X86SIMDLINEAR_HPP

#include <sgpp/parallel/datadriven/basis/common/X86SimdKernelBase.hpp>

#include <sgpp/globaldef.hpp>
#include <algorithm>
namespace sgpp {
namespace parallel {

class X86SimdLinear : public X86SimdKernelBase {
 public:
  static const KernelType kernelType = Standard;
  static inline void multImpl(sgpp::base::DataMatrix* level, sgpp::base::DataMatrix* index,
                              sgpp::base::DataMatrix* /*mask*/,    // unused for this specialization
                              sgpp::base::DataMatrix* /*offset*/,  // unused for this specialization
                              sgpp::base::DataMatrix* dataset, sgpp::base::DataVector& alpha,
                              sgpp::base::DataVector& result, const size_t start_index_grid,
                              const size_t end_index_grid, const size_t start_index_data,
                              const size_t end_index_data) {
    double* ptrLevel = level->getPointer();
    double* ptrIndex = index->getPointer();
    double* ptrAlpha = alpha.getPointer();
    double* ptrData = dataset->getPointer();
    double* ptrResult = result.getPointer();
    size_t result_size = result.getSize();
    size_t dims = dataset->getNrows();

    CHECK_ARGS_MULT(level, dataset, result, start_index_grid, end_index_grid, start_index_data,
                    end_index_data);

    for (size_t c = start_index_data; c < end_index_data;
         c += std::min<size_t>(getChunkDataPoints(), (end_index_data - c))) {
      size_t data_end = std::min<size_t>((size_t)getChunkDataPoints() + c, end_index_data);

#ifdef __ICC
#pragma ivdep
#pragma vector aligned
#endif

      for (size_t i = c; i < data_end; i++) {
        ptrResult[i] = 0.0;
      }

      for (size_t m = start_index_grid; m < end_index_grid;
           m += std::min<size_t>((size_t)getChunkGridPoints(), (end_index_grid - m))) {
#if defined(__SSE3__) && !defined(__AVX__)
        size_t grid_inc = std::min<size_t>((size_t)getChunkGridPoints(), (end_index_grid - m));

        int64_t imask = 0x7FFFFFFFFFFFFFFF;
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

#endif
#if defined(__SSE3__) && defined(__AVX__)
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
              eval_0 = _mm256_sub_pd(_mm256_mul_pd(eval_0, level), index);
              eval_1 = _mm256_sub_pd(_mm256_mul_pd(eval_1, level), index);
              eval_2 = _mm256_sub_pd(_mm256_mul_pd(eval_2, level), index);
              eval_3 = _mm256_sub_pd(_mm256_mul_pd(eval_3, level), index);
              eval_4 = _mm256_sub_pd(_mm256_mul_pd(eval_4, level), index);
              eval_5 = _mm256_sub_pd(_mm256_mul_pd(eval_5, level), index);
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

#endif
#if !defined(__SSE3__) && !defined(__AVX__)
        size_t grid_end = std::min<size_t>((size_t)getChunkGridPoints() + m, end_index_grid);

        for (size_t i = c; i < data_end; i++) {
          for (size_t j = m; j < grid_end; j++) {
            double curSupport = ptrAlpha[j];

            for (size_t d = 0; d < dims; d++) {
              double eval = ((ptrLevel[(j * dims) + d]) * (ptrData[(d * result_size) + i]));
              double index_calc = eval - (ptrIndex[(j * dims) + d]);
              double abs = fabs(index_calc);
              double last = 1.0 - abs;
              double localSupport = std::max<double>(last, 0.0);
              curSupport *= localSupport;
            }

            ptrResult[i] += curSupport;
          }
        }

#endif
      }
    }
  }

  static inline void multTransposeImpl(
      sgpp::base::DataMatrix* level, sgpp::base::DataMatrix* index,
      sgpp::base::DataMatrix* /*mask*/,    // unused for this specialization
      sgpp::base::DataMatrix* /*offset*/,  // unused for this specialization
      sgpp::base::DataMatrix* dataset, sgpp::base::DataVector& source,
      sgpp::base::DataVector& result, const size_t start_index_grid, const size_t end_index_grid,
      const size_t start_index_data, const size_t end_index_data) {
    double* ptrLevel = level->getPointer();
    double* ptrIndex = index->getPointer();
    double* ptrSource = source.getPointer();
    double* ptrData = dataset->getPointer();
    double* ptrResult = result.getPointer();
    size_t source_size = source.getSize();
    size_t dims = dataset->getNrows();

    CHECK_ARGS_MULTTRANSPOSE(level, dataset, source, start_index_grid, end_index_grid,
                             start_index_data, end_index_data);

    for (size_t k = start_index_grid; k < end_index_grid;
         k += std::min<size_t>(getChunkGridPoints(), (end_index_grid - k))) {
      size_t grid_inc = std::min<size_t>((size_t)getChunkGridPoints(), (end_index_grid - k));

#if defined(__SSE3__) && !defined(__AVX__)
      int64_t imask = 0x7FFFFFFFFFFFFFFF;
      double* fmask = reinterpret_cast<double*>(&imask);

      for (size_t i = start_index_data; i < end_index_data; i += 12) {
        for (size_t j = k; j < k + grid_inc; j++) {
          __m128d support_0 = _mm_load_pd(&(ptrSource[i]));
          __m128d support_1 = _mm_load_pd(&(ptrSource[i + 2]));
          __m128d support_2 = _mm_load_pd(&(ptrSource[i + 4]));
          __m128d support_3 = _mm_load_pd(&(ptrSource[i + 6]));
          __m128d support_4 = _mm_load_pd(&(ptrSource[i + 8]));
          __m128d support_5 = _mm_load_pd(&(ptrSource[i + 10]));

          __m128d mask = _mm_set1_pd(*fmask);
          __m128d one = _mm_set1_pd(1.0);
          __m128d zero = _mm_set1_pd(0.0);

          for (size_t d = 0; d < dims; d++) {
            __m128d eval_0 = _mm_load_pd(&(ptrData[(d * source_size) + i]));
            __m128d eval_1 = _mm_load_pd(&(ptrData[(d * source_size) + i + 2]));
            __m128d eval_2 = _mm_load_pd(&(ptrData[(d * source_size) + i + 4]));
            __m128d eval_3 = _mm_load_pd(&(ptrData[(d * source_size) + i + 6]));
            __m128d eval_4 = _mm_load_pd(&(ptrData[(d * source_size) + i + 8]));
            __m128d eval_5 = _mm_load_pd(&(ptrData[(d * source_size) + i + 10]));

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

          __m128d res_0 = _mm_setzero_pd();
          res_0 = _mm_loadl_pd(res_0, &(ptrResult[j]));

          support_0 = _mm_add_pd(support_0, support_1);
          support_2 = _mm_add_pd(support_2, support_3);
          support_4 = _mm_add_pd(support_4, support_5);
          support_0 = _mm_add_pd(support_0, support_2);
          support_0 = _mm_add_pd(support_0, support_4);

          support_0 = _mm_hadd_pd(support_0, support_0);
          res_0 = _mm_add_sd(res_0, support_0);

          _mm_storel_pd(&(ptrResult[j]), res_0);
        }
      }

#endif
#if defined(__SSE3__) && defined(__AVX__)
      int64_t imask = 0x7FFFFFFFFFFFFFFF;
      double* fmask = reinterpret_cast<double*>(&imask);

      for (size_t i = start_index_data; i < end_index_data; i += 24) {
        for (size_t j = k; j < k + grid_inc; j++) {
          __m256d support_0 = _mm256_load_pd(&(ptrSource[i]));
          __m256d support_1 = _mm256_load_pd(&(ptrSource[i + 4]));
          __m256d support_2 = _mm256_load_pd(&(ptrSource[i + 8]));
          __m256d support_3 = _mm256_load_pd(&(ptrSource[i + 12]));
          __m256d support_4 = _mm256_load_pd(&(ptrSource[i + 16]));
          __m256d support_5 = _mm256_load_pd(&(ptrSource[i + 20]));

          __m256d mask = _mm256_broadcast_sd(fmask);
          __m256d one = _mm256_set1_pd(1.0);
          __m256d zero = _mm256_set1_pd(0.0);

          for (size_t d = 0; d < dims; d++) {
            __m256d eval_0 = _mm256_load_pd(&(ptrData[(d * source_size) + i]));
            __m256d eval_1 = _mm256_load_pd(&(ptrData[(d * source_size) + i + 4]));
            __m256d eval_2 = _mm256_load_pd(&(ptrData[(d * source_size) + i + 8]));
            __m256d eval_3 = _mm256_load_pd(&(ptrData[(d * source_size) + i + 12]));
            __m256d eval_4 = _mm256_load_pd(&(ptrData[(d * source_size) + i + 16]));
            __m256d eval_5 = _mm256_load_pd(&(ptrData[(d * source_size) + i + 20]));

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
            eval_0 = _mm256_sub_pd(_mm256_mul_pd(eval_0, level), index);
            eval_1 = _mm256_sub_pd(_mm256_mul_pd(eval_1, level), index);
            eval_2 = _mm256_sub_pd(_mm256_mul_pd(eval_2, level), index);
            eval_3 = _mm256_sub_pd(_mm256_mul_pd(eval_3, level), index);
            eval_4 = _mm256_sub_pd(_mm256_mul_pd(eval_4, level), index);
            eval_5 = _mm256_sub_pd(_mm256_mul_pd(eval_5, level), index);
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

          const __m256i ldStMaskAVX = _mm256_set_epi64x(0x0000000000000000, 0x0000000000000000,
                                                        0x0000000000000000, 0xFFFFFFFFFFFFFFFF);

          support_0 = _mm256_add_pd(support_0, support_1);
          support_2 = _mm256_add_pd(support_2, support_3);
          support_4 = _mm256_add_pd(support_4, support_5);
          support_0 = _mm256_add_pd(support_0, support_2);
          support_0 = _mm256_add_pd(support_0, support_4);

          support_0 = _mm256_hadd_pd(support_0, support_0);
          __m256d tmp = _mm256_permute2f128_pd(support_0, support_0, 0x81);
          support_0 = _mm256_add_pd(support_0, tmp);

// Workaround: bug with maskload in GCC (4.6.1)
#ifdef __ICC
          __m256d res_0 = _mm256_maskload_pd(&(ptrResult[j]), ldStMaskAVX);
          res_0 = _mm256_add_pd(res_0, support_0);
          _mm256_maskstore_pd(&(ptrResult[j]), ldStMaskAVX, res_0);
#else
          double tmp_reduce;
          _mm256_maskstore_pd(&(tmp_reduce), ldStMaskAVX, support_0);
          ptrResult[j] += tmp_reduce;
#endif
        }
      }

#endif
#if !defined(__SSE3__) && !defined(__AVX__)

      for (size_t i = start_index_data; i < end_index_data; i++) {
        for (size_t j = k; j < k + grid_inc; j++) {
          double curSupport = ptrSource[i];

          for (size_t d = 0; d < dims; d++) {
            double eval = ((ptrLevel[(j * dims) + d]) * (ptrData[(d * source_size) + i]));
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
  }
};

class X86SimdLinear1 : public X86SimdKernelBase1 {
 public:
  static const KernelType kernelType = Standard;
  static inline void multImpl(sgpp::base::DataMatrix* level, sgpp::base::DataMatrix* index,
                              sgpp::base::DataMatrix* /*mask*/,    // unused for this specialization
                              sgpp::base::DataMatrix* /*offset*/,  // unused for this specialization
                              sgpp::base::DataMatrix* dataset, sgpp::base::DataVector& alpha,
                              sgpp::base::DataVector& result, const size_t start_index_grid,
                              const size_t end_index_grid, const size_t start_index_data,
                              const size_t end_index_data) {
    double* ptrLevel = level->getPointer();
    double* ptrIndex = index->getPointer();
    double* ptrAlpha = alpha.getPointer();
    double* ptrData = dataset->getPointer();
    double* ptrResult = result.getPointer();
    size_t result_size = result.getSize();
    size_t dims = dataset->getNrows();

    CHECK_ARGS_MULT(level, dataset, result, start_index_grid, end_index_grid, start_index_data,
                    end_index_data);

    for (size_t c = start_index_data; c < end_index_data;
         c += std::min<size_t>(getChunkDataPoints(), (end_index_data - c))) {
      size_t data_end = std::min<size_t>((size_t)getChunkDataPoints() + c, end_index_data);

      for (size_t i = c; i < data_end; i++) {
        ptrResult[i] = 0.0;
      }

      for (size_t m = start_index_grid; m < end_index_grid;
           m += std::min<size_t>((size_t)getChunkGridPoints(), (end_index_grid - m))) {
        size_t grid_end = std::min<size_t>((size_t)getChunkGridPoints() + m, end_index_grid);

        for (size_t i = c; i < data_end; i++) {
          for (size_t j = m; j < grid_end; j++) {
            double curSupport = ptrAlpha[j];

            for (size_t d = 0; d < dims; d++) {
              double eval = ((ptrLevel[(j * dims) + d]) * (ptrData[(d * result_size) + i]));
              double index_calc = eval - (ptrIndex[(j * dims) + d]);
              double abs = fabs(index_calc);
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

  static inline void multTransposeImpl(
      sgpp::base::DataMatrix* level, sgpp::base::DataMatrix* index,
      sgpp::base::DataMatrix* /*mask*/,    // unused for this specialization
      sgpp::base::DataMatrix* /*offset*/,  // unused for this specialization
      sgpp::base::DataMatrix* dataset, sgpp::base::DataVector& source,
      sgpp::base::DataVector& result, const size_t start_index_grid, const size_t end_index_grid,
      const size_t start_index_data, const size_t end_index_data) {
    double* ptrLevel = level->getPointer();
    double* ptrIndex = index->getPointer();
    double* ptrSource = source.getPointer();
    double* ptrData = dataset->getPointer();
    double* ptrResult = result.getPointer();
    size_t source_size = source.getSize();
    size_t dims = dataset->getNrows();

    CHECK_ARGS_MULTTRANSPOSE(level, dataset, source, start_index_grid, end_index_grid,
                             start_index_data, end_index_data);

    for (size_t k = start_index_grid; k < end_index_grid;
         k += std::min<size_t>(getChunkGridPoints(), (end_index_grid - k))) {
      size_t grid_inc = std::min<size_t>((size_t)getChunkGridPoints(), (end_index_grid - k));

      for (size_t i = start_index_data; i < end_index_data; i++) {
        for (size_t j = k; j < k + grid_inc; j++) {
          double curSupport = ptrSource[i];

          for (size_t d = 0; d < dims; d++) {
            double eval = ((ptrLevel[(j * dims) + d]) * (ptrData[(d * source_size) + i]));
            double index_calc = eval - (ptrIndex[(j * dims) + d]);
            double abs = fabs(index_calc);
            double last = 1.0 - abs;
            double localSupport = std::max<double>(last, 0.0);
            curSupport *= localSupport;
          }

          ptrResult[j] += curSupport;
        }
      }
    }
  }
};

}  // namespace parallel
}  // namespace sgpp
#endif  // X86SIMDLINEAR_HPP
