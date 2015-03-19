// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#ifndef SPX86SIMDMODLINEAR_HPP
#define SPX86SIMDMODLINEAR_HPP

#include <sgpp/parallel/datadriven/basis/common/SPX86SimdKernelBase.hpp>

#include <sgpp/globaldef.hpp>


namespace SGPP {
  namespace parallel {

    class SPX86SimdModLinear : public SPX86SimdKernelBase {
      public:
        static const KernelType kernelType = Standard;
        static inline void multImpl(
          SGPP::base::DataMatrixSP* level,
          SGPP::base::DataMatrixSP* index,
          SGPP::base::DataMatrixSP* /*mask*/, //unused for this specialization
          SGPP::base::DataMatrixSP* /*offset*/, //unused for this specialization
          SGPP::base::DataMatrixSP* dataset,
          SGPP::base::DataVectorSP& alpha,
          SGPP::base::DataVectorSP& result,
          const size_t start_index_grid,
          const size_t end_index_grid,
          const size_t start_index_data,
          const size_t end_index_data) {
          float* ptrLevel = level->getPointer();
          float* ptrIndex = index->getPointer();
          float* ptrAlpha = alpha.getPointer();
          float* ptrData = dataset->getPointer();
          float* ptrResult = result.getPointer();
          size_t result_size = result.getSize();
          size_t dims = dataset->getNrows();

          CHECK_ARGS_MULT(level, dataset, result, start_index_grid, end_index_grid, start_index_data, end_index_data);

          for (size_t c = start_index_data; c < end_index_data; c += std::min<size_t>((size_t)getChunkDataPoints(), (end_index_data - c))) {
            size_t data_end = std::min<size_t>((size_t)getChunkDataPoints() + c, end_index_data);

#ifdef __ICC
#pragma ivdep
#pragma vector aligned
#endif

            for (size_t i = c; i < data_end; i++) {
              ptrResult[i] = 0.0f;
            }

            for (size_t m = start_index_grid; m < end_index_grid; m += std::min<size_t>((size_t)getChunkGridPoints(), (end_index_grid - m))) {
#if defined(__SSE3__) && !defined(__AVX__)
              size_t grid_inc = std::min<size_t>((size_t)getChunkGridPoints(), (end_index_grid - m));

              int imask = 0x7FFFFFFF;
              float* fmask = (float*)&imask;

              for (size_t i = c; i < c + getChunkDataPoints(); i += 24) {
                for (size_t j = m; j < m + grid_inc; j++) {
                  __m128 support_0 = _mm_load1_ps(&(ptrAlpha[j]));
                  __m128 support_1 = _mm_load1_ps(&(ptrAlpha[j]));
                  __m128 support_2 = _mm_load1_ps(&(ptrAlpha[j]));
                  __m128 support_3 = _mm_load1_ps(&(ptrAlpha[j]));
                  __m128 support_4 = _mm_load1_ps(&(ptrAlpha[j]));
                  __m128 support_5 = _mm_load1_ps(&(ptrAlpha[j]));

                  __m128 one = _mm_set1_ps(1.0f);
                  __m128 two = _mm_set1_ps(2.0f);
                  __m128 zero = _mm_set1_ps(0.0f);

                  for (size_t d = 0; d < dims; d++) {
                    // special case for level 1
                    if (ptrLevel[(j * dims) + d] == 2.0f) {
                      // Nothing (multiply by one)
                    }
                    // most left basis function on every level
                    else if (ptrIndex[(j * dims) + d] == 1.0f) {
                      __m128 eval_0 = _mm_load_ps(&(ptrData[(d * result_size) + i]));
                      __m128 eval_1 = _mm_load_ps(&(ptrData[(d * result_size) + i + 4]));
                      __m128 eval_2 = _mm_load_ps(&(ptrData[(d * result_size) + i + 8]));
                      __m128 eval_3 = _mm_load_ps(&(ptrData[(d * result_size) + i + 12]));
                      __m128 eval_4 = _mm_load_ps(&(ptrData[(d * result_size) + i + 16]));
                      __m128 eval_5 = _mm_load_ps(&(ptrData[(d * result_size) + i + 20]));

                      __m128 level = _mm_load1_ps(&(ptrLevel[(j * dims) + d]));

                      eval_0 = _mm_mul_ps(eval_0, level);
                      eval_1 = _mm_mul_ps(eval_1, level);
                      eval_2 = _mm_mul_ps(eval_2, level);
                      eval_3 = _mm_mul_ps(eval_3, level);
                      eval_4 = _mm_mul_ps(eval_4, level);
                      eval_5 = _mm_mul_ps(eval_5, level);

                      eval_0 = _mm_sub_ps(two, eval_0);
                      eval_1 = _mm_sub_ps(two, eval_1);
                      eval_2 = _mm_sub_ps(two, eval_2);
                      eval_3 = _mm_sub_ps(two, eval_3);
                      eval_4 = _mm_sub_ps(two, eval_4);
                      eval_5 = _mm_sub_ps(two, eval_5);

                      eval_0 = _mm_max_ps(zero, eval_0);
                      eval_1 = _mm_max_ps(zero, eval_1);
                      eval_2 = _mm_max_ps(zero, eval_2);
                      eval_3 = _mm_max_ps(zero, eval_3);
                      eval_4 = _mm_max_ps(zero, eval_4);
                      eval_5 = _mm_max_ps(zero, eval_5);

                      support_0 = _mm_mul_ps(support_0, eval_0);
                      support_1 = _mm_mul_ps(support_1, eval_1);
                      support_2 = _mm_mul_ps(support_2, eval_2);
                      support_3 = _mm_mul_ps(support_3, eval_3);
                      support_4 = _mm_mul_ps(support_4, eval_4);
                      support_5 = _mm_mul_ps(support_5, eval_5);
                    }
                    // most right basis function on every level
                    else if (ptrIndex[(j * dims) + d] == (ptrLevel[(j * dims) + d] - 1.0f)) {
                      __m128 eval_0 = _mm_load_ps(&(ptrData[(d * result_size) + i]));
                      __m128 eval_1 = _mm_load_ps(&(ptrData[(d * result_size) + i + 4]));
                      __m128 eval_2 = _mm_load_ps(&(ptrData[(d * result_size) + i + 8]));
                      __m128 eval_3 = _mm_load_ps(&(ptrData[(d * result_size) + i + 12]));
                      __m128 eval_4 = _mm_load_ps(&(ptrData[(d * result_size) + i + 16]));
                      __m128 eval_5 = _mm_load_ps(&(ptrData[(d * result_size) + i + 20]));

                      __m128 level = _mm_load1_ps(&(ptrLevel[(j * dims) + d]));
                      __m128 index = _mm_load1_ps(&(ptrIndex[(j * dims) + d]));
#ifdef __FMA4__
                      eval_0 = _mm_msub_ps(eval_0, level, index);
                      eval_1 = _mm_msub_ps(eval_1, level, index);
                      eval_2 = _mm_msub_ps(eval_2, level, index);
                      eval_3 = _mm_msub_ps(eval_3, level, index);
                      eval_4 = _mm_msub_ps(eval_4, level, index);
                      eval_5 = _mm_msub_ps(eval_5, level, index);
#else
                      eval_0 = _mm_sub_ps(_mm_mul_ps(eval_0, level), index);
                      eval_1 = _mm_sub_ps(_mm_mul_ps(eval_1, level), index);
                      eval_2 = _mm_sub_ps(_mm_mul_ps(eval_2, level), index);
                      eval_3 = _mm_sub_ps(_mm_mul_ps(eval_3, level), index);
                      eval_4 = _mm_sub_ps(_mm_mul_ps(eval_4, level), index);
                      eval_5 = _mm_sub_ps(_mm_mul_ps(eval_5, level), index);
#endif
                      eval_0 = _mm_add_ps(one, eval_0);
                      eval_1 = _mm_add_ps(one, eval_1);
                      eval_2 = _mm_add_ps(one, eval_2);
                      eval_3 = _mm_add_ps(one, eval_3);
                      eval_4 = _mm_add_ps(one, eval_4);
                      eval_5 = _mm_add_ps(one, eval_5);

                      eval_0 = _mm_max_ps(zero, eval_0);
                      eval_1 = _mm_max_ps(zero, eval_1);
                      eval_2 = _mm_max_ps(zero, eval_2);
                      eval_3 = _mm_max_ps(zero, eval_3);
                      eval_4 = _mm_max_ps(zero, eval_4);
                      eval_5 = _mm_max_ps(zero, eval_5);

                      support_0 = _mm_mul_ps(support_0, eval_0);
                      support_1 = _mm_mul_ps(support_1, eval_1);
                      support_2 = _mm_mul_ps(support_2, eval_2);
                      support_3 = _mm_mul_ps(support_3, eval_3);
                      support_4 = _mm_mul_ps(support_4, eval_4);
                      support_5 = _mm_mul_ps(support_5, eval_5);
                    }
                    // all other basis functions
                    else {
                      __m128 eval_0 = _mm_load_ps(&(ptrData[(d * result_size) + i]));
                      __m128 eval_1 = _mm_load_ps(&(ptrData[(d * result_size) + i + 4]));
                      __m128 eval_2 = _mm_load_ps(&(ptrData[(d * result_size) + i + 8]));
                      __m128 eval_3 = _mm_load_ps(&(ptrData[(d * result_size) + i + 12]));
                      __m128 eval_4 = _mm_load_ps(&(ptrData[(d * result_size) + i + 16]));
                      __m128 eval_5 = _mm_load_ps(&(ptrData[(d * result_size) + i + 20]));

                      __m128 level = _mm_load1_ps(&(ptrLevel[(j * dims) + d]));
                      __m128 index = _mm_load1_ps(&(ptrIndex[(j * dims) + d]));
#ifdef __FMA4__
                      eval_0 = _mm_msub_ps(eval_0, level, index);
                      eval_1 = _mm_msub_ps(eval_1, level, index);
                      eval_2 = _mm_msub_ps(eval_2, level, index);
                      eval_3 = _mm_msub_ps(eval_3, level, index);
                      eval_4 = _mm_msub_ps(eval_4, level, index);
                      eval_5 = _mm_msub_ps(eval_5, level, index);
#else
                      eval_0 = _mm_sub_ps(_mm_mul_ps(eval_0, level), index);
                      eval_1 = _mm_sub_ps(_mm_mul_ps(eval_1, level), index);
                      eval_2 = _mm_sub_ps(_mm_mul_ps(eval_2, level), index);
                      eval_3 = _mm_sub_ps(_mm_mul_ps(eval_3, level), index);
                      eval_4 = _mm_sub_ps(_mm_mul_ps(eval_4, level), index);
                      eval_5 = _mm_sub_ps(_mm_mul_ps(eval_5, level), index);
#endif
                      __m128 mask = _mm_set1_ps(*fmask);

                      eval_0 = _mm_and_ps(mask, eval_0);
                      eval_1 = _mm_and_ps(mask, eval_1);
                      eval_2 = _mm_and_ps(mask, eval_2);
                      eval_3 = _mm_and_ps(mask, eval_3);
                      eval_4 = _mm_and_ps(mask, eval_4);
                      eval_5 = _mm_and_ps(mask, eval_5);

                      eval_0 = _mm_sub_ps(one, eval_0);
                      eval_1 = _mm_sub_ps(one, eval_1);
                      eval_2 = _mm_sub_ps(one, eval_2);
                      eval_3 = _mm_sub_ps(one, eval_3);
                      eval_4 = _mm_sub_ps(one, eval_4);
                      eval_5 = _mm_sub_ps(one, eval_5);

                      eval_0 = _mm_max_ps(zero, eval_0);
                      eval_1 = _mm_max_ps(zero, eval_1);
                      eval_2 = _mm_max_ps(zero, eval_2);
                      eval_3 = _mm_max_ps(zero, eval_3);
                      eval_4 = _mm_max_ps(zero, eval_4);
                      eval_5 = _mm_max_ps(zero, eval_5);

                      support_0 = _mm_mul_ps(support_0, eval_0);
                      support_1 = _mm_mul_ps(support_1, eval_1);
                      support_2 = _mm_mul_ps(support_2, eval_2);
                      support_3 = _mm_mul_ps(support_3, eval_3);
                      support_4 = _mm_mul_ps(support_4, eval_4);
                      support_5 = _mm_mul_ps(support_5, eval_5);
                    }
                  }

                  __m128 res_0 = _mm_load_ps(&(ptrResult[i]));
                  __m128 res_1 = _mm_load_ps(&(ptrResult[i + 4]));
                  __m128 res_2 = _mm_load_ps(&(ptrResult[i + 8]));
                  __m128 res_3 = _mm_load_ps(&(ptrResult[i + 12]));
                  __m128 res_4 = _mm_load_ps(&(ptrResult[i + 16]));
                  __m128 res_5 = _mm_load_ps(&(ptrResult[i + 20]));

                  res_0 = _mm_add_ps(res_0, support_0);
                  res_1 = _mm_add_ps(res_1, support_1);
                  res_2 = _mm_add_ps(res_2, support_2);
                  res_3 = _mm_add_ps(res_3, support_3);
                  res_4 = _mm_add_ps(res_4, support_4);
                  res_5 = _mm_add_ps(res_5, support_5);

                  _mm_store_ps(&(ptrResult[i]), res_0);
                  _mm_store_ps(&(ptrResult[i + 4]), res_1);
                  _mm_store_ps(&(ptrResult[i + 8]), res_2);
                  _mm_store_ps(&(ptrResult[i + 12]), res_3);
                  _mm_store_ps(&(ptrResult[i + 16]), res_4);
                  _mm_store_ps(&(ptrResult[i + 20]), res_5);
                }
              }

#endif
#if defined(__SSE3__) && defined(__AVX__)
              size_t grid_inc = std::min<size_t>((size_t)getChunkGridPoints(), (end_index_grid - m));

              int imask = 0x7FFFFFFF;
              float* fmask = (float*)&imask;

              for (size_t i = c; i < c + getChunkDataPoints(); i += 48) {
                for (size_t j = m; j < m + grid_inc; j++) {
                  __m256 support_0 = _mm256_broadcast_ss(&(ptrAlpha[j]));
                  __m256 support_1 = _mm256_broadcast_ss(&(ptrAlpha[j]));
                  __m256 support_2 = _mm256_broadcast_ss(&(ptrAlpha[j]));
                  __m256 support_3 = _mm256_broadcast_ss(&(ptrAlpha[j]));
                  __m256 support_4 = _mm256_broadcast_ss(&(ptrAlpha[j]));
                  __m256 support_5 = _mm256_broadcast_ss(&(ptrAlpha[j]));

                  __m256 one = _mm256_set1_ps(1.0f);
                  __m256 two = _mm256_set1_ps(2.0f);
                  __m256 zero = _mm256_set1_ps(0.0f);

                  for (size_t d = 0; d < dims; d++) {
                    // special case for level 1
                    if (ptrLevel[(j * dims) + d] == 2.0f) {
                      // Nothing (multiply by one)
                    }
                    // most left basis function on every level
                    else if (ptrIndex[(j * dims) + d] == 1.0f) {
                      __m256 eval_0 = _mm256_load_ps(&(ptrData[(d * result_size) + i]));
                      __m256 eval_1 = _mm256_load_ps(&(ptrData[(d * result_size) + i + 8]));
                      __m256 eval_2 = _mm256_load_ps(&(ptrData[(d * result_size) + i + 16]));
                      __m256 eval_3 = _mm256_load_ps(&(ptrData[(d * result_size) + i + 24]));
                      __m256 eval_4 = _mm256_load_ps(&(ptrData[(d * result_size) + i + 32]));
                      __m256 eval_5 = _mm256_load_ps(&(ptrData[(d * result_size) + i + 40]));

                      __m256 level = _mm256_broadcast_ss(&(ptrLevel[(j * dims) + d]));

                      eval_0 = _mm256_mul_ps(eval_0, level);
                      eval_1 = _mm256_mul_ps(eval_1, level);
                      eval_2 = _mm256_mul_ps(eval_2, level);
                      eval_3 = _mm256_mul_ps(eval_3, level);
                      eval_4 = _mm256_mul_ps(eval_4, level);
                      eval_5 = _mm256_mul_ps(eval_5, level);

                      eval_0 = _mm256_sub_ps(two, eval_0);
                      eval_1 = _mm256_sub_ps(two, eval_1);
                      eval_2 = _mm256_sub_ps(two, eval_2);
                      eval_3 = _mm256_sub_ps(two, eval_3);
                      eval_4 = _mm256_sub_ps(two, eval_4);
                      eval_5 = _mm256_sub_ps(two, eval_5);

                      eval_0 = _mm256_max_ps(zero, eval_0);
                      eval_1 = _mm256_max_ps(zero, eval_1);
                      eval_2 = _mm256_max_ps(zero, eval_2);
                      eval_3 = _mm256_max_ps(zero, eval_3);
                      eval_4 = _mm256_max_ps(zero, eval_4);
                      eval_5 = _mm256_max_ps(zero, eval_5);

                      support_0 = _mm256_mul_ps(support_0, eval_0);
                      support_1 = _mm256_mul_ps(support_1, eval_1);
                      support_2 = _mm256_mul_ps(support_2, eval_2);
                      support_3 = _mm256_mul_ps(support_3, eval_3);
                      support_4 = _mm256_mul_ps(support_4, eval_4);
                      support_5 = _mm256_mul_ps(support_5, eval_5);
                    }
                    // most right basis function on every level
                    else if (ptrIndex[(j * dims) + d] == (ptrLevel[(j * dims) + d] - 1.0f)) {
                      __m256 eval_0 = _mm256_load_ps(&(ptrData[(d * result_size) + i]));
                      __m256 eval_1 = _mm256_load_ps(&(ptrData[(d * result_size) + i + 8]));
                      __m256 eval_2 = _mm256_load_ps(&(ptrData[(d * result_size) + i + 16]));
                      __m256 eval_3 = _mm256_load_ps(&(ptrData[(d * result_size) + i + 24]));
                      __m256 eval_4 = _mm256_load_ps(&(ptrData[(d * result_size) + i + 32]));
                      __m256 eval_5 = _mm256_load_ps(&(ptrData[(d * result_size) + i + 40]));

                      __m256 level = _mm256_broadcast_ss(&(ptrLevel[(j * dims) + d]));
                      __m256 index = _mm256_broadcast_ss(&(ptrIndex[(j * dims) + d]));
#ifdef __FMA4__
                      eval_0 = _mm256_msub_ps(eval_0, level, index);
                      eval_1 = _mm256_msub_ps(eval_1, level, index);
                      eval_2 = _mm256_msub_ps(eval_2, level, index);
                      eval_3 = _mm256_msub_ps(eval_3, level, index);
                      eval_4 = _mm256_msub_ps(eval_4, level, index);
                      eval_5 = _mm256_msub_ps(eval_5, level, index);
#else
                      eval_0 = _mm256_sub_ps(_mm256_mul_ps(eval_0, level), index);
                      eval_1 = _mm256_sub_ps(_mm256_mul_ps(eval_1, level), index);
                      eval_2 = _mm256_sub_ps(_mm256_mul_ps(eval_2, level), index);
                      eval_3 = _mm256_sub_ps(_mm256_mul_ps(eval_3, level), index);
                      eval_4 = _mm256_sub_ps(_mm256_mul_ps(eval_4, level), index);
                      eval_5 = _mm256_sub_ps(_mm256_mul_ps(eval_5, level), index);
#endif
                      eval_0 = _mm256_add_ps(one, eval_0);
                      eval_1 = _mm256_add_ps(one, eval_1);
                      eval_2 = _mm256_add_ps(one, eval_2);
                      eval_3 = _mm256_add_ps(one, eval_3);
                      eval_4 = _mm256_add_ps(one, eval_4);
                      eval_5 = _mm256_add_ps(one, eval_5);

                      eval_0 = _mm256_max_ps(zero, eval_0);
                      eval_1 = _mm256_max_ps(zero, eval_1);
                      eval_2 = _mm256_max_ps(zero, eval_2);
                      eval_3 = _mm256_max_ps(zero, eval_3);
                      eval_4 = _mm256_max_ps(zero, eval_4);
                      eval_5 = _mm256_max_ps(zero, eval_5);

                      support_0 = _mm256_mul_ps(support_0, eval_0);
                      support_1 = _mm256_mul_ps(support_1, eval_1);
                      support_2 = _mm256_mul_ps(support_2, eval_2);
                      support_3 = _mm256_mul_ps(support_3, eval_3);
                      support_4 = _mm256_mul_ps(support_4, eval_4);
                      support_5 = _mm256_mul_ps(support_5, eval_5);
                    }
                    // all other basis functions
                    else {
                      __m256 eval_0 = _mm256_load_ps(&(ptrData[(d * result_size) + i]));
                      __m256 eval_1 = _mm256_load_ps(&(ptrData[(d * result_size) + i + 8]));
                      __m256 eval_2 = _mm256_load_ps(&(ptrData[(d * result_size) + i + 16]));
                      __m256 eval_3 = _mm256_load_ps(&(ptrData[(d * result_size) + i + 24]));
                      __m256 eval_4 = _mm256_load_ps(&(ptrData[(d * result_size) + i + 32]));
                      __m256 eval_5 = _mm256_load_ps(&(ptrData[(d * result_size) + i + 40]));

                      __m256 level = _mm256_broadcast_ss(&(ptrLevel[(j * dims) + d]));
                      __m256 index = _mm256_broadcast_ss(&(ptrIndex[(j * dims) + d]));
#ifdef __FMA4__
                      eval_0 = _mm256_msub_ps(eval_0, level, index);
                      eval_1 = _mm256_msub_ps(eval_1, level, index);
                      eval_2 = _mm256_msub_ps(eval_2, level, index);
                      eval_3 = _mm256_msub_ps(eval_3, level, index);
                      eval_4 = _mm256_msub_ps(eval_4, level, index);
                      eval_5 = _mm256_msub_ps(eval_5, level, index);
#else
                      eval_0 = _mm256_sub_ps(_mm256_mul_ps(eval_0, level), index);
                      eval_1 = _mm256_sub_ps(_mm256_mul_ps(eval_1, level), index);
                      eval_2 = _mm256_sub_ps(_mm256_mul_ps(eval_2, level), index);
                      eval_3 = _mm256_sub_ps(_mm256_mul_ps(eval_3, level), index);
                      eval_4 = _mm256_sub_ps(_mm256_mul_ps(eval_4, level), index);
                      eval_5 = _mm256_sub_ps(_mm256_mul_ps(eval_5, level), index);
#endif
                      __m256 mask = _mm256_set1_ps(*fmask);

                      eval_0 = _mm256_and_ps(mask, eval_0);
                      eval_1 = _mm256_and_ps(mask, eval_1);
                      eval_2 = _mm256_and_ps(mask, eval_2);
                      eval_3 = _mm256_and_ps(mask, eval_3);
                      eval_4 = _mm256_and_ps(mask, eval_4);
                      eval_5 = _mm256_and_ps(mask, eval_5);

                      eval_0 = _mm256_sub_ps(one, eval_0);
                      eval_1 = _mm256_sub_ps(one, eval_1);
                      eval_2 = _mm256_sub_ps(one, eval_2);
                      eval_3 = _mm256_sub_ps(one, eval_3);
                      eval_4 = _mm256_sub_ps(one, eval_4);
                      eval_5 = _mm256_sub_ps(one, eval_5);

                      eval_0 = _mm256_max_ps(zero, eval_0);
                      eval_1 = _mm256_max_ps(zero, eval_1);
                      eval_2 = _mm256_max_ps(zero, eval_2);
                      eval_3 = _mm256_max_ps(zero, eval_3);
                      eval_4 = _mm256_max_ps(zero, eval_4);
                      eval_5 = _mm256_max_ps(zero, eval_5);

                      support_0 = _mm256_mul_ps(support_0, eval_0);
                      support_1 = _mm256_mul_ps(support_1, eval_1);
                      support_2 = _mm256_mul_ps(support_2, eval_2);
                      support_3 = _mm256_mul_ps(support_3, eval_3);
                      support_4 = _mm256_mul_ps(support_4, eval_4);
                      support_5 = _mm256_mul_ps(support_5, eval_5);
                    }
                  }

                  __m256 res_0 = _mm256_load_ps(&(ptrResult[i]));
                  __m256 res_1 = _mm256_load_ps(&(ptrResult[i + 8]));
                  __m256 res_2 = _mm256_load_ps(&(ptrResult[i + 16]));
                  __m256 res_3 = _mm256_load_ps(&(ptrResult[i + 24]));
                  __m256 res_4 = _mm256_load_ps(&(ptrResult[i + 32]));
                  __m256 res_5 = _mm256_load_ps(&(ptrResult[i + 40]));

                  res_0 = _mm256_add_ps(res_0, support_0);
                  res_1 = _mm256_add_ps(res_1, support_1);
                  res_2 = _mm256_add_ps(res_2, support_2);
                  res_3 = _mm256_add_ps(res_3, support_3);
                  res_4 = _mm256_add_ps(res_4, support_4);
                  res_5 = _mm256_add_ps(res_5, support_5);

                  _mm256_store_ps(&(ptrResult[i]), res_0);
                  _mm256_store_ps(&(ptrResult[i + 8]), res_1);
                  _mm256_store_ps(&(ptrResult[i + 16]), res_2);
                  _mm256_store_ps(&(ptrResult[i + 24]), res_3);
                  _mm256_store_ps(&(ptrResult[i + 32]), res_4);
                  _mm256_store_ps(&(ptrResult[i + 40]), res_5);
                }
              }

#endif
#if !defined(__SSE3__) && !defined(__AVX__)
              size_t grid_end = std::min<size_t>((size_t)getChunkGridPoints() + m, end_index_grid);

              for (size_t i = c; i < data_end; i++) {
                for (size_t j = m; j < grid_end; j++) {
                  float curSupport = ptrAlpha[j];

                  for (size_t d = 0; d < dims; d++) {
                    if (ptrLevel[(j * dims) + d] == 2.0f) {
                      // nothing to do (mult with 1)
                    } else if (ptrIndex[(j * dims) + d] == 1.0f) {
                      float eval = ((ptrLevel[(j * dims) + d]) * (ptrData[(d * result_size) + i]));
                      eval = 2.0f - eval;
                      float localSupport = std::max<float>(eval, 0.0f);
                      curSupport *= localSupport;
                    } else if (ptrIndex[(j * dims) + d] == (ptrLevel[(j * dims) + d] - 1.0f)) {
                      float eval = ((ptrLevel[(j * dims) + d]) * (ptrData[(d * result_size) + i]));
                      float index_calc = eval - (ptrIndex[(j * dims) + d]);
                      float last = 1.0f + index_calc;
                      float localSupport = std::max<float>(last, 0.0f);
                      curSupport *= localSupport;
                    } else {
                      float eval = ((ptrLevel[(j * dims) + d]) * (ptrData[(d * result_size) + i]));
                      float index_calc = eval - (ptrIndex[(j * dims) + d]);
                      float abs = static_cast<float>(fabs(index_calc));
                      float last = 1.0f - abs;
                      float localSupport = std::max<float>(last, 0.0f);
                      curSupport *= localSupport;
                    }
                  }

                  ptrResult[i] += curSupport;
                }
              }

#endif
            }
          }
        }

        static inline void multTransposeImpl(
          SGPP::base::DataMatrixSP* level,
          SGPP::base::DataMatrixSP* index,
          SGPP::base::DataMatrixSP* /*mask*/, //unused for this specialization
          SGPP::base::DataMatrixSP* /*offset*/, //unused for this specialization
          SGPP::base::DataMatrixSP* dataset,
          SGPP::base::DataVectorSP& source,
          SGPP::base::DataVectorSP& result,
          const size_t start_index_grid,
          const size_t end_index_grid,
          const size_t start_index_data,
          const size_t end_index_data) {
          float* ptrLevel = level->getPointer();
          float* ptrIndex = index->getPointer();
          float* ptrSource = source.getPointer();
          float* ptrData = dataset->getPointer();
          float* ptrResult = result.getPointer();
          size_t source_size = source.getSize();
          size_t dims = dataset->getNrows();

          CHECK_ARGS_MULTTRANSPOSE(level, dataset, source, start_index_grid, end_index_grid, start_index_data, end_index_data);

          for (size_t k = start_index_grid; k < end_index_grid; k += std::min<size_t>((size_t)getChunkGridPoints(), (end_index_grid - k))) {
            size_t grid_inc = std::min<size_t>((size_t)getChunkGridPoints(), (end_index_grid - k));
#if defined(__SSE3__) && !defined(__AVX__)
            int imask = 0x7FFFFFFF;
            float* fmask = (float*)&imask;

            for (size_t i = start_index_data; i < end_index_data; i += 24) {
              for (size_t j = k; j < k + grid_inc; j++) {
                __m128 support_0 = _mm_load_ps(&(ptrSource[i]));
                __m128 support_1 = _mm_load_ps(&(ptrSource[i + 4]));
                __m128 support_2 = _mm_load_ps(&(ptrSource[i + 8]));
                __m128 support_3 = _mm_load_ps(&(ptrSource[i + 12]));
                __m128 support_4 = _mm_load_ps(&(ptrSource[i + 16]));
                __m128 support_5 = _mm_load_ps(&(ptrSource[i + 20]));

                __m128 one = _mm_set1_ps(1.0f);
                __m128 two = _mm_set1_ps(2.0f);
                __m128 zero = _mm_set1_ps(0.0f);

                for (size_t d = 0; d < dims; d++) {
                  // special case for level 1
                  if (ptrLevel[(j * dims) + d] == 2.0f) {
                    // Nothing (multiply by one)
                  }
                  // most left basis function on every level
                  else if (ptrIndex[(j * dims) + d] == 1.0f) {
                    __m128 eval_0 = _mm_load_ps(&(ptrData[(d * source_size) + i]));
                    __m128 eval_1 = _mm_load_ps(&(ptrData[(d * source_size) + i + 4]));
                    __m128 eval_2 = _mm_load_ps(&(ptrData[(d * source_size) + i + 8]));
                    __m128 eval_3 = _mm_load_ps(&(ptrData[(d * source_size) + i + 12]));
                    __m128 eval_4 = _mm_load_ps(&(ptrData[(d * source_size) + i + 16]));
                    __m128 eval_5 = _mm_load_ps(&(ptrData[(d * source_size) + i + 20]));

                    __m128 level = _mm_load1_ps(&(ptrLevel[(j * dims) + d]));

                    eval_0 = _mm_mul_ps(eval_0, level);
                    eval_1 = _mm_mul_ps(eval_1, level);
                    eval_2 = _mm_mul_ps(eval_2, level);
                    eval_3 = _mm_mul_ps(eval_3, level);
                    eval_4 = _mm_mul_ps(eval_4, level);
                    eval_5 = _mm_mul_ps(eval_5, level);

                    eval_0 = _mm_sub_ps(two, eval_0);
                    eval_1 = _mm_sub_ps(two, eval_1);
                    eval_2 = _mm_sub_ps(two, eval_2);
                    eval_3 = _mm_sub_ps(two, eval_3);
                    eval_4 = _mm_sub_ps(two, eval_4);
                    eval_5 = _mm_sub_ps(two, eval_5);

                    eval_0 = _mm_max_ps(zero, eval_0);
                    eval_1 = _mm_max_ps(zero, eval_1);
                    eval_2 = _mm_max_ps(zero, eval_2);
                    eval_3 = _mm_max_ps(zero, eval_3);
                    eval_4 = _mm_max_ps(zero, eval_4);
                    eval_5 = _mm_max_ps(zero, eval_5);

                    support_0 = _mm_mul_ps(support_0, eval_0);
                    support_1 = _mm_mul_ps(support_1, eval_1);
                    support_2 = _mm_mul_ps(support_2, eval_2);
                    support_3 = _mm_mul_ps(support_3, eval_3);
                    support_4 = _mm_mul_ps(support_4, eval_4);
                    support_5 = _mm_mul_ps(support_5, eval_5);
                  }
                  // most right basis function on every level
                  else if (ptrIndex[(j * dims) + d] == (ptrLevel[(j * dims) + d] - 1.0f)) {
                    __m128 eval_0 = _mm_load_ps(&(ptrData[(d * source_size) + i]));
                    __m128 eval_1 = _mm_load_ps(&(ptrData[(d * source_size) + i + 4]));
                    __m128 eval_2 = _mm_load_ps(&(ptrData[(d * source_size) + i + 8]));
                    __m128 eval_3 = _mm_load_ps(&(ptrData[(d * source_size) + i + 12]));
                    __m128 eval_4 = _mm_load_ps(&(ptrData[(d * source_size) + i + 16]));
                    __m128 eval_5 = _mm_load_ps(&(ptrData[(d * source_size) + i + 20]));

                    __m128 level = _mm_load1_ps(&(ptrLevel[(j * dims) + d]));
                    __m128 index = _mm_load1_ps(&(ptrIndex[(j * dims) + d]));
#ifdef __FMA4__
                    eval_0 = _mm_msub_ps(eval_0, level, index);
                    eval_1 = _mm_msub_ps(eval_1, level, index);
                    eval_2 = _mm_msub_ps(eval_2, level, index);
                    eval_3 = _mm_msub_ps(eval_3, level, index);
                    eval_4 = _mm_msub_ps(eval_4, level, index);
                    eval_5 = _mm_msub_ps(eval_5, level, index);
#else
                    eval_0 = _mm_sub_ps(_mm_mul_ps(eval_0, level), index);
                    eval_1 = _mm_sub_ps(_mm_mul_ps(eval_1, level), index);
                    eval_2 = _mm_sub_ps(_mm_mul_ps(eval_2, level), index);
                    eval_3 = _mm_sub_ps(_mm_mul_ps(eval_3, level), index);
                    eval_4 = _mm_sub_ps(_mm_mul_ps(eval_4, level), index);
                    eval_5 = _mm_sub_ps(_mm_mul_ps(eval_5, level), index);
#endif
                    eval_0 = _mm_add_ps(one, eval_0);
                    eval_1 = _mm_add_ps(one, eval_1);
                    eval_2 = _mm_add_ps(one, eval_2);
                    eval_3 = _mm_add_ps(one, eval_3);
                    eval_4 = _mm_add_ps(one, eval_4);
                    eval_5 = _mm_add_ps(one, eval_5);

                    eval_0 = _mm_max_ps(zero, eval_0);
                    eval_1 = _mm_max_ps(zero, eval_1);
                    eval_2 = _mm_max_ps(zero, eval_2);
                    eval_3 = _mm_max_ps(zero, eval_3);
                    eval_4 = _mm_max_ps(zero, eval_4);
                    eval_5 = _mm_max_ps(zero, eval_5);

                    support_0 = _mm_mul_ps(support_0, eval_0);
                    support_1 = _mm_mul_ps(support_1, eval_1);
                    support_2 = _mm_mul_ps(support_2, eval_2);
                    support_3 = _mm_mul_ps(support_3, eval_3);
                    support_4 = _mm_mul_ps(support_4, eval_4);
                    support_5 = _mm_mul_ps(support_5, eval_5);
                  }
                  // all other basis functions
                  else {
                    __m128 eval_0 = _mm_load_ps(&(ptrData[(d * source_size) + i]));
                    __m128 eval_1 = _mm_load_ps(&(ptrData[(d * source_size) + i + 4]));
                    __m128 eval_2 = _mm_load_ps(&(ptrData[(d * source_size) + i + 8]));
                    __m128 eval_3 = _mm_load_ps(&(ptrData[(d * source_size) + i + 12]));
                    __m128 eval_4 = _mm_load_ps(&(ptrData[(d * source_size) + i + 16]));
                    __m128 eval_5 = _mm_load_ps(&(ptrData[(d * source_size) + i + 20]));

                    __m128 level = _mm_load1_ps(&(ptrLevel[(j * dims) + d]));
                    __m128 index = _mm_load1_ps(&(ptrIndex[(j * dims) + d]));
#ifdef __FMA4__
                    eval_0 = _mm_msub_ps(eval_0, level, index);
                    eval_1 = _mm_msub_ps(eval_1, level, index);
                    eval_2 = _mm_msub_ps(eval_2, level, index);
                    eval_3 = _mm_msub_ps(eval_3, level, index);
                    eval_4 = _mm_msub_ps(eval_4, level, index);
                    eval_5 = _mm_msub_ps(eval_5, level, index);
#else
                    eval_0 = _mm_sub_ps(_mm_mul_ps(eval_0, level), index);
                    eval_1 = _mm_sub_ps(_mm_mul_ps(eval_1, level), index);
                    eval_2 = _mm_sub_ps(_mm_mul_ps(eval_2, level), index);
                    eval_3 = _mm_sub_ps(_mm_mul_ps(eval_3, level), index);
                    eval_4 = _mm_sub_ps(_mm_mul_ps(eval_4, level), index);
                    eval_5 = _mm_sub_ps(_mm_mul_ps(eval_5, level), index);
#endif
                    __m128 mask = _mm_set1_ps(*fmask);

                    eval_0 = _mm_and_ps(mask, eval_0);
                    eval_1 = _mm_and_ps(mask, eval_1);
                    eval_2 = _mm_and_ps(mask, eval_2);
                    eval_3 = _mm_and_ps(mask, eval_3);
                    eval_4 = _mm_and_ps(mask, eval_4);
                    eval_5 = _mm_and_ps(mask, eval_5);

                    eval_0 = _mm_sub_ps(one, eval_0);
                    eval_1 = _mm_sub_ps(one, eval_1);
                    eval_2 = _mm_sub_ps(one, eval_2);
                    eval_3 = _mm_sub_ps(one, eval_3);
                    eval_4 = _mm_sub_ps(one, eval_4);
                    eval_5 = _mm_sub_ps(one, eval_5);

                    eval_0 = _mm_max_ps(zero, eval_0);
                    eval_1 = _mm_max_ps(zero, eval_1);
                    eval_2 = _mm_max_ps(zero, eval_2);
                    eval_3 = _mm_max_ps(zero, eval_3);
                    eval_4 = _mm_max_ps(zero, eval_4);
                    eval_5 = _mm_max_ps(zero, eval_5);

                    support_0 = _mm_mul_ps(support_0, eval_0);
                    support_1 = _mm_mul_ps(support_1, eval_1);
                    support_2 = _mm_mul_ps(support_2, eval_2);
                    support_3 = _mm_mul_ps(support_3, eval_3);
                    support_4 = _mm_mul_ps(support_4, eval_4);
                    support_5 = _mm_mul_ps(support_5, eval_5);
                  }
                }

                __m128 res_0 = _mm_setzero_ps();
                res_0 = _mm_load_ss(&(ptrResult[j]));

                support_0 = _mm_add_ps(support_0, support_1);
                support_2 = _mm_add_ps(support_2, support_3);
                support_4 = _mm_add_ps(support_4, support_5);
                support_0 = _mm_add_ps(support_0, support_2);
                support_0 = _mm_add_ps(support_0, support_4);

                support_0 = _mm_hadd_ps(support_0, support_0);
                support_0 = _mm_hadd_ps(support_0, support_0);
                res_0 = _mm_add_ss(res_0, support_0);

                _mm_store_ss(&(ptrResult[j]), res_0);
              }
            }

#endif
#if defined(__SSE3__) && defined(__AVX__)
            int imask = 0x7FFFFFFF;
            float* fmask = (float*)&imask;

            for (size_t i = start_index_data; i < end_index_data; i += 48) {
              for (size_t j = k; j < k + grid_inc; j++) {
                __m256 support_0 = _mm256_load_ps(&(ptrSource[i]));
                __m256 support_1 = _mm256_load_ps(&(ptrSource[i + 8]));
                __m256 support_2 = _mm256_load_ps(&(ptrSource[i + 16]));
                __m256 support_3 = _mm256_load_ps(&(ptrSource[i + 24]));
                __m256 support_4 = _mm256_load_ps(&(ptrSource[i + 32]));
                __m256 support_5 = _mm256_load_ps(&(ptrSource[i + 40]));

                __m256 one = _mm256_set1_ps(1.0f);
                __m256 two = _mm256_set1_ps(2.0f);
                __m256 zero = _mm256_set1_ps(0.0f);

                for (size_t d = 0; d < dims; d++) {
                  // special case for level 1
                  if (ptrLevel[(j * dims) + d] == 2.0f) {
                    // Nothing (multiply by one)
                  }
                  // most left basis function on every level
                  else if (ptrIndex[(j * dims) + d] == 1.0f) {
                    __m256 eval_0 = _mm256_load_ps(&(ptrData[(d * source_size) + i]));
                    __m256 eval_1 = _mm256_load_ps(&(ptrData[(d * source_size) + i + 8]));
                    __m256 eval_2 = _mm256_load_ps(&(ptrData[(d * source_size) + i + 16]));
                    __m256 eval_3 = _mm256_load_ps(&(ptrData[(d * source_size) + i + 24]));
                    __m256 eval_4 = _mm256_load_ps(&(ptrData[(d * source_size) + i + 32]));
                    __m256 eval_5 = _mm256_load_ps(&(ptrData[(d * source_size) + i + 40]));

                    __m256 level = _mm256_broadcast_ss(&(ptrLevel[(j * dims) + d]));

                    eval_0 = _mm256_mul_ps(eval_0, level);
                    eval_1 = _mm256_mul_ps(eval_1, level);
                    eval_2 = _mm256_mul_ps(eval_2, level);
                    eval_3 = _mm256_mul_ps(eval_3, level);
                    eval_4 = _mm256_mul_ps(eval_4, level);
                    eval_5 = _mm256_mul_ps(eval_5, level);

                    eval_0 = _mm256_sub_ps(two, eval_0);
                    eval_1 = _mm256_sub_ps(two, eval_1);
                    eval_2 = _mm256_sub_ps(two, eval_2);
                    eval_3 = _mm256_sub_ps(two, eval_3);
                    eval_4 = _mm256_sub_ps(two, eval_4);
                    eval_5 = _mm256_sub_ps(two, eval_5);

                    eval_0 = _mm256_max_ps(zero, eval_0);
                    eval_1 = _mm256_max_ps(zero, eval_1);
                    eval_2 = _mm256_max_ps(zero, eval_2);
                    eval_3 = _mm256_max_ps(zero, eval_3);
                    eval_4 = _mm256_max_ps(zero, eval_4);
                    eval_5 = _mm256_max_ps(zero, eval_5);

                    support_0 = _mm256_mul_ps(support_0, eval_0);
                    support_1 = _mm256_mul_ps(support_1, eval_1);
                    support_2 = _mm256_mul_ps(support_2, eval_2);
                    support_3 = _mm256_mul_ps(support_3, eval_3);
                    support_4 = _mm256_mul_ps(support_4, eval_4);
                    support_5 = _mm256_mul_ps(support_5, eval_5);
                  }
                  // most right basis function on every level
                  else if (ptrIndex[(j * dims) + d] == (ptrLevel[(j * dims) + d] - 1.0f)) {
                    __m256 eval_0 = _mm256_load_ps(&(ptrData[(d * source_size) + i]));
                    __m256 eval_1 = _mm256_load_ps(&(ptrData[(d * source_size) + i + 8]));
                    __m256 eval_2 = _mm256_load_ps(&(ptrData[(d * source_size) + i + 16]));
                    __m256 eval_3 = _mm256_load_ps(&(ptrData[(d * source_size) + i + 24]));
                    __m256 eval_4 = _mm256_load_ps(&(ptrData[(d * source_size) + i + 32]));
                    __m256 eval_5 = _mm256_load_ps(&(ptrData[(d * source_size) + i + 40]));

                    __m256 level = _mm256_broadcast_ss(&(ptrLevel[(j * dims) + d]));
                    __m256 index = _mm256_broadcast_ss(&(ptrIndex[(j * dims) + d]));
#ifdef __FMA4__
                    eval_0 = _mm256_msub_ps(eval_0, level, index);
                    eval_1 = _mm256_msub_ps(eval_1, level, index);
                    eval_2 = _mm256_msub_ps(eval_2, level, index);
                    eval_3 = _mm256_msub_ps(eval_3, level, index);
                    eval_4 = _mm256_msub_ps(eval_4, level, index);
                    eval_5 = _mm256_msub_ps(eval_5, level, index);
#else
                    eval_0 = _mm256_sub_ps(_mm256_mul_ps(eval_0, level), index);
                    eval_1 = _mm256_sub_ps(_mm256_mul_ps(eval_1, level), index);
                    eval_2 = _mm256_sub_ps(_mm256_mul_ps(eval_2, level), index);
                    eval_3 = _mm256_sub_ps(_mm256_mul_ps(eval_3, level), index);
                    eval_4 = _mm256_sub_ps(_mm256_mul_ps(eval_4, level), index);
                    eval_5 = _mm256_sub_ps(_mm256_mul_ps(eval_5, level), index);
#endif
                    eval_0 = _mm256_add_ps(one, eval_0);
                    eval_1 = _mm256_add_ps(one, eval_1);
                    eval_2 = _mm256_add_ps(one, eval_2);
                    eval_3 = _mm256_add_ps(one, eval_3);
                    eval_4 = _mm256_add_ps(one, eval_4);
                    eval_5 = _mm256_add_ps(one, eval_5);

                    eval_0 = _mm256_max_ps(zero, eval_0);
                    eval_1 = _mm256_max_ps(zero, eval_1);
                    eval_2 = _mm256_max_ps(zero, eval_2);
                    eval_3 = _mm256_max_ps(zero, eval_3);
                    eval_4 = _mm256_max_ps(zero, eval_4);
                    eval_5 = _mm256_max_ps(zero, eval_5);

                    support_0 = _mm256_mul_ps(support_0, eval_0);
                    support_1 = _mm256_mul_ps(support_1, eval_1);
                    support_2 = _mm256_mul_ps(support_2, eval_2);
                    support_3 = _mm256_mul_ps(support_3, eval_3);
                    support_4 = _mm256_mul_ps(support_4, eval_4);
                    support_5 = _mm256_mul_ps(support_5, eval_5);
                  }
                  // all other basis functions
                  else {
                    __m256 eval_0 = _mm256_load_ps(&(ptrData[(d * source_size) + i]));
                    __m256 eval_1 = _mm256_load_ps(&(ptrData[(d * source_size) + i + 8]));
                    __m256 eval_2 = _mm256_load_ps(&(ptrData[(d * source_size) + i + 16]));
                    __m256 eval_3 = _mm256_load_ps(&(ptrData[(d * source_size) + i + 24]));
                    __m256 eval_4 = _mm256_load_ps(&(ptrData[(d * source_size) + i + 32]));
                    __m256 eval_5 = _mm256_load_ps(&(ptrData[(d * source_size) + i + 40]));

                    __m256 level = _mm256_broadcast_ss(&(ptrLevel[(j * dims) + d]));
                    __m256 index = _mm256_broadcast_ss(&(ptrIndex[(j * dims) + d]));
#ifdef __FMA4__
                    eval_0 = _mm256_msub_ps(eval_0, level, index);
                    eval_1 = _mm256_msub_ps(eval_1, level, index);
                    eval_2 = _mm256_msub_ps(eval_2, level, index);
                    eval_3 = _mm256_msub_ps(eval_3, level, index);
                    eval_4 = _mm256_msub_ps(eval_4, level, index);
                    eval_5 = _mm256_msub_ps(eval_5, level, index);
#else
                    eval_0 = _mm256_sub_ps(_mm256_mul_ps(eval_0, level), index);
                    eval_1 = _mm256_sub_ps(_mm256_mul_ps(eval_1, level), index);
                    eval_2 = _mm256_sub_ps(_mm256_mul_ps(eval_2, level), index);
                    eval_3 = _mm256_sub_ps(_mm256_mul_ps(eval_3, level), index);
                    eval_4 = _mm256_sub_ps(_mm256_mul_ps(eval_4, level), index);
                    eval_5 = _mm256_sub_ps(_mm256_mul_ps(eval_5, level), index);
#endif
                    __m256 mask = _mm256_set1_ps(*fmask);

                    eval_0 = _mm256_and_ps(mask, eval_0);
                    eval_1 = _mm256_and_ps(mask, eval_1);
                    eval_2 = _mm256_and_ps(mask, eval_2);
                    eval_3 = _mm256_and_ps(mask, eval_3);
                    eval_4 = _mm256_and_ps(mask, eval_4);
                    eval_5 = _mm256_and_ps(mask, eval_5);

                    eval_0 = _mm256_sub_ps(one, eval_0);
                    eval_1 = _mm256_sub_ps(one, eval_1);
                    eval_2 = _mm256_sub_ps(one, eval_2);
                    eval_3 = _mm256_sub_ps(one, eval_3);
                    eval_4 = _mm256_sub_ps(one, eval_4);
                    eval_5 = _mm256_sub_ps(one, eval_5);

                    eval_0 = _mm256_max_ps(zero, eval_0);
                    eval_1 = _mm256_max_ps(zero, eval_1);
                    eval_2 = _mm256_max_ps(zero, eval_2);
                    eval_3 = _mm256_max_ps(zero, eval_3);
                    eval_4 = _mm256_max_ps(zero, eval_4);
                    eval_5 = _mm256_max_ps(zero, eval_5);

                    support_0 = _mm256_mul_ps(support_0, eval_0);
                    support_1 = _mm256_mul_ps(support_1, eval_1);
                    support_2 = _mm256_mul_ps(support_2, eval_2);
                    support_3 = _mm256_mul_ps(support_3, eval_3);
                    support_4 = _mm256_mul_ps(support_4, eval_4);
                    support_5 = _mm256_mul_ps(support_5, eval_5);
                  }
                }

                const __m256i ldStMaskSPAVX = _mm256_set_epi32(0x00000000, 0x00000000, 0x00000000, 0x00000000, 0x00000000, 0x00000000, 0x00000000, 0xFFFFFFFF);

                support_0 = _mm256_add_ps(support_0, support_1);
                support_2 = _mm256_add_ps(support_2, support_3);
                support_4 = _mm256_add_ps(support_4, support_5);
                support_0 = _mm256_add_ps(support_0, support_2);
                support_0 = _mm256_add_ps(support_0, support_4);

                support_0 = _mm256_hadd_ps(support_0, support_0);
                __m256 tmp = _mm256_permute2f128_ps(support_0, support_0, 0x81);
                support_0 = _mm256_add_ps(support_0, tmp);
                support_0 = _mm256_hadd_ps(support_0, support_0);

                // Workaround: bug with maskload in GCC (4.6.1)
#ifdef __ICC
                __m256 res_0 = _mm256_maskload_ps(&(ptrResult[j]), ldStMaskSPAVX);
                res_0 = _mm256_add_ps(res_0, support_0);
                _mm256_maskstore_ps(&(ptrResult[j]), ldStMaskSPAVX, res_0);
#else
                float tmp_reduce;
                _mm256_maskstore_ps(&(tmp_reduce), ldStMaskSPAVX, support_0);
                ptrResult[j] += tmp_reduce;
#endif
              }
            }

#endif
#if !defined(__SSE3__) && !defined(__AVX__)

            for (size_t i = start_index_data; i < end_index_data; i++) {
              for (size_t j = k; j < k + grid_inc; j++) {
                float curSupport = ptrSource[i];

                for (size_t d = 0; d < dims; d++) {
                  if (ptrLevel[(j * dims) + d] == 2.0f) {
                    // nothing to do (mult with 1)
                  } else if (ptrIndex[(j * dims) + d] == 1.0f) {
                    float eval = ((ptrLevel[(j * dims) + d]) * (ptrData[(d * source_size) + i]));
                    eval = 2.0f - eval;
                    float localSupport = std::max<float>(eval, 0.0f);
                    curSupport *= localSupport;
                  } else if (ptrIndex[(j * dims) + d] == (ptrLevel[(j * dims) + d] - 1.0f)) {
                    float eval = ((ptrLevel[(j * dims) + d]) * (ptrData[(d * source_size) + i]));
                    float index_calc = eval - (ptrIndex[(j * dims) + d]);
                    float last = 1.0f + index_calc;
                    float localSupport = std::max<float>(last, 0.0f);
                    curSupport *= localSupport;
                  } else {
                    float eval = ((ptrLevel[(j * dims) + d]) * (ptrData[(d * source_size) + i]));
                    float index_calc = eval - (ptrIndex[(j * dims) + d]);
                    float abs = static_cast<float>(fabs(index_calc));
                    float last = 1.0f - abs;
                    float localSupport = std::max<float>(last, 0.0f);
                    curSupport *= localSupport;
                  }
                }

                ptrResult[j] += curSupport;
              }
            }

#endif
          }
        }
    };

  }
}

#endif // SPX86SIMDMODLINEAR_HPP