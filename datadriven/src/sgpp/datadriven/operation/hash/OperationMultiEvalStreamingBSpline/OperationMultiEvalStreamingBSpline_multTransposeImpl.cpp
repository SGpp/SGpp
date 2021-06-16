// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <sgpp/datadriven/operation/hash/OperationMultiEvalStreamingBSpline/OperationMultiEvalStreamingBSpline.hpp>

#include <sgpp/globaldef.hpp>

#if defined(__SSE3__) && !defined(__AVX__)
#include <pmmintrin.h>
#endif

#if defined(__SSE3__) && defined(__AVX__)
#include <immintrin.h>
#endif

#if defined(__MIC__)
#include <immintrin.h>
#endif

namespace sgpp {
namespace datadriven {

void OperationMultiEvalStreamingBSpline::multTransposeImpl(
    sgpp::base::DataMatrix* level, sgpp::base::DataMatrix* index, sgpp::base::DataMatrix* dataset,
    sgpp::base::DataVector& source, sgpp::base::DataVector& result, const size_t start_index_grid,
    const size_t end_index_grid, const size_t start_index_data, const size_t end_index_data) {
  double* ptrLevel = level->getPointer();
  double* ptrIndex = index->getPointer();
  double* ptrSource = source.getPointer();
  double* ptrData = dataset->getPointer();
  double* ptrResult = result.getPointer();
  size_t sourceSize = source.getSize();
  size_t dims = dataset->getNrows();

#if defined(__SSE3__) && defined(__AVX2__) && !defined(__AVX512F__)

  for (size_t k = start_index_grid; k < end_index_grid; k++) {
    for (size_t i = start_index_data; i < end_index_data; i = i + 8) {
      // double support = ptrSource[i];
      __m256d support_0 = _mm256_load_pd(&ptrSource[i]);
      __m256d support_1 = _mm256_load_pd(&ptrSource[i + 4]);

      __m256d zero = _mm256_set1_pd(0.0);
      __m256d one = _mm256_set1_pd(1.0);
      __m256d two = _mm256_set1_pd(2.0);
      __m256d three = _mm256_set1_pd(3.0);
      __m256d four = _mm256_set1_pd(4.0);

      for (size_t d = 0; d < dims; d++) {
        __m256d levels = _mm256_broadcast_sd(&ptrLevel[(k * dims) + d]);
        __m256d indexs = _mm256_broadcast_sd(&ptrIndex[(k * dims) + d]);

        __m256d evals_0;
        __m256d evals_1;

        __m256d ptrData_0 = _mm256_load_pd(&ptrData[(d * sourceSize) + i]);
        __m256d ptrData_1 = _mm256_load_pd(&ptrData[(d * sourceSize) + i + 4]);

        __m256d x_0 = _mm256_fmsub_pd(ptrData_0, levels, indexs);
        __m256d x_1 = _mm256_fmsub_pd(ptrData_1, levels, indexs);

        x_0 = _mm256_add_pd(x_0, two);
        x_1 = _mm256_add_pd(x_1, two);

        __m256d n_0_0;
        __m256d n_0_1;
        __m256d n_0_0_temp;
        __m256d n_0_1_temp;

        __m256d result0_0;
        __m256d result0_1;

        // x>=0
        n_0_0 = _mm256_cmp_pd(x_0, zero, 13);
        n_0_1 = _mm256_cmp_pd(x_1, zero, 13);

        // x<1
        n_0_0_temp = _mm256_cmp_pd(x_0, one, 1);
        n_0_1_temp = _mm256_cmp_pd(x_1, one, 1);

        n_0_0 = _mm256_and_pd(n_0_0_temp, n_0_0);
        n_0_1 = _mm256_and_pd(n_0_1_temp, n_0_1);

        n_0_0 = _mm256_and_pd(n_0_0, one);
        n_0_1 = _mm256_and_pd(n_0_1, one);

        __m256d pow_x_2_0 = _mm256_mul_pd(x_0, x_0);
        __m256d pow_x_2_1 = _mm256_mul_pd(x_1, x_1);

        __m256d pow_x_3_0 = _mm256_mul_pd(pow_x_2_0, x_0);
        __m256d pow_x_3_1 = _mm256_mul_pd(pow_x_2_1, x_1);

        result0_0 = _mm256_mul_pd(pow_x_3_0, _mm256_set1_pd(1.0 / 6.0));
        result0_1 = _mm256_mul_pd(pow_x_3_1, _mm256_set1_pd(1.0 / 6.0));

        evals_0 = _mm256_mul_pd(result0_0, n_0_0);
        evals_1 = _mm256_mul_pd(result0_1, n_0_1);

        // n_1
        // x>=1
        n_0_0 = _mm256_cmp_pd(x_0, one, 13);
        n_0_1 = _mm256_cmp_pd(x_1, one, 13);

        // x<2
        n_0_0_temp = _mm256_cmp_pd(x_0, two, 1);
        n_0_1_temp = _mm256_cmp_pd(x_1, two, 1);
        n_0_0 = _mm256_and_pd(n_0_0_temp, n_0_0);
        n_0_1 = _mm256_and_pd(n_0_1_temp, n_0_1);

        n_0_0 = _mm256_and_pd(n_0_0, one);
        n_0_1 = _mm256_and_pd(n_0_1, one);

        result0_0 = _mm256_fmadd_pd(_mm256_set1_pd(-0.5), x_0, two);
        result0_1 = _mm256_fmadd_pd(_mm256_set1_pd(-0.5), x_1, two);
        result0_0 = _mm256_fmadd_pd(result0_0, x_0, _mm256_set1_pd(-2.0));
        result0_1 = _mm256_fmadd_pd(result0_1, x_1, _mm256_set1_pd(-2.0));
        result0_0 = _mm256_fmadd_pd(result0_0, x_0, _mm256_set1_pd(2.0 / 3.0));
        result0_1 = _mm256_fmadd_pd(result0_1, x_1, _mm256_set1_pd(2.0 / 3.0));

        evals_0 = _mm256_fmadd_pd(result0_0, n_0_0, evals_0);
        evals_1 = _mm256_fmadd_pd(result0_1, n_0_1, evals_1);

        // n2
        // x>=2
        n_0_0 = _mm256_cmp_pd(x_0, two, 13);
        n_0_1 = _mm256_cmp_pd(x_1, two, 13);

        // x<3
        n_0_0_temp = _mm256_cmp_pd(x_0, three, 1);
        n_0_1_temp = _mm256_cmp_pd(x_1, three, 1);
        n_0_0 = _mm256_and_pd(n_0_0_temp, n_0_0);
        n_0_1 = _mm256_and_pd(n_0_1_temp, n_0_1);

        n_0_0 = _mm256_and_pd(n_0_0, one);
        n_0_1 = _mm256_and_pd(n_0_1, one);

        result0_0 = _mm256_fmadd_pd(_mm256_set1_pd(0.5), x_0, _mm256_set1_pd(-4.0));
        result0_1 = _mm256_fmadd_pd(_mm256_set1_pd(0.5), x_1, _mm256_set1_pd(-4.0));
        result0_0 = _mm256_fmadd_pd(result0_0, x_0, _mm256_set1_pd(10.0));
        result0_1 = _mm256_fmadd_pd(result0_1, x_1, _mm256_set1_pd(10.0));
        result0_0 = _mm256_fmadd_pd(result0_0, x_0, _mm256_set1_pd(-22.0 / 3.0));
        result0_1 = _mm256_fmadd_pd(result0_1, x_1, _mm256_set1_pd(-22.0 / 3.0));

        evals_0 = _mm256_fmadd_pd(result0_0, n_0_0, evals_0);
        evals_1 = _mm256_fmadd_pd(result0_1, n_0_1, evals_1);

        // n3
        // x>=3
        n_0_0 = _mm256_cmp_pd(x_0, three, 13);
        n_0_1 = _mm256_cmp_pd(x_1, three, 13);

        // x<4
        n_0_0_temp = _mm256_cmp_pd(x_0, four, 1);
        n_0_1_temp = _mm256_cmp_pd(x_1, four, 1);
        n_0_0 = _mm256_and_pd(n_0_0_temp, n_0_0);
        n_0_1 = _mm256_and_pd(n_0_1_temp, n_0_1);

        n_0_0 = _mm256_and_pd(n_0_0, one);
        n_0_1 = _mm256_and_pd(n_0_1, one);

        result0_0 = _mm256_fmadd_pd(_mm256_set1_pd(-1.0 / 6.0), x_0, two);
        result0_1 = _mm256_fmadd_pd(_mm256_set1_pd(-1.0 / 6.0), x_1, two);
        result0_0 = _mm256_fmadd_pd(result0_0, x_0, _mm256_set1_pd(-8.0));
        result0_1 = _mm256_fmadd_pd(result0_1, x_1, _mm256_set1_pd(-8.0));
        result0_0 = _mm256_fmadd_pd(result0_0, x_0, _mm256_set1_pd(32.0 / 3.0));
        result0_1 = _mm256_fmadd_pd(result0_1, x_1, _mm256_set1_pd(32.0 / 3.0));

        evals_0 = _mm256_fmadd_pd(result0_0, n_0_0, evals_0);
        evals_1 = _mm256_fmadd_pd(result0_1, n_0_1, evals_1);

        support_0 = _mm256_mul_pd(support_0, evals_0);
        support_1 = _mm256_mul_pd(support_1, evals_1);
      }
      //
      // ptrResult[k] += support;

      support_0 = _mm256_add_pd(support_0, support_1);

      support_0 = _mm256_hadd_pd(support_0, support_0);
      __m256d temp = _mm256_permute2f128_pd(support_0, support_0, 0x81);
      support_0 = _mm256_add_pd(support_0, temp);

      double support_temp[4];

      _mm256_store_pd(&support_temp[0], support_0);

      ptrResult[k] += support_temp[0];
    }
  }

#elif defined(__SSE3__) && defined(__AVX__) && !defined(__AVX512F__)

  for (size_t k = start_index_grid; k < end_index_grid; k++) {
    for (size_t i = start_index_data; i < end_index_data; i = i + 4) {
      // double support = ptrSource[i];
      __m256d support = _mm256_load_pd(&ptrSource[i]);

      __m256d zero = _mm256_set1_pd(0.0);
      __m256d one = _mm256_set1_pd(1.0);
      __m256d two = _mm256_set1_pd(2.0);
      __m256d three = _mm256_set1_pd(3.0);
      __m256d four = _mm256_set1_pd(4.0);

      for (size_t d = 0; d < dims; d++) {
        __m256d levels = _mm256_broadcast_sd(&ptrLevel[(k * dims) + d]);
        __m256d indexs = _mm256_broadcast_sd(&ptrIndex[(k * dims) + d]);

        __m256d ptrDatas = _mm256_load_pd(&ptrData[(d * sourceSize) + i]);

        __m256d x = _mm256_sub_pd(_mm256_mul_pd(ptrDatas, levels), indexs);
        x = _mm256_add_pd(x, two);

        __m256d n;
        __m256d result;
        __m256d eval;
        __m256d pow_x_2 = _mm256_mul_pd(x, x);
        __m256d pow_x_3 = _mm256_mul_pd(pow_x_2, x);

        // x>=0
        n = _mm256_cmp_pd(x, zero, 13);
        // x<1
        n = _mm256_and_pd(_mm256_cmp_pd(x, one, 1), n);
        n = _mm256_and_pd(n, one);

        result = _mm256_mul_pd(n, _mm256_mul_pd(pow_x_3, _mm256_set1_pd(1.0 / 6.0)));
        eval = _mm256_mul_pd(result, n);

        // x>=1
        n = _mm256_cmp_pd(x, one, 13);
        // x<2
        n = _mm256_and_pd(_mm256_cmp_pd(x, two, 1), n);
        n = _mm256_and_pd(n, one);

        result = _mm256_mul_pd(_mm256_set1_pd(-0.5), pow_x_3);
        result = _mm256_add_pd(result, _mm256_mul_pd(two, pow_x_2));  // 10
        result = _mm256_sub_pd(result, _mm256_mul_pd(x, _mm256_set1_pd(2.0)));
        result = _mm256_add_pd(result, _mm256_set1_pd(2.0 / 3.0));
        result = _mm256_mul_pd(result, n);
        eval = _mm256_add_pd(result, eval);

        // x>=2
        n = _mm256_cmp_pd(x, two, 13);
        // x<3
        n = _mm256_and_pd(_mm256_cmp_pd(x, three, 1), n);
        n = _mm256_and_pd(n, one);

        result = _mm256_mul_pd(_mm256_set1_pd(0.5), pow_x_3);
        result = _mm256_sub_pd(result, _mm256_mul_pd(four, pow_x_2));
        result = _mm256_add_pd(result, _mm256_mul_pd(_mm256_set1_pd(10.0), x));
        result = _mm256_sub_pd(result, _mm256_set1_pd(22.0 / 3.0));
        result = _mm256_mul_pd(n, result);
        eval = _mm256_add_pd(result, eval);

        // x>=3
        n = _mm256_cmp_pd(x, three, 13);
        // x<4
        n = _mm256_and_pd(_mm256_cmp_pd(x, four, 1), n);
        n = _mm256_and_pd(n, one);
        result = _mm256_mul_pd(_mm256_set1_pd(-1.0 / 6.0), pow_x_3);
        result = _mm256_add_pd(result, _mm256_mul_pd(two, pow_x_2));
        result = _mm256_sub_pd(result, _mm256_mul_pd(_mm256_set1_pd(8.0), x));
        result = _mm256_add_pd(result, _mm256_set1_pd(32.0 / 3.0));
        result = _mm256_mul_pd(result, n);
        eval = _mm256_add_pd(result, eval);

        support = _mm256_mul_pd(support, eval);
      }

      support = _mm256_hadd_pd(support, support);
      __m256d temp = _mm256_permute2f128_pd(support, support, 0x81);
      support = _mm256_add_pd(support, temp);

      double support_temp[4];

      _mm256_store_pd(&support_temp[0], support);

      ptrResult[k] += support_temp[0];
    }
  }

#else
  // For-Schleife über Grid
  for (size_t k = start_index_grid; k < end_index_grid; k++) {
    for (size_t i = start_index_data; i < end_index_data; i++) {
      double support = ptrSource[i];

      for (size_t d = 0; d < dims; d++) {
        double level = ptrLevel[(k * dims) + d];
        double index = ptrIndex[(k * dims) + d];
        double x = ptrData[(d * sourceSize) + i] * level - index + 2.0;

        double eval = 0.0;

        if ((x < 0.0) || (x >= 4.0)) {
          eval = 0.0;
        } else if (x < 1.0) {
          eval = 1.0 / 6.0 * x * x * x;
        } else if (x < 2.0) {
          eval = -0.5 * x * x * x + 2.0 * x * x - 2.0 * x + 2.0 / 3.0;
        } else if (x < 3.0) {
          eval = 0.5 * x * x * x - 4.0 * x * x + 10.0 * x - 22.0 / 3.0;
        } else {
          eval = -1.0 / 6.0 * x * x * x + 2.0 * x * x - 8.0 * x + 32.0 / 3.0;
        }
        support *= eval;
      }
      ptrResult[k] += support;
    }
  }

#endif
}
}  // namespace datadriven
}  // namespace sgpp

