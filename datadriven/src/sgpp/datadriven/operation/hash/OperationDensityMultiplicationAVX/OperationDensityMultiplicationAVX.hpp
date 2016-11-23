// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#ifndef OPERATIONDENSITYMULTIPLICATIONAVX_H
#define OPERATIONDENSITYMULTIPLICATIONAVX_H

#include <x86intrin.h>
#include <malloc.h>
#include <sgpp/datadriven/operation/hash/OperationDensityOCLMultiPlatform/OperationDensityOCL.hpp>
#include <sgpp/base/grid/Grid.hpp>
#include <sgpp/base/grid/GridStorage.hpp>
#include <sgpp/base/tools/SGppStopwatch.hpp>
#include <sgpp/base/exception/operation_exception.hpp>
#include <sgpp/globaldef.hpp>
#include <sgpp/base/operation/hash/OperationMatrix.hpp>

#include <chrono>
#include <vector>
#include <string>
#include <exception>

namespace sgpp {
namespace datadriven {
namespace DensityAVX {
class OperationDensityMultiplicationAVX : public DensityOCLMultiPlatform::OperationDensity {
 private:
  size_t actual_gridsize;
  size_t used_gridsize;
  size_t dimensions;

  double *result;
  int *gridpoints;
  double *positions;
  double *hs;
  double *hs_inverse;

  double *alpha;

 public:
  OperationDensityMultiplicationAVX(base::Grid& grid, double lambda) {
    result = NULL;
    // Store grid into int array and add values until it is divisible through 128
    // (required for opencl comparison)
    actual_gridsize = grid.getSize();
    used_gridsize = actual_gridsize + 128 - actual_gridsize % 128;
    dimensions = grid.getDimension();
    gridpoints = new int[used_gridsize * 2 * dimensions];
    sgpp::base::GridStorage& gridStorage = grid.getStorage();
    for (size_t i = 0; i < actual_gridsize; ++i) {
      sgpp::base::HashGridPoint &point = gridStorage.getPoint(i);
      for (size_t d = 0; d < dimensions; d++) {
        gridpoints[i * 2 * dimensions + d * 2] = (point.getIndex(d));
        gridpoints[i * 2 * dimensions + d * 2 + 1] = (point.getLevel(d));
      }
    }
    // padding
    for (size_t i = actual_gridsize; i < used_gridsize; ++i) {
      for (size_t d = 0; d < dimensions; d++) {
        gridpoints[i * 2 * dimensions + d * 2] = 0;
        gridpoints[i * 2 * dimensions + d * 2 + 1] = 0;
      }
    }

    // Store grid point positions and hs
    positions = new double[used_gridsize * dimensions];
    hs = new double[used_gridsize * dimensions];
    hs_inverse = new double[used_gridsize * dimensions];

    for (size_t dim = 0; dim < dimensions; dim++) {
      for (size_t i = 0; i < used_gridsize; ++i) {
        hs_inverse[dim * used_gridsize + i] = ((int)1 << gridpoints[i*2*dimensions + 2*dim + 1]);
        hs[dim * used_gridsize + i]
            = (static_cast<double>(1.0 / ((int)1 << gridpoints[i*2*dimensions + 2*dim + 1])));
        positions[dim * used_gridsize + i]
            = (static_cast<double>(gridpoints[i*2*dimensions + 2*dim]) /
               static_cast<double>(1 << gridpoints[i*2*dimensions + 2*dim + 1]));
      }
    }
    // padding
    // for (size_t i = actual_gridsize; i < used_gridsize; ++i) {
    //   for (size_t d = 0; d < dimensions; d++) {
    //     positions[d * used_gridsize + i] = 0.0;
    //     hs[d * used_gridsize + i] = 0.0;
    //     hs_inverse[d * used_gridsize + i] = 0.0;
    //   }
    // }
  }
  virtual ~OperationDensityMultiplicationAVX() {
    delete [] gridpoints;
    delete [] positions;
    delete [] hs;
    delete [] hs_inverse;
  }

  /// Execute one matrix-vector multiplication with the density matrix
  virtual void mult(base::DataVector& alpha, base::DataVector& result) {
    if (this->result == NULL)
      this->result = result.getPointer();
    else
      throw std::logic_error("Result pointer already bounded!");
    this->alpha = alpha.getPointer();
    start_partial_mult(alpha.getPointer(), 0, alpha.getSize());
    finish_partial_mult(this->result, 0, result.getSize());
  }
  /// Execute a partial (startindex to startindex+chunksize) multiplication with the density matrix
  virtual void start_partial_mult(double *alpha, int start_id, int chunksize) {
    std::cerr << "Starting AVX mult with ..." << used_gridsize << std::endl;
    // Copy into SIMD vectors
    const size_t blocksize = 128;
    #pragma omp parallel for
    for (size_t workitem = 0; workitem < actual_gridsize; workitem+=4) {
      __m256d zellenintegral[blocksize];
      __m256d *workitem_positions = new __m256d[dimensions];
      __m256d *workitem_hs = new __m256d[dimensions];
      __m256d *workitem_hs_inverse = new __m256d[dimensions];
      __m256d current_workitems;
      __m256d tmp;
      __m256d tmp_alternative;
      __m256d distance;
      __m256d current_positions;
      __m256d current_hs;
      __m256d current_hs_inverse;
      __m256d bitmask;
      __m256d bitmask2;
      current_workitems = _mm256_setzero_pd();
      tmp = _mm256_setzero_pd();
      // load workitem positions
      for (size_t d = 0; d < dimensions; ++d) {
        workitem_positions[d] =
            _mm256_load_pd(&positions[d * used_gridsize + workitem]);
        workitem_hs[d] =
            _mm256_load_pd(&hs[d * used_gridsize + workitem]);
        workitem_hs_inverse[d] =
            _mm256_load_pd(&hs_inverse[d * used_gridsize + workitem]);
      }

      for (size_t point = 0; point < used_gridsize / blocksize; point++) {
        for (size_t i = 0; i < blocksize; ++i) {
          zellenintegral[i] = _mm256_set_pd(1.0, 1.0, 1.0, 1.0);
        }
        for (size_t dim = 0; dim < dimensions; dim++) {
          for (size_t i = 0; i < blocksize; ++i) {
            // load gridpoint positions;
            current_positions = _mm256_set1_pd(positions[dim * used_gridsize +
                                                         point * blocksize + i]);
            current_hs = _mm256_set1_pd(hs[dim * used_gridsize + point * blocksize + i]);
            current_hs_inverse = _mm256_set1_pd(hs_inverse[dim * used_gridsize +
                                                          point * blocksize + i]);
            tmp = _mm256_set_pd(1.0, 1.0, 1.0, 1.0);
            // Calculate distance
            distance = _mm256_sub_pd(current_positions, workitem_positions[dim]);
            distance = _mm256_max_pd(_mm256_sub_pd(_mm256_setzero_pd(), distance), distance);
            tmp = _mm256_fnmadd_pd(distance, workitem_hs_inverse[dim],
                                   _mm256_set1_pd(1.0));
            tmp = _mm256_mul_pd(tmp, current_hs);
            tmp = _mm256_max_pd(tmp,_mm256_setzero_pd());
            tmp_alternative = _mm256_fnmadd_pd(distance, current_hs_inverse,
                                               _mm256_set1_pd(1.0));
            tmp_alternative = _mm256_max_pd(tmp_alternative,_mm256_setzero_pd());
            tmp = _mm256_fmadd_pd(tmp_alternative, workitem_hs[dim], tmp);
            bitmask = _mm256_cmp_pd(workitem_hs[dim], current_hs, _CMP_EQ_OQ);
            tmp_alternative = _mm256_mul_pd(_mm256_set1_pd(1.0/3.0), tmp);
            bitmask2 = _mm256_cmp_pd(workitem_hs[dim], current_hs, _CMP_NEQ_OQ);
            tmp = _mm256_and_pd(tmp, bitmask2);
            tmp = _mm256_add_pd(tmp, _mm256_and_pd(tmp_alternative, bitmask));
            zellenintegral[i] = _mm256_mul_pd(zellenintegral[i], tmp);
          }
        }
        for (size_t i = 0; i < blocksize; ++i) {
        current_workitems =
            _mm256_fmadd_pd(zellenintegral[i],
                            _mm256_set1_pd(alpha[point * blocksize + i]),
                            current_workitems);
        }
      }
      if (workitem + 4 < actual_gridsize)
        _mm256_storeu_pd(&result[workitem], current_workitems);

      delete [] workitem_positions;
      delete [] workitem_hs;
      delete [] workitem_hs_inverse;
    }
  }
  void print_avx_register(__m256d reg) {
    double tmp_result[4];
    _mm256_storeu_pd(tmp_result, reg);
    std::cout << tmp_result[0] << " " << tmp_result[1] << " " << tmp_result[2] << " "
              << tmp_result[3] << std::endl;
  }
  /// Just a dummy function
  virtual void finish_partial_mult(double *result, int start_id, int chunksize) {
  }
  /// Generates the right hand side vector for the density equation
  virtual void generateb(base::DataMatrix &dataset, sgpp::base::DataVector &b,
                         size_t start_id = 0,  size_t chunksize = 0) {
  }
  virtual void start_rhs_generation(base::DataMatrix &dataset,
                                    size_t start_id,  size_t chunksize) {
  }
  virtual void finalize_rhs_generation(sgpp::base::DataVector &b,
                                       size_t start_id,  size_t chunksize) {
  }
};

}  // namespace DensityAVX
}  // namespace datadriven
}  // namespace sgpp
#endif /* OPERATIONDENSITYMULTIPLICATIONAVX_H */
