// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#ifdef __AVX2__

#ifndef OPERATIONDENSITYMULTIPLICATIONAVX_H
#define OPERATIONDENSITYMULTIPLICATIONAVX_H

#include <x86intrin.h>
#include <malloc.h>
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
#define AVX_MULT_MAX_DIMENSIONS 32
#define AVX_MULT_INNER_BLOCKSIZE 128
    class OperationDensityMultiplicationAVX {
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

  const int blocksize = 128;

 public:
  OperationDensityMultiplicationAVX(base::Grid& grid, double lambda) {
    result = nullptr;
    // Store grid into int array and add values until it is divisible through 128
    // (required for opencl comparison)
    actual_gridsize = grid.getSize();
    dimensions = grid.getDimension();
    if (dimensions > AVX_MULT_MAX_DIMENSIONS) {
      std::string msg = std::string("Dimensions of grid too large (>") +
                    std::to_string(AVX_MULT_MAX_DIMENSIONS) +
                    std::string(").\nFixable by changing AVX_MULT_MAX_DIMENSIONS to a ") +
                    std::string("larger value (class OperationDensityMultiplicationAVX)");
      throw std::logic_error(msg);
    }
    used_gridsize =
        actual_gridsize + AVX_MULT_INNER_BLOCKSIZE - actual_gridsize % AVX_MULT_INNER_BLOCKSIZE;
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

    for (size_t i = 0; i < used_gridsize; ++i) {
      for (size_t dim = 0; dim < dimensions; dim++) {
        hs_inverse[i * dimensions + dim] = (1 << gridpoints[i*2*dimensions + 2*dim + 1]);
        hs[i * dimensions + dim] =
            (static_cast<double> (1.0 /(1 << gridpoints[i*2*dimensions + 2*dim + 1])));
        positions[i * dimensions + dim]
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
    if (this->result == nullptr)
      this->result = result.getPointer();
    else
      throw std::logic_error("Result pointer already bounded!");
    this->alpha = alpha.getPointer();
    start_partial_mult(0, alpha.getSize());
    finish_partial_mult(this->result, 0, result.getSize());
  }
  /// Execute a partial (startindex to startindex+chunksize) multiplication with the density matrix
  virtual void start_partial_mult(int start_id, int chunksize) {
    std::cerr << "Starting AVX mult with ..." << used_gridsize << std::endl;
    // size_t counter = 0;
    // size_t sicherheit = 0;
    #pragma omp parallel for
    for (size_t workitem = 0; workitem < used_gridsize; workitem+=AVX_MULT_INNER_BLOCKSIZE) {
      double last_positions[AVX_MULT_MAX_DIMENSIONS];
      __m256d last_integral[AVX_MULT_INNER_BLOCKSIZE / 4 * AVX_MULT_MAX_DIMENSIONS];
      __m256d workitem_positions[AVX_MULT_INNER_BLOCKSIZE / 4 * AVX_MULT_MAX_DIMENSIONS];
      __m256d workitem_hs[AVX_MULT_INNER_BLOCKSIZE / 4 * AVX_MULT_MAX_DIMENSIONS];
      __m256d workitem_hs_inverse[AVX_MULT_INNER_BLOCKSIZE / 4 * AVX_MULT_MAX_DIMENSIONS];
      __m256d zellenintegrals[AVX_MULT_INNER_BLOCKSIZE / 4];
      // load workitem positions
      for (size_t i = 0; i < AVX_MULT_INNER_BLOCKSIZE / 4; ++i) {
        for (size_t d = 0; d < dimensions; ++d) {
          workitem_positions[i * dimensions + d] =
              _mm256_set_pd(positions[(workitem + i * 4 + 3) * dimensions + d],
                            positions[(workitem + i * 4 + 2) * dimensions + d],
                            positions[(workitem + i * 4 + 1) * dimensions + d],
                            positions[(workitem + i * 4 + 0) * dimensions + d]);
          workitem_hs[i * dimensions + d] =
              _mm256_set_pd(hs[(workitem + i * 4 + 3) * dimensions + d],
                            hs[(workitem + i * 4 + 2) * dimensions + d],
                            hs[(workitem + i * 4 + 1) * dimensions + d],
                            hs[(workitem + i * 4 + 0) * dimensions + d]);
          workitem_hs_inverse[i * dimensions + d] =
              _mm256_set_pd(hs_inverse[(workitem + i * 4 + 3) * dimensions + d],
                            hs_inverse[(workitem + i * 4 + 2) * dimensions + d],
                            hs_inverse[(workitem + i * 4 + 1) * dimensions + d],
                            hs_inverse[(workitem + i * 4 + 0) * dimensions + d]);
        }
      }

      for (size_t point = 0; point < used_gridsize / AVX_MULT_INNER_BLOCKSIZE; point++) {
        for (size_t dim = 0; dim < dimensions; ++dim) {
          last_positions[dim] = -1.0;
        }
        for (size_t i = 0; i < AVX_MULT_INNER_BLOCKSIZE; ++i) {
          for (size_t local_item = 0; local_item < AVX_MULT_INNER_BLOCKSIZE / 4; local_item++) {
            zellenintegrals[local_item] = _mm256_set_pd(1.0, 1.0, 1.0, 1.0);
          }
          // __m256d zellenintegral = _mm256_set_pd(1.0, 1.0, 1.0, 1.0);
          // __m256d zellenintegral_unrolled = _mm256_set_pd(1.0, 1.0, 1.0, 1.0);
          for (size_t dim = 0; dim < dimensions; dim++) {
            // sicherheit++;
            if (positions[(point * AVX_MULT_INNER_BLOCKSIZE + i) * dimensions + dim] ==
                last_positions[dim]) {
              for (size_t local_item = 0; local_item < AVX_MULT_INNER_BLOCKSIZE / 4;
                   local_item+=2) {
                const size_t unrolled_local_item = local_item + 1;
                zellenintegrals[local_item] =
                    _mm256_mul_pd(zellenintegrals[local_item],
                                  last_integral[local_item * dimensions + dim]);
                zellenintegrals[unrolled_local_item] =
                    _mm256_mul_pd(zellenintegrals[unrolled_local_item],
                                  last_integral[unrolled_local_item * dimensions + dim]);
              }
              // counter++;
              continue;
            }
            // load gridpoint positions;
            const __m256d current_positions =
                _mm256_set1_pd(positions[(point * AVX_MULT_INNER_BLOCKSIZE + i)
                                         * dimensions + dim]);
            last_positions[dim] =
                positions[(point * AVX_MULT_INNER_BLOCKSIZE + i) * dimensions + dim];
            const __m256d current_hs =
                _mm256_set1_pd(hs[(point * AVX_MULT_INNER_BLOCKSIZE + i) * dimensions + dim]);
            const __m256d current_hs_inverse =
                _mm256_set1_pd(hs_inverse[(point * AVX_MULT_INNER_BLOCKSIZE + i)
                                          * dimensions + dim]);
            for (size_t local_item = 0; local_item < AVX_MULT_INNER_BLOCKSIZE / 4; local_item+=2) {
              const size_t unrolled_local_item = local_item + 1;
              __m256d tmp = _mm256_set_pd(1.0, 1.0, 1.0, 1.0);
              __m256d tmp_unrolled = _mm256_set_pd(1.0, 1.0, 1.0, 1.0);
              // Calculate distance
              __m256d distance =
                  _mm256_sub_pd(current_positions,
                                workitem_positions[local_item * dimensions + dim]);
              __m256d distance_unrolled =
                  _mm256_sub_pd(current_positions, workitem_positions[unrolled_local_item *
                                                                      dimensions + dim]);
              distance = _mm256_max_pd(_mm256_sub_pd(_mm256_setzero_pd(), distance), distance);
              distance_unrolled =
                  _mm256_max_pd(_mm256_sub_pd(_mm256_setzero_pd(),
                                              distance_unrolled), distance_unrolled);
              tmp = _mm256_fnmadd_pd(distance, workitem_hs_inverse[local_item * dimensions + dim],
                                     _mm256_set1_pd(1.0));
              tmp_unrolled =
                  _mm256_fnmadd_pd(distance_unrolled,
                                   workitem_hs_inverse[unrolled_local_item * dimensions + dim],
                                   _mm256_set1_pd(1.0));
              tmp = _mm256_mul_pd(tmp, current_hs);
              tmp_unrolled = _mm256_mul_pd(tmp_unrolled, current_hs);
              tmp = _mm256_max_pd(tmp, _mm256_setzero_pd());
              tmp_unrolled = _mm256_max_pd(tmp_unrolled, _mm256_setzero_pd());
              __m256d tmp_alternative = _mm256_fnmadd_pd(distance, current_hs_inverse,
                                                         _mm256_set1_pd(1.0));
              __m256d tmp_alternative_unrolled =
                  _mm256_fnmadd_pd(distance_unrolled, current_hs_inverse,
                                   _mm256_set1_pd(1.0));
              tmp_alternative = _mm256_max_pd(tmp_alternative, _mm256_setzero_pd());
              tmp_alternative_unrolled =
                  _mm256_max_pd(tmp_alternative_unrolled, _mm256_setzero_pd());
              tmp = _mm256_fmadd_pd(tmp_alternative,
                                    workitem_hs[local_item * dimensions + dim], tmp);
              tmp_unrolled = _mm256_fmadd_pd(tmp_alternative_unrolled,
                                             workitem_hs[unrolled_local_item * dimensions + dim],
                                             tmp_unrolled);
              const __m256d bitmask = _mm256_cmp_pd(workitem_hs[local_item * dimensions + dim],
                                                    current_hs, _CMP_EQ_OQ);
              const __m256d bitmask_unrolled =
                  _mm256_cmp_pd(workitem_hs[unrolled_local_item * dimensions + dim],
                                current_hs, _CMP_EQ_OQ);
              tmp_alternative = _mm256_mul_pd(_mm256_set1_pd(1.0/3.0), tmp);
              tmp_alternative_unrolled = _mm256_mul_pd(_mm256_set1_pd(1.0/3.0), tmp_unrolled);
              const __m256d bitmask2 = _mm256_cmp_pd(workitem_hs[local_item * dimensions + dim],
                                                     current_hs, _CMP_NEQ_OQ);
              const __m256d bitmask2_unrolled =
                  _mm256_cmp_pd(workitem_hs[unrolled_local_item * dimensions + dim],
                                current_hs, _CMP_NEQ_OQ);
              tmp = _mm256_and_pd(tmp, bitmask2);
              tmp_unrolled = _mm256_and_pd(tmp_unrolled, bitmask2_unrolled);
              tmp = _mm256_add_pd(tmp, _mm256_and_pd(tmp_alternative, bitmask));
              tmp_unrolled = _mm256_add_pd(tmp_unrolled,
                                           _mm256_and_pd(tmp_alternative_unrolled,
                                                         bitmask_unrolled));
              last_integral[local_item * dimensions + dim] = tmp;
              last_integral[unrolled_local_item * dimensions + dim] = tmp_unrolled;
              zellenintegrals[local_item] = _mm256_mul_pd(zellenintegrals[local_item], tmp);
              zellenintegrals[unrolled_local_item] =
                  _mm256_mul_pd(zellenintegrals[unrolled_local_item], tmp_unrolled);
            }
          }
          for (size_t local_item = 0; local_item < AVX_MULT_INNER_BLOCKSIZE / 4; local_item+=2) {
            const size_t unrolled_local_item = local_item + 1;
            __m256d currentworkitem = _mm256_load_pd(&result[workitem + local_item * 4]);
            __m256d currentworkitem_unrolled =
                _mm256_load_pd(&result[workitem + unrolled_local_item * 4]);
            currentworkitem =
                _mm256_fmadd_pd(zellenintegrals[local_item],
                                _mm256_set1_pd(alpha[point * AVX_MULT_INNER_BLOCKSIZE + i]),
                                currentworkitem);
            currentworkitem_unrolled =
                _mm256_fmadd_pd(zellenintegrals[unrolled_local_item],
                                _mm256_set1_pd(alpha[point * AVX_MULT_INNER_BLOCKSIZE + i]),
                                currentworkitem_unrolled);

            _mm256_store_pd(&result[workitem + local_item * 4], currentworkitem);
            _mm256_store_pd(&result[workitem + unrolled_local_item * 4], currentworkitem_unrolled);
          }
        }
      }
    }

    // size_t maximum = sicherheit;
    // std::cout << "Skipped " << counter << " out of " << maximum << " ["
    //           << double(counter) / double(maximum) * 100.0 << "%]" << std::endl;
    // std::cout << "Calculated " << maximum - counter << std::endl;
  }
  void print_avx_register(__m256d reg) {
    double tmp_result[4];
    _mm256_storeu_pd(tmp_result, reg);
    std::cout << tmp_result[0] << " " << tmp_result[1] << " " << tmp_result[2] << " "
              << tmp_result[3] << std::endl;
  }
  /// Just a dummy function
  virtual void initialize_alpha(double *alpha) {
  }
  /// Just a dummy function
  virtual void initialize_dataset(base::DataMatrix &dataset) {
  }
  /// Just a dummy function
  virtual void finish_partial_mult(double *result, int start_id, int chunksize) {
  }
  /// Generates the right hand side vector for the density equation
  virtual void generateb(base::DataMatrix &dataset, sgpp::base::DataVector &b,
                         size_t start_id = 0,  size_t chunksize = 0) {
  }
  virtual void start_rhs_generation(size_t start_id,  size_t chunksize) {
  }
  virtual void finalize_rhs_generation(sgpp::base::DataVector &b,
                                       size_t start_id,  size_t chunksize) {
  }
};

}  // namespace DensityAVX
}  // namespace datadriven
}  // namespace sgpp
#endif /* OPERATIONDENSITYMULTIPLICATIONAVX_H */

#endif
