// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#ifndef OPERATIONLTWODOTPRODUCTVECTORIZEDLINEARBOUNDARY_HPP
#define OPERATIONLTWODOTPRODUCTVECTORIZEDLINEARBOUNDARY_HPP

#include <sgpp/base/operation/hash/OperationMatrix.hpp>
#include <sgpp/base/datatypes/DataMatrix.hpp>
#include <sgpp/base/grid/Grid.hpp>
#include <sgpp/base/tools/SGppStopwatch.hpp>
#include <sgpp/parallel/tools/TypesParallel.hpp>

#include <vector>

#if defined(__SSE4_2__) || defined(__AVX__)
#include <immintrin.h>
#endif
#if defined(__FMA4__)
#include <x86intrin.h>
#endif

#ifdef __USEAVX128__
#undef __AVX__
#endif

#include <sgpp/globaldef.hpp>

namespace sgpp {
namespace parallel {

/**
 * Implementation for linear functions of LTwoDotProduct Operation, linear grids with boundaries
 *
 */
class OperationLTwoDotProductVectorizedLinearBoundary : public sgpp::base::OperationMatrix {
 private:
  sgpp::base::GridStorage* storage;
  sgpp::base::DataMatrix* level_;
  sgpp::base::DataMatrix* level_int_;
  sgpp::base::DataMatrix* index_;
  sgpp::base::DataVector* lcl_q_;
  sgpp::base::DataVector* constants_;
  sgpp::base::DataVector* alpha_padded_;

  sgpp::base::DataVector* result_boundary_filtered_;

  sgpp::base::DataMatrix* level_boundary_filtered_;
  sgpp::base::DataMatrix* level_int_boundary_filtered_;
  sgpp::base::DataMatrix* index_boundary_filtered_;

  sgpp::base::DataVector** gradient_temp;
  sgpp::base::DataVector** l2dot_temp;

#if defined(STORE_PDE_MATRIX_BOUNDARY)
  sgpp::base::DataMatrix* operation_result_matrix_;
  bool operation_result_generated_;
#endif
  std::vector<std::size_t> i_boundary_filtered;

  int process_count;
  int process_index;

  std::vector<int> all_i_start;
  std::vector<int> all_i_size;

  std::vector<int> send_start;
  std::vector<int> send_size;

  std::vector<int> recv_start;
  std::vector<int> recv_size;

  void init_constants();
  void init_grid_storage();

  double l2dot_dirichlet(size_t i, size_t j, size_t dim);

  void mult_dirichlet(sgpp::base::DataVector& alpha, sgpp::base::DataVector& result);
  double all_time;
  double all_iterations;
  sgpp::base::SGppStopwatch stopWatch;

 public:
  /**
   * Construtor of OperationLTwoProductLinear
   *
   * @param storage Pointer to the grid's gridstorage obejct
   */
  explicit OperationLTwoDotProductVectorizedLinearBoundary(sgpp::base::GridStorage* storage);

  /**
   * Destructor
   */
  virtual ~OperationLTwoDotProductVectorizedLinearBoundary();

  virtual void mult(sgpp::base::DataVector& alpha, sgpp::base::DataVector& result);
  virtual void reset();
};
}  // namespace parallel
}  // namespace sgpp

#endif /* OPERATIONLAPLACEVECTORIZEDLINEARBOUNDARY_HPP */