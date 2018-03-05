// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#ifndef MORTONORDER_HPP
#define MORTONORDER_HPP

#include <sgpp/base/datatypes/DataMatrix.hpp>
#include <sgpp/base/datatypes/DataVector.hpp>

#include <sgpp/globaldef.hpp>

#include <stdint.h>
#include <vector>

///@cond DOXY_IGNORE // NOLINT()
namespace sgpp {
namespace datadriven {
namespace OpMultiEvalCudaDetail {

/// Class for re-arranging Datasets along a Morton order curve
class MortonOrder {
 public:
  /// Generates the identic permutation list
  explicit MortonOrder(size_t size);
  /// Generates the permutation list according to the given dataset. The computation is done on GPU
  explicit MortonOrder(sgpp::base::DataMatrix& matrix);

  /// Re-arrange a DataMatrix object inplace
  void orderDataMatrix(sgpp::base::DataMatrix& matrix) const;
  /// Re-arrange a DataMatrix object with other destination
  void orderDataMatrix(sgpp::base::DataMatrix& matrix, double* dest) const;
  /// Re-arrange a DataVector object inplace
  void orderDataVector(sgpp::base::DataVector& data) const;
  /// Re-arrange a DataVector object with other destination
  void orderDataVector(sgpp::base::DataVector& data, double* dest) const;
  /// Restores the original order of a DataMatrix object inplace
  void restoreDataMatrix(sgpp::base::DataMatrix& matrix) const;
  /// Restores the original order of a DataMatrix object with other source
  void restoreDataMatrix(sgpp::base::DataMatrix& matrix, double* src) const;
  /// Restores the original order of a DataVector object inplace
  void restoreDataVector(sgpp::base::DataVector& data) const;
  /// Restores the original order of a DataVector object with other source
  void restoreDataVector(sgpp::base::DataVector& data, double* src) const;

  /// Check if permutation is identity
  bool isIdentity() const;

 protected:
  std::vector<size_t> permutation;
};

}  // namespace OpMultiEvalCudaDetail
}  // namespace datadriven
}  // namespace sgpp
///@endcond // NOLINT()

#endif  // MORTONORDER_HPP
