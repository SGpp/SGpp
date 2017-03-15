#ifndef MORTONORDER_HPP
#define MORTONORDER_HPP

#include <sgpp/base/datatypes/DataVector.hpp>
#include <sgpp/base/datatypes/DataMatrix.hpp>

#include <sgpp/globaldef.hpp>

#include <stdint.h>
#include <vector>

namespace sgpp {
namespace datadriven {

/// Class for re-arranging Datasets along a Morton order curve
class MortonOrder {
 public:
  /// Generates the identic permutation list
  MortonOrder (size_t size);
  /// Generates the permutation list according to the given dataset. The computation is done on GPU
  MortonOrder (sgpp::base::DataMatrix& matrix);

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
 protected:
  std::vector<size_t> permutation;
};

}  // datadriven
}  // sgpp

#endif // MORTONORDER_HPP
