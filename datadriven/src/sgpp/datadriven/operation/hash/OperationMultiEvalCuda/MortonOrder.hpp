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
  /// Generates the permutation list according to the given dataset. The computation is done on the GPU
  MortonOrder (sgpp::base::DataMatrix& matrix);
  
  /// Re-arrange a DataMatrix object
  void orderDataMatrix(sgpp::base::DataMatrix& matrix) const;
  /// Re-arrange a DataVector object
  void orderDataVector(sgpp::base::DataVector& data) const;
  /// Restores the original order of a DataMatrix object
  void restoreDataMatrix(sgpp::base::DataMatrix& matrix) const;
  /// Restores the original order of a DataVector object
  void restoreDataVector(sgpp::base::DataVector& data) const;
 protected:
  std::vector<size_t> permutation;
};

}  // datadriven
}  // sgpp

#endif // MORTONORDER_HPP
