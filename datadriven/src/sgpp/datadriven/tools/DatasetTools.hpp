/*
 * DatasetTools.hpp
 *
 *  Created on: Jan 20, 2016
 *      Author: pfandedd
 */

#include <sgpp/globaldef.hpp>

#include <sgpp/base/datatypes/DataMatrix.hpp>
#include <sgpp/base/datatypes/DataVector.hpp>

namespace SGPP {
namespace datadriven {

class DatasetTools {
 public:
  static void splitset(base::DataMatrix& dataset, base::DataVector& datasetValues,
                       size_t kFold,
                       std::vector<base::DataMatrix>& trainingSets,
                       std::vector<base::DataVector>& trainingSetValues,
                       std::vector<base::DataMatrix>& testSets,
                       std::vector<base::DataVector>& testSetValues,
                       bool verbose = false);

};

}
}
