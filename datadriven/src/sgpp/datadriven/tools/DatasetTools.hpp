// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <sgpp/globaldef.hpp>

#include <sgpp/base/datatypes/DataMatrix.hpp>
#include <sgpp/base/datatypes/DataVector.hpp>

#include <vector>

namespace sgpp {
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

}  // namespace datadriven
}  // namespace sgpp

