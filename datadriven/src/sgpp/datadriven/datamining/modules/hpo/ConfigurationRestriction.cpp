// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <sgpp/base/algorithm/AlgorithmMultipleEvaluation.hpp>
#include <sgpp/datadriven/datamining/modules/hpo/ConfigurationRestriction.hpp>
#include <sgpp/base/operation/hash/common/basis/LinearBasis.hpp>
#include <sgpp/base/datatypes/DataMatrix.hpp>

#include <sgpp/globaldef.hpp>
#include <list>

namespace sgpp {
namespace datadriven {



int ConfigurationRestriction::getBias(){
  return bias;
}

std::list<ConfigurationBit> ConfigurationRestriction::getConfigBits(){ // std::list<ConfigurationBit>
  return parameters;
}

}  // namespace datadriven
}  // namespace sgpp
