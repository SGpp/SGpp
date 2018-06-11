// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <sgpp/base/algorithm/AlgorithmMultipleEvaluation.hpp>
#include <sgpp/datadriven/datamining/modules/hpo/harmonica/ConfigurationRestriction.hpp>
#include <sgpp/datadriven/datamining/modules/hpo/harmonica/ConfigurationBit.hpp>
#include <sgpp/base/operation/hash/common/basis/LinearBasis.hpp>
#include <sgpp/base/datatypes/DataMatrix.hpp>

#include <sgpp/globaldef.hpp>
#include <list>

namespace sgpp {
namespace datadriven {

ConfigurationRestriction::ConfigurationRestriction(std::vector<ConfigurationBit *> &parameters,
                                                   int bias)
    : parameters(parameters), bias(bias) {
}

void ConfigurationRestriction::reduceOpenBits() {
  openBits--;
}

void ConfigurationRestriction::resolve() {
  size_t idx = 100;
  int tmp = bias;
  for (size_t i = 0; i < parameters.size(); ++i) {
    int v = parameters[i]->getValue();
    if (v == 0) {
      idx = i;
    } else {
      tmp = tmp * v;
    }
  }
  parameters[idx]->setValue(tmp);
}

int ConfigurationRestriction::getOpenBits() {
  return openBits;
}

bool ConfigurationRestriction::check() {
  int tmp = 1;
  for (auto &bit : parameters) {
    tmp = tmp * bit->getValue();
  }
  return tmp == bias;
}

void ConfigurationRestriction::reset() {
  openBits = static_cast<int>(parameters.size());
}

}  // namespace datadriven
}  // namespace sgpp
