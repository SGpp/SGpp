// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <sgpp/datadriven/datamining/modules/hpo/harmonica/ConfigurationBit.hpp>
#include <sgpp/base/operation/hash/common/basis/LinearBasis.hpp>

namespace sgpp {
namespace datadriven {

void ConfigurationBit::addConstraint(ConfigurationRestriction *constraint) {
  constraints.push_back(constraint);
}

void ConfigurationBit::removeLastConstraint() {
  constraints.pop_back();
}

void ConfigurationBit::reset() {
  value = 0;
}

void ConfigurationBit::setValue(int input) {
  value = input;
  for (auto &constraint : constraints) {
    constraint->reduceOpenBits();
  }
}

int ConfigurationBit::getValue() {
  return value;
}

std::string ConfigurationBit::getName() {
  return name;
}
}  // namespace datadriven
}  // namespace sgpp
