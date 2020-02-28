// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <sgpp/datadriven/datamining/modules/hpo/harmonica/ConfigurationRestriction.hpp>

#include <sgpp/base/algorithm/AlgorithmMultipleEvaluation.hpp>
#include <sgpp/datadriven/datamining/modules/hpo/harmonica/ConfigurationBit.hpp>

#include <vector>
#include <string>

namespace sgpp {
namespace datadriven {

ConfigurationRestriction::ConfigurationRestriction(std::vector<ConfigurationBit *> &parameters,
                                                   int bias)
    : parameters(parameters), bias(bias) {
}

void ConfigurationRestriction::reduceOpenBits() {
  openBits--;
  if (openBits == 2) {
    findComplex();
  }
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

void ConfigurationRestriction::findComplex() {
  std::string otherbit;
  size_t idx = 100;
  int effectivebias = bias;
  for (size_t i = 0; i < parameters.size(); ++i) {
    if (parameters[i]->getValue() == 0) {
      if (otherbit.empty()) {
        otherbit = parameters[i]->getName();
      } else {
        idx = i;
      }
    } else {
      effectivebias *= parameters[i]->getValue();
    }
  }
  parameters[idx]->findComplexinner(otherbit, effectivebias);
}

void ConfigurationRestriction::resolveComplex(const std::string &id,
                                              const std::string &idtwo, int otherbias) {
  bool resolve = false;
  size_t idx = 100;
  for (size_t i = 0; i < parameters.size(); ++i) {
    if (parameters[i]->getName() == id) {
      resolve = true;
    } else if (parameters[i]->getName() != idtwo) {
      idx = i;
    }
  }
  if (resolve) {
    parameters[idx]->setValue(bias*otherbias);
  }
}

}  // namespace datadriven
}  // namespace sgpp
