// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <sgpp/datadriven/datamining/modules/hpo/parameters/HyperParameter.hpp>

#include <vector>
#include <iostream>
#include <string>


namespace sgpp {
namespace datadriven {

void HyperParameter::makeConfigBits(std::vector<ConfigurationBit *> &configBits) {
  bits.reserve(nBits);
  for (size_t i = 0; i < nBits; i++) {
    std::string bitName = name + std::to_string(i);
    bits.emplace_back(bitName);
    configBits.push_back(&bits[i]);
  }
}
} /* namespace datadriven */
} /* namespace sgpp */
