// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <sgpp/datadriven/datamining/modules/fitting/FitterConfigurationClassification.hpp>
#include <sgpp/datadriven/algorithm/DBMatDatabase.hpp>

namespace sgpp {
namespace datadriven {

FitterConfiguration* FitterConfigurationClassification::clone() const {
  return new FitterConfigurationClassification(*this);
}

void FitterConfigurationClassification::setupDefaults() {
  FitterConfigurationDensityEstimation::setupDefaults();
}
} /* namespace datadriven */
} /* namespace sgpp */
