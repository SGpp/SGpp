// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <sgpp/combigrid/combigrid/AdaptiveSerialCombiGridVariableCoefficients.hpp>
/* namespace combigrid */

void combigrid::AdaptiveSerialCombiGridVariableCoefficients::changeCoefficients(
  std::vector<double> newCoef) {
  //  COMBIGRID_OUT("Setting new coefficients in combischeme.");

  combischeme_->setCoef(newCoef);
  combikernel_->setCoef(newCoef);
}

void combigrid::AdaptiveSerialCombiGridVariableCoefficients::changeCoefficients(
  int i, double newCoef) {
  //  COMBIGRID_OUT("Setting new coefficient in combischeme.");

  combischeme_->setCoef(i, newCoef);
  combikernel_->setCoef(i, newCoef);
}