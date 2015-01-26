/* ****************************************************************************
* Copyright (C) 2012 Technische Universitaet Muenchen                         *
* This file is part of the SG++ project. For conditions of distribution and   *
* use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
**************************************************************************** */
// @author Christoph Kowitz (kowitz@in.tum.de)

#include "combigrid/combigrid/AdaptiveSerialCombiGridVariableCoefficients.hpp"
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
