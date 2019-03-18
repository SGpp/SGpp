/* Copyright (C) 2008-today The SG++ project
 * This file is part of the SG++ project. For conditions of distribution and
 * use, please see the copyright notice provided with SG++ or at
 * sgpp.sparsegrids.org
 * ModelFittingBaseSingleGrid.cpp
 *
 *  Created on: May 22, 2018
 *      Author: dominik
 */

#include <sgpp/datadriven/datamining/modules/fitting/ModelFittingBaseSingleGrid.hpp>

#include <sgpp/base/exception/application_exception.hpp>

#include <string>

namespace sgpp {
namespace datadriven {

using base::DataVector;
using base::Grid;
using base::application_exception;

ModelFittingBaseSingleGrid::ModelFittingBaseSingleGrid()
    : ModelFittingBase(), grid{nullptr}, alpha{} {}

Grid& ModelFittingBaseSingleGrid::getGrid() {
  if (grid != nullptr) {
    return *grid;
  } else {
    throw application_exception("No grid was fitted yet");
  }
}

DataVector& ModelFittingBaseSingleGrid::getSurpluses() {
  return alpha;
}

std::string ModelFittingBaseSingleGrid::storeFitter() {
  std::string output;
  output = output + "Grid: \n" + getGrid().serialize() + "\n";
  output = output + "Alphas \n" + getSurpluses().toString() + "\n";
  return output;
}

} /* namespace datadriven */
} /* namespace sgpp */
