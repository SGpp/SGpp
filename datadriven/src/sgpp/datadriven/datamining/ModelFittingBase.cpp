/*
 * ModelFittingBase.cpp
 *
 *  Created on: Feb 8, 2016
 *      Author: franzefn
 */

#include "../datamining/ModelFittingBase.hpp"

using namespace SGPP::base;

namespace SGPP {
namespace datadriven {

ModelFittingBase::ModelFittingBase() : grid(nullptr), alpha(0) {
}

ModelFittingBase::~ModelFittingBase() {
}

void ModelFittingBase::evaluate(DataMatrix& samples, DataVector& result) {
  uint32_t numSamples = samples.getNrows();
  uint32_t numDims = samples.getNcols();

  DataVector sample(numDims);
  for (size_t i = 0; i < numSamples; i++) {
    samples.getRow(i, sample);
    result[i] = evaluate(sample);
  }
}

std::shared_ptr<SGPP::base::Grid> ModelFittingBase::getGrid() {
  return grid;
}

std::shared_ptr<SGPP::base::DataVector> ModelFittingBase::getSurpluses() {
  return alpha;
}

} /* namespace datadriven */
} /* namespace SGPP */
