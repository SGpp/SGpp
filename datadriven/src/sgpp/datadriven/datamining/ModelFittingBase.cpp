/*
 * ModelFittingBase.cpp
 *
 *  Created on: Feb 8, 2016
 *      Author: franzefn
 */

#include "ModelFittingBase.hpp"

using namespace SGPP::base;

namespace SGPP {
  namespace datadriven {

    ModelFittingBase::ModelFittingBase(DataMiningConfiguration config) : config(config), grid(nullptr), alpha(0) {
    }

    ModelFittingBase::~ModelFittingBase() {
    }

    SGPP::float_t ModelFittingBase::evaluate(DataVector& sample) {
      return 0.0;
    }

    void ModelFittingBase::evaluate(DataMatrix& samples, DataVector& result) {
    }

    std::shared_ptr<SGPP::base::Grid> ModelFittingBase::getGrid() {
      return grid;
    }

    std::shared_ptr<SGPP::base::DataVector> ModelFittingBase::getSurpluses() {
      return alpha;
    }

  } /* namespace datadriven */
} /* namespace SGPP */
