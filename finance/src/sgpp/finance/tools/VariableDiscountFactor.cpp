// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <sgpp/finance/tools/VariableDiscountFactor.hpp>

#include <sgpp/globaldef.hpp>

#include <string>

namespace sgpp {
namespace finance {

VariableDiscountFactor::VariableDiscountFactor(sgpp::base::GridStorage* storage, int dim_r)
    : myBoundingBox(storage->getBoundingBox()), storage(storage), dim_r(dim_r) {}

VariableDiscountFactor::~VariableDiscountFactor() {}

void VariableDiscountFactor::getDiscountFactor(sgpp::base::DataVector& factor, double T) {
  double tmp;

  for (size_t i = 0; i < storage->getSize(); i++) {
    std::string coords = storage->getCoordinates((*storage)[i]).toString();
    std::stringstream coordsStream(coords);
    double dblFuncValues[2];

    for (size_t j = 0; j < 2; j++) {
      coordsStream >> tmp;
      dblFuncValues[j] = tmp;
    }

    // std::cout<<dblFuncValues[1]<<std::endl;
    // factor.set(i, exp((-1.0)*dblFuncValues[1]*T));
    factor.set(i, exp((-1.0) * dblFuncValues[this->dim_r] * T));
  }
}
}  // namespace finance
}  // namespace sgpp
