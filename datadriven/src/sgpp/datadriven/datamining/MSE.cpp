// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include "MSE.hpp"

#include <sgpp/globaldef.hpp>

namespace SGPP {
namespace datadriven {

MSE::MSE() {
  // TODO(valeriy) Auto-generated constructor stub
}

MSE::~MSE() {
  // TODO(valeriy) Auto-generated destructor stub
}

double MSE::operator()(DataVector& predictedValues, DataVector& trueValues) {
  DataVector tmp(predictedValues);
  tmp.sub(trueValues);
  double error = tmp.l2Norm();
  return error * error / static_cast<double>(tmp.getSize());
}

} /* namespace datadriven */
} /* namespace SGPP */
