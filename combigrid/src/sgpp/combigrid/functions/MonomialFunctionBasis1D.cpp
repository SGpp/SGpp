// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <sgpp/combigrid/functions/MonomialFunctionBasis1D.hpp>
#include <sgpp/combigrid/utils/Utils.hpp>

namespace sgpp {
namespace combigrid {

MonomialFunctionBasis1D::~MonomialFunctionBasis1D() {}

double MonomialFunctionBasis1D::evaluate(size_t basisIndex, double xValue) {
  return pow(xValue, basisIndex);
}

} /* namespace combigrid */
} /* namespace sgpp */
