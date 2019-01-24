// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#pragma once

#include <sgpp/base/tools/GaussLegendreQuadRule1D.hpp>

#include <functional>

namespace sgpp {
namespace datadriven {

/**
 * calculates \int_a^b f(x) dx
 */
double gaussQuad(std::function<double(double)> f, double a, double b, size_t quadOrder);
// the initialization of GaussLegendreQuadRule1D is slow. Do not use this frequently
// If needed often initialize GaussLegendreQuadRule1D in the beginning and hand coordinates and
// weights through and copy the summation from here to the corresponding place

}  // namespace datadriven
}  // namespace sgpp
