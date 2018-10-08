// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#pragma once

#include <sgpp/base/tools/GaussLegendreQuadRule1D.hpp>

#include <functional>

namespace sgpp {
namespace optimization {

/**
 * calculates \int_a^b f(x) dx
 */
double quad(std::function<double(double)> f, double a, double b, size_t quadOrder);

}  // namespace optimization
}  // namespace sgpp
