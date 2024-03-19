// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#pragma once

#include <cmath>
#include <functional>

namespace sgpp {
namespace combigrid {
double squaredSummation(double a, double b) { return std::sqrt(a * a + b * b); }
std::function<double(double, double)> squaredSummation_f = squaredSummation;
}  // namespace combigrid
}  // namespace sgpp
