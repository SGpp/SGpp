// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#pragma once

#include <sgpp/base/datatypes/DataMatrix.hpp>
#include <sgpp/base/datatypes/DataVector.hpp>
#include <sgpp/combigrid/definitions.hpp>
#include <sgpp/combigrid/operation/Configurations.hpp>

#include <algorithm>
#include <vector>

#include "../../../../../base/src/sgpp/base/operation/hash/common/basis/NakBsplineBasis.hpp"
#include "../../../../../base/src/sgpp/base/operation/hash/common/basis/NakBsplineModifiedBasis.hpp"

namespace sgpp {
namespace combigrid {

/**
 * evaluates a modified not a knot Bspline on an expUniformGrid, given by its degree, index and the
 *knot sequence it is defined on in x.
 *@param x       evaluation point
 *@param degree     B-spline degree
 *@param i       index of B-spline
 *@param points    points of the 1D grid
 *@return        value of non-uniform B-spline in x
 */
double expUniformModifiedNakBspline(double const& x, size_t const& degree, size_t i,
                                    std::vector<double> points);

}  // namespace combigrid
}  // namespace sgpp
