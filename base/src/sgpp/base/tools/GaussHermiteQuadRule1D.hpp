// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#ifndef GAUSSHERMITEQUADRULE1D_HPP_
#define GAUSSHERMITEQUADRULE1D_HPP_

#include "QuadRule1D.hpp"
#include <sgpp/globaldef.hpp>

namespace SGPP {
namespace base {

class GaussHermiteQuadRule1D : public QuadRule1D {
public:
    GaussHermiteQuadRule1D();
    virtual ~GaussHermiteQuadRule1D();

    /**
     * load gauss quadrature points for standard normal weight function. The points
     * and the weights are generated with numpy.polynomial.hermite.hermgauss,
     * the coordinates are scaled by sqrt(2), the weights are normalized to 1.
     */

    /**
     * the coordinates are scaled by sqrt(2) and then normalized with respect
     * to a given mean and standard deviation. The weights are normalized
     * to 1.
     *
     * @param level
     * @param coordinates
     * @param weights
     * @param mean
     * @param stdd
     */
    void getLevelPointsAndWeightsNormalized(size_t level,
            base::DataVector& coordinates, base::DataVector& weights,
            float_t mean = 0.0f, float_t stdd = 1.0f);

};

} /* namespace base */
} /* namespace SGPP */

#endif /* GAUSSHERMITEQUADRULE1D_HPP_ */
