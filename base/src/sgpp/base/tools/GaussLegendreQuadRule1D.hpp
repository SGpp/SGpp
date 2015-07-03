// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#ifndef GAUSSLEGENDREQUADRULE1D_HPP_
#define GAUSSLEGENDREQUADRULE1D_HPP_

#include "QuadRule1D.hpp"
#include <sgpp/globaldef.hpp>

namespace SGPP {
namespace base {

class GaussLegendreQuadRule1D : public QuadRule1D {
public:
    /**
     * load gauss quadrature points for uniform weight function. The points
     * and the weights are generated with numpy.polynomial.legendre.leggauss.
     * the weights are additionally normalized to 1.
     */
    GaussLegendreQuadRule1D();
    virtual ~GaussLegendreQuadRule1D();

    /**
     * the coordinates are normalized to [0, 1].
     *
     * @param level level of quadrature, is equal to the number of quadrature points
     * @param coordinates returns the x-coordinates in [0, 1]
     * @param weights returns the corresponding weights (scaled by 0.5)
     */
    void getLevelPointsAndWeightsNormalized(size_t level,
            base::DataVector& coordinates, base::DataVector& weights);
};

} /* namespace base */
} /* namespace SGPP */

#endif /* GAUSSLEGENDREQUADRULE1D_HPP_ */
