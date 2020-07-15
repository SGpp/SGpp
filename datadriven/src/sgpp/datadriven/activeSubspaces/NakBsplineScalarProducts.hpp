// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#pragma once
// #ifdef USE_EIGEN

#include <map>
#include <sgpp/base/grid/Grid.hpp>
#include <sgpp/base/grid/type/ModNakBsplineGrid.hpp>
#include <sgpp/base/operation/hash/common/basis/NakBsplineBasis.hpp>
#include <sgpp/base/operation/hash/common/basis/NakBsplineBoundaryBasis.hpp>
#include <sgpp/base/operation/hash/common/basis/NakBsplineExtendedBasis.hpp>
#include <sgpp/base/operation/hash/common/basis/NakBsplineModifiedBasis.hpp>
#include <sgpp/base/tools/Distribution.hpp>
#include <sgpp/base/tools/DistributionUniform.hpp>
#include <sgpp/datadriven/activeSubspaces/GaussQuadrature.hpp>
#include <tuple>
#include <utility>

namespace sgpp {
namespace datadriven {

/**
 * Calculates and stores the scalar products of not a knot B-spline functions
 * Attention: Currently all scalar product calculations are saved to the same innerProducts entity.
 * When using different probability density functions, they get mixed up!
 * Initialize a NakBsplineScalarProducts for each pdf.
 */
class NakBsplineScalarProducts {
 public:
  /**
   * Constructor
   *
   * @param gridType1          	type of the first basis
   * @param gridType2			type of the second basis
   * @param degree1				degree of the firts basis
   * @param degree2 			degree of the second basis
   * @param quadOrder			order for the quadrature
   */
  NakBsplineScalarProducts(sgpp::base::GridType gridType1, sgpp::base::GridType gridType2,
                           size_t degree1, size_t degree2, size_t quadOrder)
      : gridType1(gridType1),
        gridType2(gridType2),
        degree1(degree1),
        degree2(degree2),
        quadOrder(quadOrder) {
    basis1 = initializeBasis(gridType1, degree1);
    basis2 = initializeBasis(gridType2, degree2);
    base::GaussLegendreQuadRule1D gauss;
    gauss.getLevelPointsAndWeightsNormalized(quadOrder, coordinates, weights);
  }

  /**
   * initializes a nak B spline basis according to gidType and degree
   * @param gridType	type of the basis
   * @param degree		degree of the basis
   * @return 			the nak Bsplien basis
   */
  std::shared_ptr<sgpp::base::SBasis> initializeBasis(sgpp::base::GridType gridType, size_t degree);

  /**
   * calculates the one dimensional integral \int f*g dx where f and g are B-spline basis
   * functions or first derivatives of B-spline basis functions
   *
   * @param level1 	level of the first B-spline
   * @param index1 	index of the first B-spline
   * @param dx1 	evaluate B-spline if False, evaluate d/dx B-spline if True
   * @param level2 	level of the second B-spline
   * @param index2 	index of the second B-spline
   * @param dx2 	evaluate B-spline if False, evaluate d/dx B-spline if True
   *
   * @return  integral (derivative of) first basis function * (derivative of) second basis
   * function
   */
  double basisScalarProduct(unsigned int level1, unsigned int index1, bool dx1, unsigned int level2,
                            unsigned int index2, bool dx2);

  /**
   * calculates the one dimensional integral \int f*g \rho dx where f and g are B-spline basis
   * functions or first derivatives of B-spline basis functions and \rho is a probability
   * DensityFunction
   *
   * @param level1 	level of the first B-spline
   * @param index1 	index of the first B-spline
   * @param dx1 	evaluate B-spline if False, evaluate d/dx B-spline if True
   * @param level2 	level of the second B-spline
   * @param index2 	index of the second B-spline
   * @param dx2 	evaluate B-spline if False, evaluate d/dx B-spline if True
   * @param pdf		probability density function
   *
   * @return  integral (derivative of) first basis function * (derivative of) second basis
   * function * probability density function
   */
  double weightedBasisScalarProduct(unsigned int level1, unsigned int index1, bool dx1,
                                    unsigned int level2, unsigned int index2, bool dx2,
                                    std::shared_ptr<sgpp::base::Distribution> pdf);

  /**
   * calculates the scalar product of two sparse grid nak B-spline interpolants given through their
   * grids and coefficients
   *
   * @param grid1	grid of the first interpolant
   * @param coeff1	coefficients of the first intepolant
   * @param grid2 	grid of the second interpolant
   * @param coeff2	coefficients of the second interpolant
   *
   * @return the scalar product \int \sum coeff1_i b_i(x) \sum coeff2_jb_j(x) dx
   *
   */
  double calculateScalarProduct(std::shared_ptr<sgpp::base::Grid> grid1,
                                sgpp::base::DataVector coeff1,
                                std::shared_ptr<sgpp::base::Grid> grid2,
                                sgpp::base::DataVector coeff2);

  /**
   * calculates the scalar product of two sparse grid nak B-spline interpolants given through their
   * grids and coefficients w.r.t a probability density function \rho
   *
   * @param grid1	grid of the first interpolant
   * @param coeff1	coefficients of the first intepolant
   * @param grid2 	grid of the second interpolant
   * @param coeff2	coefficients of the second interpolant
   * @param pdf		probability density function
   *
   * @return the scalar product \int \sum coeff1_i b_i(x) \sum coeff2_jb_j(x) \rho dx
   *
   */
  double calculateWeightedScalarProduct(std::shared_ptr<sgpp::base::Grid> grid1,
                                        sgpp::base::DataVector coeff1,
                                        std::shared_ptr<sgpp::base::Grid> grid2,
                                        sgpp::base::DataVector coeff2,
                                        std::shared_ptr<sgpp::base::Distribution> pdf);

 private:
  // type of first and second basis
  sgpp::base::GridType gridType1;
  sgpp::base::GridType gridType2;
  // degrees of first and second basis
  size_t degree1;
  size_t degree2;
  // order for the quadrature
  size_t quadOrder;
  // quadrature coordinates
  sgpp::base::DataVector coordinates;
  // quadrature weights
  sgpp::base::DataVector weights;
  // instances of first and second basis
  std::shared_ptr<sgpp::base::SBasis> basis1;
  std::shared_ptr<sgpp::base::SBasis> basis2;
  // tuple used as hash to store scalar products in innerProducts
  typedef std::tuple<size_t, size_t, bool, size_t, size_t, bool> asMatrixHashType;
  // hash storage for scalar products. Holds all calculated scalar products s.t. they do not have to
  // calculated again if the same combination of indices, levels and dx is queried
  std::map<asMatrixHashType, double> innerProducts;

  /**
   * used to get the support segments of a not a knot B-spline basis functions.
   *
   * @param level	level of the B-spline basis function
   * @param index	index of the B-spline basis function
   * @param degree	degree of the B-spline basis function
   *
   * @return the indices of the segments of the not a knot B-spline basis functions support
   */
  sgpp::base::DataVector nakBSplineSupport(size_t level, size_t index, size_t degree);
};

}  // namespace datadriven
}  // namespace sgpp

// #endif /* USE_EIGEN */
