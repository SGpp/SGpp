// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#pragma once
//#ifdef USE_EIGEN

#include <sgpp/base/grid/type/NakBsplineGrid.hpp>
#include <sgpp/base/operation/hash/common/basis/NakBsplineBasis.hpp>
#include <sgpp/base/operation/hash/common/basis/NakBsplineBoundaryBasis.hpp>
#include <sgpp/base/operation/hash/common/basis/NakBsplineModifiedBasis.hpp>
#include <sgpp/optimization/activeSubspaces/ASMatrix.hpp>
#include <sgpp/optimization/activeSubspaces/GaussQuadrature.hpp>
#include <sgpp/optimization/function/scalar/InterpolantScalarFunction.hpp>
#include <sgpp/optimization/function/scalar/InterpolantScalarFunctionGradient.hpp>
#include <sgpp/optimization/sle/solver/Armadillo.hpp>
#include <sgpp/optimization/sle/system/HierarchisationSLE.hpp>
#include <sgpp/optimization/tools/RandomNumberGenerator.hpp>

namespace sgpp {
namespace optimization {

class ASMatrixNakBspline : public ASMatrix {
 public:
  ASMatrixNakBspline(WrapperScalarFunction objectiveFunc, sgpp::base::GridType gridType,
                     size_t degree)
      : ASMatrix(objectiveFunc), gridType(gridType), degree(degree) {
    if (gridType == sgpp::base::GridType::NakBspline) {
      grid = std::make_shared<sgpp::base::NakBsplineGrid>(numDim, degree);
    } else if (gridType == sgpp::base::GridType::NakBsplineBoundary) {
      std::cout << "numDim " << numDim << std::endl;
      grid = std::make_shared<sgpp::base::NakBsplineBoundaryGrid>(numDim, degree);
    } else if (gridType == sgpp::base::GridType::NakBsplineModified) {
      grid = std::make_shared<sgpp::base::NakBsplineModifiedGrid>(numDim, degree);
    } else {
      throw sgpp::base::generation_exception("ASMatrixNakBspline: gridType not supported.");
    }
  }
  void buildRegularInterpolant(size_t level);
  void createMatrix(size_t numPoints);
  void createMatrixMonteCarlo(size_t numPoints);
  void createMatrixGauss();

  // auxiliary routines

  /**
   * calculates the entry C_{i,j} of the matrix,
   * C_{i,j} = int df/dxi df/dxj dx
   *
   * @param i row
   * @param j column
   * @return matrix entry C_{i,j}
   */
  double matrixEntryGauss(size_t i, size_t j);

  /**
   * calculates \int d/dx_i b_k(x) d/dx_j b_l(x) dx
   * This is done by evaluating one dimensional integrals of b_{k_d} b_{l_d},
   * d/dx_i b_{k_i} b_{l_i}, b_{k_j} d/dx_j b_{l_j}, d/dx_i b_{k_i} d/dx_i b_{l_i}
   *
   *@param i index of matrix row
   *@param j index of matrix column
   *@param k index of first basis function
   *@param l index of second basis function
   *
   * @return integral \int d/dx_i b_k(x) d/dx_j b_l(x) dx
   */
  double scalarProductDxbiDxbj(size_t i, size_t j, size_t k, size_t l);

  /**
   * calculates the one diensional integral \int f*g dx where f and g are B-spline basis functions
   * or first derivatives of B-spline basis functions
   *
   * @param level1 level of the first B-spline
   * @param index1 index of the first B-spline
   * @param dx1 evaluate B-spline if False, evaluate d/dx B-spline if True
   * @param level2 level of the second B-spline
   * @param index2 index of the second B-spline
   * @param dx2 evaluate B-spline if False, evaluate d/dx B-spline if True
   *
   * @return  integral (derivative of) first basis function * (derivative of) second basis function
   */
  double univariateScalarProduct(size_t level1, size_t index1, bool dx1, size_t level2,
                                 size_t index2, bool dx2);

  sgpp::base::DataVector nakBSplineSupport(size_t level, size_t index);

 private:
  sgpp::base::GridType gridType;
  size_t degree;
  sgpp::base::DataVector coefficients;
  std::shared_ptr<sgpp::base::Grid> grid;
  bool interpolantFlag = 0;
};

}  // namespace optimization
}  // namespace sgpp

//#endif /* USE_EIGEN */
