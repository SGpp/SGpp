// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#pragma once
#ifdef USE_EIGEN
#include <eigen3/Eigen/Dense>

#include <sgpp/base/datatypes/DataMatrix.hpp>
#include <sgpp/base/datatypes/DataVector.hpp>
#include <sgpp/base/grid/type/ModNakBsplineGrid.hpp>
#include <sgpp/base/grid/type/NakBsplineBoundaryGrid.hpp>
#include <sgpp/base/grid/type/NakBsplineExtendedGrid.hpp>
#include <sgpp/base/grid/type/NakBsplineGrid.hpp>
#include <sgpp/base/operation/hash/common/basis/NakBsplineBasis.hpp>
#include <sgpp/base/operation/hash/common/basis/NakBsplineBoundaryBasis.hpp>
#include <sgpp/base/operation/hash/common/basis/NakBsplineModifiedBasis.hpp>

namespace sgpp {
namespace datadriven {
/**
 * converts a SG++ DataVector to Eigen vector
 *
 * @param v SG++ DataVector
 * @return Eigen library vector containing the elements of v
 */
Eigen::VectorXd DataVectorToEigen(sgpp::base::DataVector v);

/**
 * converts an Eigen vector to a SG++ DataVector
 *
 * @param v Eigen library vector
 * @return  SG++ DataVector containing the elements of v
 */
sgpp::base::DataVector EigenToDataVector(Eigen::VectorXd e);

/**
 * converts an Eigen matrix to a SG++ DataMatrix
 *
 * @param m Eigen library matrix
 * @return  SG++ DataMatrix containing the elements of m
 */
sgpp::base::DataMatrix EigenToDataMatrix(Eigen::MatrixXd m);

/**
 * converts a SG++ DataMatrix to an Eigen::MatrixXd
 *
 * @param m SG++ DataMatrix containing the elements of m
 * @return  Eigen library matrix
 */
Eigen::MatrixXd DataMatrixToEigen(sgpp::base::DataMatrix d);

/**
 * calculates the regression for given data by  approximating the least squares optimal
 * coefficients with Tikhonov regularization
 *
 * @param grid				grid
 * @param degree			basis function degree
 * @param evaluationPoints  set of points in the original space of the objective function
 * @param functionValues	the objective function evaluated at evaluationPoints
 * @param mse				reference to return eman sqaure error
 * @param errorPerBasis		reference to return error per basis function
 * @param lambda			regularization parameter
 */
sgpp::base::DataVector EigenRegression(std::shared_ptr<sgpp::base::Grid> grid, size_t degree,
                                       Eigen::MatrixXd evaluationPoints,
                                       sgpp::base::DataVector functionValues, double& mse,
                                       sgpp::base::DataVector& errorPerBasis, double lambda = 1e-6);

}  // namespace datadriven
}  // namespace sgpp

#endif /* USE_EIGEN */
