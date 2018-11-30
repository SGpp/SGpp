// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#pragma once
//#ifdef USE_EIGEN

#include <sgpp/base/datatypes/DataMatrix.hpp>
#include <sgpp/base/datatypes/DataVector.hpp>

#include <eigen3/Eigen/Dense>
#include <eigen3/Eigen/Eigenvalues>

namespace sgpp {
namespace optimization {
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
sgpp::base::DataMatrix DataMatrixToEigen(Eigen::MatrixXd m);

}  // namespace optimization
}  // namespace sgpp

//#endif /* USE_EIGEN */
