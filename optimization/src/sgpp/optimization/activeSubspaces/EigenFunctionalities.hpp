// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#pragma once
//#ifdef USE_EIGEN

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

}  // namespace optimization
}  // namespace sgpp

//#endif /* USE_EIGEN */
