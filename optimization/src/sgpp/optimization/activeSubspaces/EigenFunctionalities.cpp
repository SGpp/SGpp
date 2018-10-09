// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

//#ifdef USE_EIGEN

#include <sgpp/optimization/activeSubspaces/EigenFunctionalities.hpp>

namespace sgpp {
namespace optimization {

Eigen::VectorXd DataVectorToEigen(sgpp::base::DataVector v) {
  return Eigen::Map<Eigen::VectorXd, Eigen::Unaligned>(v.data(), v.size());
}

sgpp::base::DataVector EigenToDataVector(Eigen::VectorXd e) {
  sgpp::base::DataVector v;
  v.resize(e.size());
  Eigen::VectorXd::Map(&v[0], e.size()) = e;
  return v;
}

}  // namespace optimization
}  // namespace sgpp

//#endif /* USE_EIGEN */
