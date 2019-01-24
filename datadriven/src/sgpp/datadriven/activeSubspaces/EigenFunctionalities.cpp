// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

//#ifdef USE_EIGEN

#include "../../../../../datadriven/src/sgpp/datadriven/activeSubspaces/EigenFunctionalities.hpp"

namespace sgpp {
namespace datadriven {

Eigen::VectorXd DataVectorToEigen(sgpp::base::DataVector v) {
  return Eigen::Map<Eigen::VectorXd, Eigen::Unaligned>(v.data(), v.size());
}

sgpp::base::DataVector EigenToDataVector(Eigen::VectorXd e) {
  sgpp::base::DataVector v;
  v.resize(e.size());
  Eigen::VectorXd::Map(&v[0], e.size()) = e;
  return v;
}

// This is inefficient!
sgpp::base::DataMatrix EigenToDataMatrix(Eigen::MatrixXd m) {
  sgpp::base::DataMatrix d(m.rows(), m.cols());
  for (size_t i = 0; i < static_cast<size_t>(m.rows()); i++) {
    for (size_t j = 0; j < static_cast<size_t>(m.cols()); j++) {
      d.set(i, j, m(i, j));
    }
  }
  return d;
}

// This too is inefficient!
Eigen::MatrixXd DataMatrixToEigen(sgpp::base::DataMatrix d) {
  Eigen::MatrixXd m(d.getNrows(), d.getNcols());
  for (size_t i = 0; i < d.getNrows(); i++) {
    for (size_t j = 0; j < d.getNcols(); j++) {
      m(i, j) = d.get(i, j);
    }
  }
  return m;
}

}  // namespace datadriven
}  // namespace sgpp

//#endif /* USE_EIGEN */
