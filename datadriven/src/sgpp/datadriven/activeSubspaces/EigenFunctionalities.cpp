// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

// #ifdef USE_EIGEN

#include <sgpp/datadriven/activeSubspaces/EigenFunctionalities.hpp>

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

// ToDo(rehmemk) This is inefficient!
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

sgpp::base::DataVector EigenRegression(std::shared_ptr<sgpp::base::Grid> grid, size_t degree,
                                       Eigen::MatrixXd evaluationPoints,
                                       sgpp::base::DataVector functionValues, double& mse,
                                       sgpp::base::DataVector& errorPerBasis, double lambda) {
  std::unique_ptr<sgpp::base::SBasis> basis;
  sgpp::base::GridType gridType = grid->getType();
  if (gridType == sgpp::base::GridType::NakBspline) {
    basis = std::make_unique<sgpp::base::SNakBsplineBase>(degree);
  } else if (gridType == sgpp::base::GridType::NakBsplineBoundary) {
    basis = std::make_unique<sgpp::base::SNakBsplineBoundaryBase>(degree);
  } else if (gridType == sgpp::base::GridType::NakBsplineModified) {
    basis = std::make_unique<sgpp::base::SNakBsplineModifiedBase>(degree);
  } else if (gridType == sgpp::base::GridType::NakBsplineExtended) {
    basis = std::make_unique<sgpp::base::SNakBsplineExtendedBase>(degree);
  } else {
    throw sgpp::base::generation_exception("EigenFunctionalities: gridType not supported.");
  }

  sgpp::base::GridStorage& gridStorage = grid->getStorage();
  // regressionMatrix(i,j) = b_j (y_i) = b_j (W1T * evaluationPoint_i)
  Eigen::MatrixXd regressionMatrix(evaluationPoints.cols(), gridStorage.getSize());
  for (unsigned int i = 0; i < evaluationPoints.cols(); i++) {
    Eigen::VectorXd point = evaluationPoints.col(i);

    for (size_t j = 0; j < gridStorage.getSize(); j++) {
      double basisEval = 1;
      for (size_t t = 0; t < gridStorage.getDimension(); t++) {
        double basisEval1D =
            basis->eval(gridStorage.getPointLevel(j, t), gridStorage.getPointIndex(j, t), point(t));
        if (basisEval1D == 0) {
          basisEval = 0;
          break;
        } else {
          basisEval *= basisEval1D;
        }
      }
      regressionMatrix(i, j) = basisEval;
    }
  }
  Eigen::VectorXd functionValues_Eigen = DataVectorToEigen(functionValues);
  // three different ways to solve least squares with Eigen
  // https://eigen.tuxfamily.org/dox/group__LeastSquares.html

  Eigen::setNbThreads(4);
  /* 1. SVD, slowest but most accurate*/
  //  Eigen::VectorXd alpha_Eigen = regressionMatrix.bdcSvd(Eigen::ComputeThinU |
  //  Eigen::ComputeThinV)
  //                                    .solve(functionValues_Eigen);
  /* 2. QR medium speed, medium accuracy*/
  //  Eigen::VectorXd alpha_Eigen =
  //  regressionMatrix.colPivHouseholderQr().solve(functionValues_Eigen);
  /* 3. normal equations, fastest but least accurate. Terrible for ill conditioned matrices*/
  //  Eigen::VectorXd alpha_Eigen = (regressionMatrix.transpose() * regressionMatrix)
  //                                    .ldlt()
  //                                    .solve(regressionMatrix.transpose() * functionValues_Eigen);
  /* 4. = 3. with Tikhonov regularization*/
  //  Eigen::MatrixXd identity =
  //      Eigen::MatrixXd::Identity(regressionMatrix.cols(), regressionMatrix.cols());
  //  double lambda = 1e-6;
  //  Eigen::VectorXd alpha_Eigen =
  //      (regressionMatrix.transpose() * regressionMatrix / gridStorage.getSize() + lambda *
  //      identity)
  //          .ldlt()
  //          .solve(regressionMatrix.transpose() * functionValues_Eigen / gridStorage.getSize());
  /* 5. = 1. with Tikhonov*/
  Eigen::MatrixXd identity =
      Eigen::MatrixXd::Identity(regressionMatrix.cols(), regressionMatrix.cols());
  Eigen::VectorXd alpha_Eigen =
      (regressionMatrix.transpose() * regressionMatrix / gridStorage.getSize() + lambda * identity)
          .bdcSvd(Eigen::ComputeThinU | Eigen::ComputeThinV)
          .solve(regressionMatrix.transpose() * functionValues_Eigen / gridStorage.getSize());

  sgpp::base::DataVector alpha = EigenToDataVector(alpha_Eigen);

  // calculate MSE and residual
  Eigen::VectorXd error = regressionMatrix * alpha_Eigen;
  for (unsigned int i = 0; i < error.size(); i++) {
    error(i) *= error(i);
  }
  mse = error.sum() / error.size();
  Eigen::VectorXd errorPerBasis_Eigen = regressionMatrix.transpose() * error;
  errorPerBasis = EigenToDataVector(errorPerBasis_Eigen);
  errorPerBasis.componentwise_mult(alpha);

  //  std::cout << "residual = " << (regressionMatrix * alpha_Eigen - functionValues_Eigen).norm()
  //            << "\n";

  return alpha;
}

}  // namespace datadriven
}  // namespace sgpp

// #endif /* USE_EIGEN */
