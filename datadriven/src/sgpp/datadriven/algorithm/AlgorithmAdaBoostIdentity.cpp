// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <sgpp/datadriven/algorithm/AlgorithmAdaBoostIdentity.hpp>
#include <sgpp/datadriven/algorithm/DMWeightMatrix.hpp>
#include <sgpp/base/operation/BaseOpFactory.hpp>

#include <sgpp/globaldef.hpp>

namespace sgpp {
namespace datadriven {

AlgorithmAdaBoostIdentity::AlgorithmAdaBoostIdentity(sgpp::base::Grid& SparseGrid,
                                                     base::GridType gridType,
                                                     sgpp::base::level_t gridLevel,
                                                     sgpp::base::DataMatrix& trainData,
                                                     sgpp::base::DataVector& trainDataClass,
                                                     size_t NUM, double lambda, size_t IMAX,
                                                     double eps, size_t IMAX_final,
                                                     double eps_final, double firstLabel,
                                                     double secondLabel, double threshold,
                                                     double maxLambda, double minLambda,
                                                     size_t searchNum, bool refine,
                                                     size_t refineMode, size_t refineNum,
                                                     size_t numberOfAda, double percentOfAda,
                                                     size_t mode)
    : AlgorithmAdaBoostBase(SparseGrid, gridType, gridLevel, trainData, trainDataClass, NUM, lambda,
                            IMAX, eps, IMAX_final, eps_final, firstLabel, secondLabel, threshold,
                            maxLambda, minLambda, searchNum, refine, refineMode, refineNum,
                            numberOfAda, percentOfAda, mode) {
}

AlgorithmAdaBoostIdentity::~AlgorithmAdaBoostIdentity() {
}

void AlgorithmAdaBoostIdentity::alphaSolver(double& lambda, sgpp::base::DataVector& weight,
                                            sgpp::base::DataVector& alpha, bool final) {
  std::unique_ptr < sgpp::base::OperationMatrix
      > C(sgpp::op_factory::createOperationIdentity(*this->grid));
  sgpp::datadriven::DMWeightMatrix WMatrix(*this->grid, *this->data, *C, lambda, weight);
  sgpp::base::DataVector rhs(alpha.getSize());
  WMatrix.generateb(*this->classes, rhs);

  if (final) {
    sgpp::solver::ConjugateGradients myCG(this->imax_final, this->epsilon_final);
    myCG.solve(WMatrix, alpha, rhs, false, false, -1.0);
  } else {
    sgpp::solver::ConjugateGradients myCG(this->imax, this->epsilon);
    myCG.solve(WMatrix, alpha, rhs, false, false, -1.0);
  }
}

}  // namespace datadriven
}  // namespace sgpp
