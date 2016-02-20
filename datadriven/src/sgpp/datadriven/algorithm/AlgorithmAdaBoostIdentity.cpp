// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <sgpp/datadriven/algorithm/AlgorithmAdaBoostIdentity.hpp>
#include <sgpp/datadriven/algorithm/DMWeightMatrix.hpp>
#include <sgpp/base/operation/BaseOpFactory.hpp>

#include <sgpp/globaldef.hpp>


namespace SGPP {
namespace datadriven {

AlgorithmAdaBoostIdentity::AlgorithmAdaBoostIdentity(SGPP::base::Grid&
    SparseGrid, size_t gridType, SGPP::base::HashGenerator::level_t gridLevel,
    SGPP::base::DataMatrix& trainData, SGPP::base::DataVector& trainDataClass,
    size_t NUM, float_t lambda, size_t IMAX, float_t eps, size_t IMAX_final,
    float_t eps_final, float_t firstLabel, float_t secondLabel, float_t threshold,
    float_t maxLambda, float_t minLambda, size_t searchNum, bool refine,
    size_t refineMode, size_t refineNum, size_t numberOfAda, float_t percentOfAda,
    size_t mode) : AlgorithmAdaBoostBase(SparseGrid, gridType, gridLevel, trainData,
          trainDataClass, NUM, lambda, IMAX, eps, IMAX_final, eps_final, firstLabel,
          secondLabel, threshold, maxLambda, minLambda, searchNum, refine, refineMode,
          refineNum, numberOfAda, percentOfAda, mode) {
}

AlgorithmAdaBoostIdentity::~AlgorithmAdaBoostIdentity() {
}

void AlgorithmAdaBoostIdentity::alphaSolver(float_t& lambda,
    SGPP::base::DataVector& weight, SGPP::base::DataVector& alpha, bool final) {
  std::unique_ptr<SGPP::base::OperationMatrix> C(
      SGPP::op_factory::createOperationIdentity(*this->grid));
  SGPP::datadriven::DMWeightMatrix WMatrix(*this->grid, *this->data, *C, lambda,
      weight);
  SGPP::base::DataVector rhs(alpha.getSize());
  WMatrix.generateb(*this->classes, rhs);

  if (final) {
    SGPP::solver::ConjugateGradients myCG(this->imax_final, this->epsilon_final);
    myCG.solve(WMatrix, alpha, rhs, false, false, -1.0);
  } else {
    SGPP::solver::ConjugateGradients myCG(this->imax, this->epsilon);
    myCG.solve(WMatrix, alpha, rhs, false, false, -1.0);
  }
}

}  // namespace datadriven
}  // namespace SGPP

