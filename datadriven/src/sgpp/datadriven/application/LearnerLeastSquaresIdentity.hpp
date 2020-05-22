// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#pragma once

#include <sgpp/datadriven/application/LearnerBase.hpp>
#include <sgpp/base/operation/hash/OperationMultipleEval.hpp>
#include <sgpp/datadriven/operation/hash/DatadrivenOperationCommon.hpp>

#include <sgpp/globaldef.hpp>

#include <utility>
#include <string>
#include <vector>

namespace sgpp {
namespace datadriven {

/**
 * This class implements standard sparse grid regression
 * with an Identity matrix as regularization operator.
 *
 * Furthermore this Learner provides support for several
 * vectorization approaches covering GPUs, CPUs and coprocessors.
 */
class LearnerLeastSquaresIdentity : public sgpp::datadriven::LearnerBase {
 private:
  std::vector<std::pair<size_t, double> > ExecTimeOnStep;

  // sgpp::base::OperationMultipleEval* kernel = nullptr;

  sgpp::datadriven::OperationMultipleEvalConfiguration implementationConfiguration;

 protected:
  std::unique_ptr<sgpp::datadriven::DMSystemMatrixBase> createDMSystem(
      sgpp::base::DataMatrix& trainDataset, double lambda) override;

  void postProcessing(const sgpp::base::DataMatrix& trainDataset,
                      const sgpp::solver::SLESolverType& solver,
                      const size_t numNeededIterations) override;

 public:
  using LearnerBase::predict;

  /**
   * Constructor
   *
   * @param isRegression set to true if a regression task should be executed
   * @param isVerbose set to true in order to allow console output
   */
  explicit LearnerLeastSquaresIdentity(const bool isRegression, const bool isVerbose = true);

  /**
   * Destructor
   */
  ~LearnerLeastSquaresIdentity() override;

  void predict(sgpp::base::DataMatrix& testDataset,
               sgpp::base::DataVector& classesComputed) override;

  void multTranspose(sgpp::base::DataMatrix& dataset, sgpp::base::DataVector& multiplier,
                     sgpp::base::DataVector& result) override;

  double testRegular(const sgpp::base::RegularGridConfiguration& gridConfig,
                     sgpp::base::DataMatrix& testDataset);

  std::vector<std::pair<size_t, double> > getRefinementExecTimes();

  void setImplementation(
      sgpp::datadriven::OperationMultipleEvalConfiguration operationConfiguration) {
    this->implementationConfiguration = operationConfiguration;
  }
};

}  // namespace datadriven
}  // namespace sgpp
