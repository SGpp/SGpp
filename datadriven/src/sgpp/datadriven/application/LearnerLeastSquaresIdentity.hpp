// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#ifndef LEARNERLEASTSQUARESIDENTITY_HPP
#define LEARNERLEASTSQUARESIDENTITY_HPP

#include <sgpp/datadriven/application/LearnerBase.hpp>
#include <sgpp/base/operation/hash/OperationMultipleEval.hpp>
#include <sgpp/datadriven/operation/hash/simple/DatadrivenOperationCommon.hpp>

#include <sgpp/globaldef.hpp>


namespace SGPP {

namespace datadriven {

/**
 * This class implements standard sparse grid regression
 * with an Identity matrix as regularization operator.
 *
 * Furthermore this Learner provides support for several
 * vectorization approaches covering GPUs, CPUs and coprocessors.
 */
class LearnerLeastSquaresIdentity: public SGPP::datadriven::LearnerBase {
 private:
  std::vector<std::pair<size_t, float_t> > ExecTimeOnStep;

  SGPP::base::OperationMultipleEval* kernel = nullptr;

  SGPP::datadriven::OperationMultipleEvalConfiguration
  implementationConfiguration;
 protected:

  virtual SGPP::datadriven::DMSystemMatrixBase* createDMSystem(
    SGPP::base::DataMatrix& trainDataset, float_t lambda);

  virtual void postProcessing(const SGPP::base::DataMatrix& trainDataset,
                              const SGPP::solver::SLESolverType& solver,
                              const size_t numNeededIterations);

 public:
  /**
   * Constructor
   *
   * @param isRegression set to true if a regression task should be executed
   * @param isVerbose set to true in order to allow console output
   */
  LearnerLeastSquaresIdentity(const bool isRegression,
                              const bool isVerbose = true);

  /**
   * Constructor
   *
   * @param tGridFilename path to file that contains a serialized grid
   * @param tAlphaFilename path to file that contains the grid's coefficients
   * @param isRegression set to true if a regression task should be executed
   * @param isVerbose set to true in order to allow console output
   */
  LearnerLeastSquaresIdentity(const std::string tGridFilename,
                              const std::string tAlphaFilename,
                              const bool isRegression, const bool isVerbose = true);

  /**
   * Destructor
   */
  virtual ~LearnerLeastSquaresIdentity();

  virtual SGPP::base::DataVector predict(SGPP::base::DataMatrix& testDataset);

  float_t testRegular(const SGPP::base::RegularGridConfiguration& GridConfig,
                      SGPP::base::DataMatrix& testDataset);

  std::vector<std::pair<size_t, float_t> > getRefinementExecTimes();

  void setImplementation(SGPP::datadriven::OperationMultipleEvalConfiguration
                         operationConfiguration) {
    this->implementationConfiguration = operationConfiguration;
  }
};

}

}

#endif /* LEARNERLEASTSQUARESIDENTITY_HPP */
