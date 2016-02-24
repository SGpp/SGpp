// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#ifndef LEARNERVECTORIZEDIDENTITY_HPP
#define LEARNERVECTORIZEDIDENTITY_HPP

#include <sgpp/datadriven/application/LearnerBase.hpp>

#include <sgpp/parallel/tools/TypesParallel.hpp>

#include <sgpp/globaldef.hpp>

#include <string>
namespace SGPP {

namespace parallel {

/**
 * This class implements standard sparse grid regression
 * with an Identity matrix as regularization operator.
 *
 * Furthermore this Learner provides support for several
 * vectorization approaches covering GPUs, CPUs and coprocessors.
 */
class LearnerVectorizedIdentity : public SGPP::datadriven::LearnerBase {
 protected:
  /// vectorization selector
  VectorizationType vecType_;

  MPIType mpiType_;

  virtual SGPP::datadriven::DMSystemMatrixBase* createDMSystem(SGPP::base::DataMatrix& trainDataset,
                                                               double lambda);

  virtual void postProcessing(const SGPP::base::DataMatrix& trainDataset,
                              const SGPP::solver::SLESolverType& solver,
                              const size_t numNeededIterations);

 public:
  /**
   * Constructor
   *
   * @param vecType selection of vectorization to employ
   * @param isRegression set to true if a regression task should be executed
   * @param isVerbose set to true in order to allow console output
   */
  LearnerVectorizedIdentity(const VectorizationType vecType, const bool isRegression,
                            const bool isVerbose = true);
  LearnerVectorizedIdentity(const VectorizationType vecType, const MPIType mpiType,
                            const bool isRegression, const bool isVerbose = true);

  /**
   * Constructor
   *
   * @param tGridFilename path to file that contains a serialized grid
   * @param tAlphaFilename path to file that contains the grid's coefficients
   * @param vecType selection of vectorization to employ
   * @param isRegression set to true if a regression task should be executed
   * @param isVerbose set to true in order to allow console output
   */
  LearnerVectorizedIdentity(const std::string tGridFilename, const std::string tAlphaFilename,
                            const VectorizationType vecType, const bool isRegression,
                            const bool isVerbose = true);

  /**
   * Destructor
   */
  virtual ~LearnerVectorizedIdentity();

  virtual SGPP::base::DataVector predict(SGPP::base::DataMatrix& testDataset);
};

}  // namespace parallel
}  // namespace SGPP

#endif /* LEARNERVECTORIZEDIDENTITY_HPP */
