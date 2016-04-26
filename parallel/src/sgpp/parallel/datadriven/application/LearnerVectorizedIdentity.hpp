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
namespace sgpp {

namespace parallel {

/**
 * This class implements standard sparse grid regression
 * with an Identity matrix as regularization operator.
 *
 * Furthermore this Learner provides support for several
 * vectorization approaches covering GPUs, CPUs and coprocessors.
 */
class LearnerVectorizedIdentity : public sgpp::datadriven::LearnerBase {
 protected:
  /// vectorization selector
  VectorizationType vecType_;

  MPIType mpiType_;

  virtual sgpp::datadriven::DMSystemMatrixBase* createDMSystem(sgpp::base::DataMatrix& trainDataset,
                                                               double lambda);

  virtual void postProcessing(const sgpp::base::DataMatrix& trainDataset,
                              const sgpp::solver::SLESolverType& solver,
                              const size_t numNeededIterations);

 public:
  using LearnerBase::predict;

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

  virtual sgpp::base::DataVector predict(sgpp::base::DataMatrix& testDataset);
};

}  // namespace parallel
}  // namespace sgpp

#endif /* LEARNERVECTORIZEDIDENTITY_HPP */
