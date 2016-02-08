// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#ifndef LEARNER_HPP
#define LEARNER_HPP

#include <sgpp/datadriven/application/LearnerBase.hpp>
#include <sgpp/pde/application/RegularizationConfiguration.hpp>
#include <sgpp/globaldef.hpp>


namespace SGPP {

namespace datadriven {

/**
 * This class implements standard sparse grid regression
 * with an arbitrary regularization operator
 */
class Learner : public LearnerBase {
 protected:
  /// regularization mode
  SGPP::pde::RegularizationType CMode_;
  /// regularization operator
  SGPP::base::OperationMatrix* C_;

  /// construct system matrix
  virtual SGPP::datadriven::DMSystemMatrixBase* createDMSystem(
    SGPP::base::DataMatrix& trainDataset, float_t lambda);

 public:
  /**
   * Constructor
   *
   * @param regularization enum that gives the regurlarization method
   * @param isRegression flag for regression
   * @param isVerbose flag for verbose output
   */
  Learner(SGPP::pde::RegularizationType& regularization, const bool isRegression,
          const bool isVerbose = true);

  /**
   * Constructor
   *
   * @param tGridFilename path to file that contains a serialized grid
   * @param tAlphaFilename path to file that contains the grid's coefficients
   * @param regularization enum that gives the regurlarization method
   * @param isRegression set to true if a regression task should be executed
   * @param isVerbose set to true in order to allow console output
   */
  Learner(const std::string tGridFilename, const std::string tAlphaFilename,
          SGPP::pde::RegularizationType& regularization,
          const bool isRegression, const bool isVerbose = true);

  /**
   * Destructor
   */
  virtual ~Learner();
};

}

}

#endif /* LEARNER_HPP */
