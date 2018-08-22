// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#pragma once

#include <sgpp/datadriven/application/LearnerBase.hpp>
#include <sgpp/datadriven/configuration/RegularizationConfiguration.hpp>
#include <sgpp/globaldef.hpp>

#include <string>

namespace sgpp {
namespace datadriven {

/**
 * This class implements standard sparse grid regression
 * with an arbitrary regularization operator
 */
class Learner : public LearnerBase {
 protected:
  /// regularization mode
  sgpp::datadriven::RegularizationType CMode;

  /// construct system matrix
  virtual std::unique_ptr<sgpp::datadriven::DMSystemMatrixBase> createDMSystem(
      sgpp::base::DataMatrix& trainDataset, double lambda);

 public:
  /**
   * Constructor
   *
   * @param regularization enum that gives the regurlarization method
   * @param isRegression flag for regression
   * @param isVerbose flag for verbose output
   */
  Learner(sgpp::datadriven::RegularizationType& regularization, const bool isRegression,
          const bool isVerbose = true);

  /**
   * Destructor
   */
  virtual ~Learner();
};

}  // namespace datadriven
}  // namespace sgpp
