// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <sgpp/base/operation/BaseOpFactory.hpp>

#include <sgpp/pde/operation/PdeOpFactory.hpp>

#include <sgpp/datadriven/application/Learner.hpp>
#include <sgpp/datadriven/algorithm/DMSystemMatrix.hpp>

#include <sgpp/globaldef.hpp>


namespace SGPP {
  namespace datadriven {

    Learner::Learner(SGPP::datadriven::LearnerRegularizationType& regularization, const bool isRegression, const bool verbose)
      : LearnerBase(isRegression, verbose), CMode_(regularization), C_(NULL) {
    }

    Learner::Learner(const std::string tGridFilename, const std::string tAlphaFilename, SGPP::datadriven::LearnerRegularizationType& regularization,
                     const bool isRegression, const bool verbose)
      : LearnerBase(tGridFilename, tAlphaFilename, isRegression, verbose), CMode_(regularization), C_(NULL) {
    }

    Learner::~Learner() {
      if (C_ != NULL)
        delete C_;
    }

    SGPP::datadriven::DMSystemMatrixBase* Learner::createDMSystem(SGPP::base::DataMatrix& trainDataset, float_t lambda) {
      if (this->grid_ == NULL)
        return NULL;

      // Clean up, if needed
      if (C_ != NULL)
        delete C_;

      if (this->CMode_ == Laplace) {
        C_ = SGPP::op_factory::createOperationLaplace(*this->grid_);
      } else if (this->CMode_ == Identity) {
        C_ = SGPP::op_factory::createOperationIdentity(*this->grid_);
      } else {
        // should not happen
      }

      return new SGPP::datadriven::DMSystemMatrix(*(this->grid_), trainDataset, *C_, lambda);
    }

  }

}