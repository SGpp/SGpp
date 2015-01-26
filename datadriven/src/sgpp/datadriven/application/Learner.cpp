/******************************************************************************
* Copyright (C) 2012 Technische Universitaet Muenchen                         *
* This file is part of the SG++ project. For conditions of distribution and   *
* use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
******************************************************************************/
// @author Alexander Heinecke (Alexander.Heinecke@mytum.de)

#include <sgpp/base/operation/BaseOpFactory.hpp>

#include <sgpp/pde/operation/PdeOpFactory.hpp>

#include <sgpp/datadriven/application/Learner.hpp>
#include <sgpp/datadriven/algorithm/DMSystemMatrix.hpp>

namespace sg {
  namespace datadriven {

    Learner::Learner(sg::datadriven::LearnerRegularizationType& regularization, const bool isRegression, const bool verbose)
      : LearnerBase(isRegression, verbose), CMode_(regularization), C_(NULL) {
    }

    Learner::Learner(const std::string tGridFilename, const std::string tAlphaFilename, sg::datadriven::LearnerRegularizationType& regularization,
                     const bool isRegression, const bool verbose)
      : LearnerBase(tGridFilename, tAlphaFilename, isRegression, verbose), CMode_(regularization), C_(NULL) {
    }

    Learner::~Learner() {
      if (C_ != NULL)
        delete C_;
    }

    sg::datadriven::DMSystemMatrixBase* Learner::createDMSystem(sg::base::DataMatrix& trainDataset, double lambda) {
      if (this->grid_ == NULL)
        return NULL;

      // Clean up, if needed
      if (C_ != NULL)
        delete C_;

      if (this->CMode_ == Laplace) {
        C_ = sg::op_factory::createOperationLaplace(*this->grid_);
      } else if (this->CMode_ == Identity) {
        C_ = sg::op_factory::createOperationIdentity(*this->grid_);
      } else {
        // should not happen
      }

      return new sg::datadriven::DMSystemMatrix(*(this->grid_), trainDataset, *C_, lambda);
    }

  }

}
