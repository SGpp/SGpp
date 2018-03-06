// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <sgpp/base/operation/BaseOpFactory.hpp>

#include <sgpp/pde/operation/PdeOpFactory.hpp>

#include <sgpp/base/operation/hash/OperationMatrix.hpp>

#include <sgpp/datadriven/application/Learner.hpp>
#include <sgpp/datadriven/algorithm/DMSystemMatrix.hpp>

#include <sgpp/globaldef.hpp>

#include <string>

namespace sgpp {
namespace datadriven {

Learner::Learner(sgpp::datadriven::RegularizationType& regularization, const bool isRegression,
                 const bool verbose)
    : LearnerBase(isRegression, verbose),
      CMode(regularization)
//, C_(NULL)
{}

// Learner::Learner(const std::string tGridFilename, const std::string tAlphaFilename,
//                 sgpp::datadriven::RegularizationType& regularization, const bool isRegression,
//                 const bool verbose)
//    : LearnerBase(tGridFilename, tAlphaFilename, isRegression, verbose),
//      CMode_(regularization)
//// ,C_(NULL)
//{}

Learner::~Learner() {
  //  if (C_ != NULL) delete C_;
}

std::unique_ptr<sgpp::datadriven::DMSystemMatrixBase> Learner::createDMSystem(
    sgpp::base::DataMatrix& trainDataset, double lambda) {
  base::OperationMatrix* C = nullptr;
  if (this->grid == NULL) return NULL;

  // Clean up, if needed
  //  if (C_ != NULL) delete C_;

  if (this->CMode == datadriven::RegularizationType::Laplace) {
    C = sgpp::op_factory::createOperationLaplace(*this->grid);
  } else if (this->CMode == datadriven::RegularizationType::Identity) {
    C = sgpp::op_factory::createOperationIdentity(*this->grid);
  } else {
    // should not happen
  }

  return std::make_unique<sgpp::datadriven::DMSystemMatrix>(*(this->grid), trainDataset, C, lambda);
}

}  // namespace datadriven
}  // namespace sgpp
