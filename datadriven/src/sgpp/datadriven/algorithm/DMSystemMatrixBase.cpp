// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <sgpp/datadriven/algorithm/DMSystemMatrixBase.hpp>

#include <sgpp/globaldef.hpp>

namespace sgpp {
namespace datadriven {

DMSystemMatrixBase::DMSystemMatrixBase(sgpp::base::DataMatrix& trainData, double lambda)
    : dataset_(trainData),
      lambda_(lambda),
      completeTimeMult_(0.0),
      computeTimeMult_(0.0),
      completeTimeMultTrans_(0.0),
      computeTimeMultTrans_(0.0) {
  myTimer_ = new sgpp::base::SGppStopwatch();
}

DMSystemMatrixBase::~DMSystemMatrixBase() { delete myTimer_; }

void DMSystemMatrixBase::prepareGrid() {}

void DMSystemMatrixBase::resetTimers() {
  completeTimeMult_ = 0.0;
  computeTimeMult_ = 0.0;
  completeTimeMultTrans_ = 0.0;
  computeTimeMultTrans_ = 0.0;
}

void DMSystemMatrixBase::getTimers(double& timeMult, double& computeMult, double& timeMultTrans,
                                   double& computeMultTrans) {
  timeMult = completeTimeMult_;
  computeMult = computeTimeMult_;
  timeMultTrans = completeTimeMultTrans_;
  computeMultTrans = computeTimeMultTrans_;
}

}  // namespace datadriven
}  // namespace sgpp
