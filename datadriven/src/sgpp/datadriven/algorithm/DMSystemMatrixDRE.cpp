// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <sgpp/datadriven/algorithm/DMSystemMatrixDRE.hpp>

#include <sgpp/globaldef.hpp>

namespace sgpp {
namespace datadriven {

DMSystemMatrixDRE::DMSystemMatrixDRE(sgpp::base::DataMatrix& trainDataP,
                                     sgpp::base::DataMatrix& trainDataQ, double lambda)
    : datasetP_(trainDataP),
      datasetQ_(trainDataQ),
      lambda_(lambda),
      completeTimeMult_(0.0),
      computeTimeMult_(0.0),
      completeTimeMultTrans_(0.0),
      computeTimeMultTrans_(0.0) {
  myTimer_ = new sgpp::base::SGppStopwatch();
}

DMSystemMatrixDRE::~DMSystemMatrixDRE() { delete myTimer_; }

void DMSystemMatrixDRE::prepareGrid() {}

void DMSystemMatrixDRE::resetTimers() {
  completeTimeMult_ = 0.0;
  computeTimeMult_ = 0.0;
  completeTimeMultTrans_ = 0.0;
  computeTimeMultTrans_ = 0.0;
}

void DMSystemMatrixDRE::getTimers(double& timeMult, double& computeMult, double& timeMultTrans,
                                  double& computeMultTrans) {
  timeMult = completeTimeMult_;
  computeMult = computeTimeMult_;
  timeMultTrans = completeTimeMultTrans_;
  computeMultTrans = computeTimeMultTrans_;
}

}  // namespace datadriven
}  // namespace sgpp
