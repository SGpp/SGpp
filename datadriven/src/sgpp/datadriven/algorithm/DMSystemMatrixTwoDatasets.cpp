// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <sgpp/globaldef.hpp>
#include <sgpp/datadriven/algorithm/DMSystemMatrixTwoDatasets.hpp>

namespace sgpp {
namespace datadriven {

DMSystemMatrixTwoDatasets::DMSystemMatrixTwoDatasets(sgpp::base::DataMatrix& trainDataP,
                                                     sgpp::base::DataMatrix& trainDataQ,
                                                     double lambda)
    : datasetP_(trainDataP),
      datasetQ_(trainDataQ),
      lambda_(lambda),
      completeTimeMult_(0.0),
      computeTimeMult_(0.0),
      completeTimeMultTrans_(0.0),
      computeTimeMultTrans_(0.0) {
  myTimer_ = new sgpp::base::SGppStopwatch();
}

DMSystemMatrixTwoDatasets::~DMSystemMatrixTwoDatasets() { delete myTimer_; }

void DMSystemMatrixTwoDatasets::prepareGrid() {}

void DMSystemMatrixTwoDatasets::resetTimers() {
  completeTimeMult_ = 0.0;
  computeTimeMult_ = 0.0;
  completeTimeMultTrans_ = 0.0;
  computeTimeMultTrans_ = 0.0;
}

void DMSystemMatrixTwoDatasets::getTimers(double& timeMult, double& computeMult,
                                          double& timeMultTrans, double& computeMultTrans) {
  timeMult = completeTimeMult_;
  computeMult = computeTimeMult_;
  timeMultTrans = completeTimeMultTrans_;
  computeMultTrans = computeTimeMultTrans_;
}

}  // namespace datadriven
}  // namespace sgpp
