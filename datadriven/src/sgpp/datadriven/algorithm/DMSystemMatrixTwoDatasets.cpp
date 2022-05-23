// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <sgpp/datadriven/algorithm/DMSystemMatrixTwoDatasets.hpp>
#include <sgpp/globaldef.hpp>

namespace sgpp {
namespace datadriven {

DMSystemMatrixTwoDatasets::DMSystemMatrixTwoDatasets(sgpp::base::DataMatrix& trainDataP,
                                                     sgpp::base::DataMatrix& trainDataQ,
                                                     double lambda)
    : datasetP(trainDataP),
      datasetQ(trainDataQ),
      lambda(lambda),
      completeTimeMult(0.0),
      computeTimeMult(0.0),
      completeTimeMultTrans(0.0),
      computeTimeMultTrans(0.0) {
  myTimer = new sgpp::base::SGppStopwatch();
}

DMSystemMatrixTwoDatasets::~DMSystemMatrixTwoDatasets() { delete myTimer; }

void DMSystemMatrixTwoDatasets::prepareGrid() {}

void DMSystemMatrixTwoDatasets::resetTimers() {
  completeTimeMult = 0.0;
  computeTimeMult = 0.0;
  completeTimeMultTrans = 0.0;
  computeTimeMultTrans = 0.0;
}

void DMSystemMatrixTwoDatasets::getTimers(double& timeMult, double& computeMult,
                                          double& timeMultTrans, double& computeMultTrans) {
  timeMult = completeTimeMult;
  computeMult = computeTimeMult;
  timeMultTrans = completeTimeMultTrans;
  computeMultTrans = computeTimeMultTrans;
}

}  // namespace datadriven
}  // namespace sgpp
