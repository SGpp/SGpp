// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at 
// sgpp.sparsegrids.org

#include <sgpp/datadriven/algorithm/DMSystemMatrixBaseSP.hpp>

#include <sgpp/globaldef.hpp>


namespace SGPP {
  namespace datadriven {

    DMSystemMatrixBaseSP::DMSystemMatrixBaseSP(SGPP::base::DataMatrixSP& trainData, float lambda)
      : dataset_(&trainData), lambda_(lambda), completeTimeMult_(0.0), computeTimeMult_(0.0),
        completeTimeMultTrans_(0.0), computeTimeMultTrans_(0.0) {
      myTimer_ = new SGPP::base::SGppStopwatch();
    }

    DMSystemMatrixBaseSP::~DMSystemMatrixBaseSP() {
      delete myTimer_;
    }

    void DMSystemMatrixBaseSP::rebuildLevelAndIndex() {
    }

    void DMSystemMatrixBaseSP::resetTimers() {
      completeTimeMult_ = 0.0;
      computeTimeMult_ = 0.0;
      completeTimeMultTrans_ = 0.0;
      computeTimeMultTrans_ = 0.0;
    }

    void DMSystemMatrixBaseSP::getTimers(double& timeMult, double& computeMult, double& timeMultTrans, double& computeMultTrans) {
      timeMult = completeTimeMult_;
      computeMult = computeTimeMult_;
      timeMultTrans = completeTimeMultTrans_;
      computeMultTrans = computeTimeMultTrans_;
    }

  }

}