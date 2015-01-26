/******************************************************************************
* Copyright (C) 2012 Technische Universitaet Muenchen                         *
* This file is part of the SG++ project. For conditions of distribution and   *
* use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
******************************************************************************/
// @author Alexander Heinecke (Alexander.Heinecke@mytum.de)

#include <sgpp/datadriven/algorithm/DMSystemMatrixBase.hpp>

namespace sg {
namespace datadriven {

DMSystemMatrixBase::DMSystemMatrixBase(sg::base::DataMatrix& trainData, double lambda)
    : dataset_(&trainData), lambda_(lambda), completeTimeMult_(0.0), computeTimeMult_(0.0),
      completeTimeMultTrans_(0.0), computeTimeMultTrans_(0.0) {
    myTimer_ = new sg::base::SGppStopwatch();
}

DMSystemMatrixBase::~DMSystemMatrixBase() {
    delete myTimer_;
}

void DMSystemMatrixBase::rebuildLevelAndIndex() {
}

void DMSystemMatrixBase::resetTimers() {
    completeTimeMult_ = 0.0;
    computeTimeMult_ = 0.0;
    completeTimeMultTrans_ = 0.0;
    computeTimeMultTrans_ = 0.0;
}

void DMSystemMatrixBase::getTimers(double& timeMult, double& computeMult, double& timeMultTrans, double& computeMultTrans) {
    timeMult = completeTimeMult_;
    computeMult = computeTimeMult_;
    timeMultTrans = completeTimeMultTrans_;
    computeMultTrans = computeTimeMultTrans_;
}

}

}
