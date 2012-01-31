/******************************************************************************
* Copyright (C) 2012 Technische Universitaet Muenchen                         *
* This file is part of the SG++ project. For conditions of distribution and   *
* use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
******************************************************************************/
// @author Alexander Heinecke (Alexander.Heinecke@mytum.de)

#include "datadriven/algorithm/BaseDMSystemMatrixSP.hpp"

namespace sg
{
namespace datadriven
{

BaseDMSystemMatrixSP::BaseDMSystemMatrixSP(sg::base::DataMatrixSP& trainData, float lambda)
	: dataset_(&trainData), lambda_(lambda), completeTimeMult_(0.0), computeTimeMult_(0.0),
	  completeTimeMultTrans_(0.0), computeTimeMultTrans_(0.0)
{
	myTimer_ = new sg::base::SGppStopwatch();
}

BaseDMSystemMatrixSP::~BaseDMSystemMatrixSP()
{
	delete myTimer_;
}

void BaseDMSystemMatrixSP::rebuildLevelAndIndex()
{
}

void BaseDMSystemMatrixSP::resetTimers()
{
	completeTimeMult_ = 0.0;
	computeTimeMult_ = 0.0;
	completeTimeMultTrans_ = 0.0;
	computeTimeMultTrans_ = 0.0;
}

void BaseDMSystemMatrixSP::getTimers(double& timeMult, double& computeMult, double& timeMultTrans, double& computeMultTrans)
{
	timeMult = completeTimeMult_;
	computeMult = computeTimeMult_;
	timeMultTrans = completeTimeMultTrans_;
	computeMultTrans = computeTimeMultTrans_;
}

}

}
