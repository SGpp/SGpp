/******************************************************************************
* Copyright (C) 2012 Technische Universitaet Muenchen                         *
* This file is part of the SG++ project. For conditions of distribution and   *
* use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
******************************************************************************/
// @author Alexander Heinecke (Alexander.Heinecke@mytum.de)

#include "parallel/datadriven/application/LearnerVectorized.hpp"
#include "parallel/datadriven/algorithm/DMSystemMatrixVectorizedIdentity.hpp"
#include "parallel/datadriven/tools/LearnerVectorizedPerformanceCalculator.hpp"

namespace sg
{

namespace parallel
{

LearnerVectorized::LearnerVectorized(const VectorizationType vecType, const bool isRegression, const bool verbose)
	: sg::datadriven::LearnerBase(isRegression, verbose), vecType_(vecType)
{
}

LearnerVectorized::LearnerVectorized(const std::string tGridFilename, const std::string tAlphaFilename,
		const VectorizationType vecType, const bool isRegression, const bool verbose)
	: sg::datadriven::LearnerBase(tGridFilename, tAlphaFilename, isRegression, verbose), vecType_(vecType)
{
}

LearnerVectorized::~LearnerVectorized()
{
}

sg::datadriven::DMSystemMatrixBase* LearnerVectorized::createDMSystem(sg::base::DataMatrix trainDataset, double lambda)
{
	if (this->grid_ == NULL)
		return NULL;

	return new sg::parallel::DMSystemMatrixVectorizedIdentity(*(this->grid_), trainDataset, lambda, this->vecType_);
}

void LearnerVectorized::postProcessing(const sg::base::DataMatrix& trainDataset, const sg::solver::SLESolverType& solver,
		const size_t numNeededIterations)
{
	LearnerVectorizedPerformanceCalculator PerfCalc;

	LearnerVectorizedPerformance currentPerf = PerfCalc.getGFlopAndGByte(*this->grid_,
			trainDataset.getNrows(), solver, numNeededIterations, sizeof(double));

	this->GFlop_ = currentPerf.GFlop_;
	this->GByte_ = currentPerf.GByte_;

	// Caluate GFLOPS and GBytes/s and write them to console
	if (this->isVerbose_)
	{
    	std::cout << std::endl;
        std::cout << "Current GFlop/s: " << this->GFlop_/this->execTime_ << std::endl;
        std::cout << "Current GByte/s: " << this->GByte_/this->execTime_ << std::endl;
        std::cout << std::endl;
	}
}

}

}
