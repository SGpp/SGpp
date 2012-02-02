/******************************************************************************
* Copyright (C) 2012 Technische Universitaet Muenchen                         *
* This file is part of the SG++ project. For conditions of distribution and   *
* use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
******************************************************************************/
// @author Alexander Heinecke (Alexander.Heinecke@mytum.de)

#include "parallel/datadriven/application/LearnerVectorizedSP.hpp"
#include "parallel/datadriven/algorithm/DMSystemMatrixSPVectorizedIdentity.hpp"
#include "parallel/datadriven/tools/LearnerVectorizedPerformanceCalculator.hpp"

namespace sg
{

namespace parallel
{

LearnerVectorizedSP::LearnerVectorizedSP(const VectorizationType vecType, const bool isRegression, const bool verbose)
	: sg::datadriven::LearnerBaseSP(isRegression, verbose), vecType_(vecType)
{
}

LearnerVectorizedSP::LearnerVectorizedSP(const std::string tGridFilename, const std::string tAlphaFilename,
		const VectorizationType vecType, const bool isRegression, const bool verbose)
	: sg::datadriven::LearnerBaseSP(tGridFilename, tAlphaFilename, isRegression, verbose), vecType_(vecType)
{
}

LearnerVectorizedSP::~LearnerVectorizedSP()
{
}

sg::datadriven::DMSystemMatrixBaseSP* LearnerVectorizedSP::createDMSystem(sg::base::DataMatrixSP& trainDataset, float lambda)
{
	if (this->grid_ == NULL)
		return NULL;

	return new sg::parallel::DMSystemMatrixSPVectorizedIdentity(*(this->grid_), trainDataset, lambda, this->vecType_);
}

void LearnerVectorizedSP::postProcessing(const sg::base::DataMatrixSP& trainDataset, const sg::solver::SLESolverType& solver,
		const size_t numNeededIterations)
{
	LearnerVectorizedPerformance currentPerf = LearnerVectorizedPerformanceCalculator::getGFlopAndGByte(*this->grid_,
			trainDataset.getNrows(), solver, numNeededIterations, sizeof(float));

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
