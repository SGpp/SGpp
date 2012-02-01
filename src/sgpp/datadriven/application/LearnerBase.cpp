/******************************************************************************
* Copyright (C) 2012 Technische Universitaet Muenchen                         *
* This file is part of the SG++ project. For conditions of distribution and   *
* use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
******************************************************************************/
// @author Alexander Heinecke (Alexander.Heinecke@mytum.de)

#include "datadriven/application/Learner.hpp"

namespace sg
{

namespace datadriven
{

LearnerBase::LearnerBase(bool isRegression, bool verbose) : alpha_(NULL), grid_(NULL), verbose_(verbose), isRegression_(isRegression), isTrained_(false)
{
}

LearnerBase::LearnerBase(std::string tGridFilename, std::string tAlphaFilename, bool isRegression, bool verbose) : alpha_(NULL), grid_(NULL), verbose_(verbose), isRegression_(isRegression), isTrained_(false)
{
	// @TODO (heinecke)
}

LearnerBase::~LearnerBase()
{
	if (alpha_ != NULL)
		delete alpha_;

	if (grid_ != NULL)
		delete grid_;
}

void LearnerBase::createInitialGrid(sg::base::RegularGridConfiguration& GridConfig)
{
	if (GridConfig.type_ == sg::base::LinearTrapezoidBoundary)
	{
		grid_ = new sg::base::LinearTrapezoidBoundaryGrid(GridConfig.dim_);
	}
	else if (GridConfig.type_ == sg::base::ModLinear)
	{
		grid_ = new sg::base::ModLinearGrid(GridConfig.dim_);
	}
	else if (GridConfig.type_ == sg::base::Linear)
	{
		grid_ = new sg::base::LinearGrid(GridConfig.dim_);
	}
	else
	{
		std::cout << std::endl << "An unsupported grid type was chosen! Exiting...." << std::endl << std::endl;
		grid_ = NULL;
		return;
	}

	// Generate regular Grid with LEVELS Levels
	sg::base::GridGenerator* myGenerator = grid_->createGridGenerator();
	myGenerator->regular(GridConfig.level_);
	delete myGenerator;
}

void LearnerBase::postProcessing()
{
}

void LearnerBase::trainGrid(sg::base::DataMatrix& trainDataset, sg::base::DataVector& classes,
		sg::solver::SLESolverConfiguration& SolverConfigRefine, sg::solver::SLESolverConfiguration& SolverConfigFinal,
		sg::base::AdpativityConfiguration& AdaptConfig, sg::datadriven::DMSystemMatrixBase& SLESystem,
		bool testAccDuringAdapt)
{
    double execTime = 0.0;
    double oldAcc = 0.0;

	// Generate CG to solve System
	sg::solver::SLESolver* myCG;

	if (SolverConfigRefine.type_ != sg::solver::CG)
	{
		myCG = new sg::solver::ConjugateGradients(SolverConfigRefine.maxIterations_, SolverConfigRefine.eps_);
	}
	else if (SolverConfigRefine.type_ != sg::solver::BiCGSTAB)
	{
		myCG = new sg::solver::BiCGStab(SolverConfigRefine.maxIterations_, SolverConfigRefine.eps_);
	}
	else
	{
		// @TODO: Error
	}

	if (verbose_)
		std::cout << "Starting Learning...." << std::endl;

    // execute adaptsteps
    sg::base::SGppStopwatch* myStopwatch = new sg::base::SGppStopwatch();

    for (size_t i = 0; i < AdaptConfig.numRefinements_+1; i++)
    {
    	if (verbose_)
    		std::cout << std::endl << "Doing refinement: " << i << std::endl;

    	myStopwatch->start();

    	// Do Refinements
    	if (i > 0)
    	{
    		sg::base::SurplusRefinementFunctor* myRefineFunc = new sg::base::SurplusRefinementFunctor(alpha_, AdaptConfig.noPoints_, AdaptConfig.threshold_);
    		grid_->createGridGenerator()->refine(myRefineFunc);
    		delete myRefineFunc;

    		SLESystem.rebuildLevelAndIndex();

    		if (verbose_)
    			std::cout << "New Grid Size: " << grid_->getSize() << std::endl;

    		alpha_->resizeZero(grid_->getSize());
    	}
    	else
    	{
    		std::cout << "Grid Size: " << grid_->getSize() << std::endl;
    	}

    	sg::base::DataVector b(alpha_->getSize());
    	SLESystem.generateb(classes, b);

    	if (i == AdaptConfig.numRefinements_)
    	{
    		myCG->setMaxIterations(SolverConfigFinal.maxIterations_);
    		myCG->setEpsilon(SolverConfigFinal.eps_);
    	}
    	myCG->solve(SLESystem, *alpha_, b, true, false, 0.0);

        execTime += myStopwatch->stop();

        if (verbose_)
        {
        	std::cout << "Needed Iterations: " << myCG->getNumberIterations() << std::endl;
        	std::cout << "Final residuum: " << myCG->getResiduum() << std::endl;
        }

        postProcessing();

        if(testAccDuringAdapt || i == AdaptConfig.numRefinements_)
        {
 			double acc = getAccuracy(trainDataset, classes);

 			if (verbose_)
 			{
 				if (isRegression_)
 					std::cout << "MSE (train): " << acc << std::endl;
 				else
 					std::cout << "Acc (train): " << acc << std::endl;
 			}

 			if (isRegression_)
 			{
				if ((i > 0) && (oldAcc <= acc))
				{
					if (verbose_)
						std::cout << "The grid is becoming worse --> stop learning" << std::endl;

					break;
				}
 			}
 			else
 			{
				if ((i > 0) && (oldAcc >= acc))
				{
					if (verbose_)
						std::cout << "The grid is becoming worse --> stop learning" << std::endl;

					break;
				}
 			}

			oldAcc = acc;
        }
    }

    if (verbose_)
    {
    	std::cout << "Finished Training!" << std::endl << std::endl;
    	std::cout << "Training took: " << execTime << " seconds" << std::endl << std::endl;
    }

    isTrained_ = true;

    delete myStopwatch;
    delete myCG;
}

void LearnerBase::trainGrid(sg::base::DataMatrix& trainDataset, sg::base::DataVector& classes,
		sg::solver::SLESolverConfiguration& SolverConfig, sg::datadriven::DMSystemMatrixBase& SLESystem)
{
	sg::base::AdpativityConfiguration AdaptConfig;

	AdaptConfig.maxLevelType_ = false;
	AdaptConfig.noPoints_ = 0;
	AdaptConfig.numRefinements_ = 0;
	AdaptConfig.percent_ = 0.0;
	AdaptConfig.threshold_ = 0.0;

	trainGrid(trainDataset, classes, SolverConfig, SolverConfig, AdaptConfig, SLESystem, false);
}

sg::base::DataVector LearnerBase::test(sg::base::DataMatrix& testDataset)
{
	sg::base::DataVector classesComputed(testDataset.getNrows());

	sg::base::OperationMultipleEval* MultEval = sg::op_factory::createOperationMultipleEval(*grid_, &testDataset);
	MultEval->mult(*alpha_, classesComputed);
	delete MultEval;

	return classesComputed;
}

void LearnerBase::store(std::string tGridFilename, std::string tAlphaFilename)
{
	// @TODO (heinecke)
}

double LearnerBase::getAccuracy(sg::base::DataMatrix& testDataset, const sg::base::DataVector& classesReference, const double threshold = 0.0)
{
	// evaluate test dataset
	sg::base::DataVector classesComputed = test(testDataset);

	return getAccuracy(classesComputed, classesReference, threshold);
}

double LearnerBase::getAccuracy(const sg::base::DataVector& classesComputed, const sg::base::DataVector& classesReference, const double threshold = 0.0)
{
	double result = -1.0;

	if (classesComputed.getSize() != classesReference.getSize())
	{
		// @TODO (heinecke) throw exception
	}

	if (isRegression_)
	{
		sg::base::DataVector tmp(classesComputed);
	    tmp.sub(classesReference);
	    tmp.sqr();
	    result = tmp.sum();
	    result /= static_cast<double>(tmp.getSize());
	}
	else
	{
		size_t correct = 0;

		for(size_t i = 0; i < classesComputed.getSize(); i++)
		{
			if( (classesComputed[i] >= threshold && classesReference[i] >= 0.0)
					|| (classesComputed[i] < threshold && classesReference[i] < 0.0) )
			{
				correct++;
			}
		}

		result = static_cast<double>(correct)/static_cast<double>(classesComputed.getSize());
	}

	return result;
}

ClassificatorQuality LearnerBase::getCassificatorQuality(sg::base::DataMatrix& testDataset, const sg::base::DataVector& classesReference, const double threshold = 0.0)
{
	// evaluate test dataset
	sg::base::DataVector classesComputed = test(testDataset);

	return getCassificatorQuality(classesComputed, classesReference, threshold);
}

ClassificatorQuality LearnerBase::getCassificatorQuality(const sg::base::DataVector& classesComputed, const sg::base::DataVector& classesReference, const double threshold = 0.0)
{
	ClassificatorQuality result;

	if (isRegression_)
	{
		// @TODO (heinecke) throw exception
	}

	result.truePositive_ = 0;
	result.trueNegative_ = 0;
	result.falsePositive_ = 0;
	result.falseNegative_ = 0;

	for(size_t i = 0; i < classesComputed.getSize(); i++)
	{
		if( (classesComputed[i] >= threshold && classesReference[i] >= 0) )
		{
			result.truePositive_++;
		}
		else if( (classesComputed[i] < threshold && classesReference[i] < 0) )
		{
			result.trueNegative_++;
		}
		else if( (classesComputed[i] >= threshold && classesReference[i] < 0) )
		{
			result.falsePositive_++;
		}
		else // ( (classesComputed[i] < threshold && classesReference[i] >= 0) )
		{
			result.falseNegative_++;
		}
	}

	return result;
}

bool LearnerBase::getIsRegression() const
{
	return isRegression_;
}

bool LearnerBase::getVerbose() const
{
	return verbose_;
}

void LearnerBase::setVerbose(bool verbose)
{
	verbose_ = verbose;
}

}

}
