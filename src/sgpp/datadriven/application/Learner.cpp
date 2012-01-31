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

Learner::Learner(bool is Regression, bool verbose) : alpha_(NULL), grid_(NULL), verbose_(verbose), isRegression_(isRegression)
{
}

Learner::~Learner()
{
	if (alpha_ != NULL)
		delete alpha_;

	if (grid_ != NULL)
		delete grid_;
}

void Learner::createInitialGrid(sg::base::RegularGridConfiguration& GridConfig)
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

void Learner::postProcessing()
{
}

void Learner::trainGrid(sg::base::DataMatrix& trainDataset, sg::base::DataVector& classes, sg::base::DataMatrix& testDataset,
		sg::solver::SLESolverConfiguration& SolverConfigRefine, sg::solver::SLESolverConfiguration& SolverConfigFinal,
		sg::base::AdpativityConfiguration& AdaptConfig,sg::datadriven::BaseDMSystemMatrix& SLESystem,
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
    	std::cout << "Finished Learning!" << std::endl;

    delete myStopwatch;
    delete myCG;
}

void Learner::trainGrid(sg::base::DataMatrix& testDataset, sg::base::DataVector& classes,
		sg::solver::SLESolverConfiguration& SolverConfig, sg::datadriven::BaseDMSystemMatrix& SLESystem)
{

}

void Learner::train(sg::base::DataMatrix& testDataset, sg::base::DataVector& classes,
		sg::base::RegularGridConfiguration& GridConfig, sg::base::AdpativityConfiguration& AdaptConfig,
		sg::solver::SLESolverConfiguration& SolverConfigRefine, sg::solver::SLESolverConfiguration& SolverConfigFinal,
		bool testAccDuringAdapt)
{

}

void Learner::train(sg::base::DataMatrix& testDataset, sg::base::RegularGridConfiguration& GridConfig,
		sg::solver::SLESolverConfiguration& SolverConfig)
{

}

void sg::base::DataVector Learner::test(sg::base::DataMatrix& testDataset)
{
	sg::base::DataVector classesComputed(testDataset.getNroes());

	return classesComputed;
}

void Learner::store(std::string tGridFilename, std::string tAlphaFilename)
{

}

}

}
