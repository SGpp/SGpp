/* ****************************************************************************
* Copyright (C) 2009 Technische Universitaet Muenchen                         *
* This file is part of the SG++ project. For conditions of distribution and   *
* use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
**************************************************************************** */
#ifndef LEARNERDENSITYBASED_HPP_
#define LEARNERDENSITYBASED_HPP_

#include "datadriven/application/LearnerBase.hpp"
#include "datadriven/application/Learner.hpp"
#include "base/grid/type/LinearGrid.hpp"
#include "base/grid/type/LinearTrapezoidBoundaryGrid.hpp"
#include "base/grid/type/ModLinearGrid.hpp"

namespace sg
{

namespace datadriven
{


class LearnerDensityBased: public sg::datadriven::LearnerBase
{
protected:
	//Mapping from class index to class number:
	std::map<int, double> index_to_class_;
	//Stores the coefficients for every class
	std::vector<sg::base::DataVector> alphaVec_;
	/// regularization mode
	sg::datadriven::LearnerRegularizationType CMode_;
        // prior of data
  std::vector<double> prior;
  // vectors of grids
  std::vector<sg::base::Grid*> gridVec_;
  // vector of regterms
  std::vector<sg::base::OperationMatrix*> CVec_;
public:
	LearnerDensityBased(sg::datadriven::LearnerRegularizationType&, const bool isRegression, const bool isVerbose = true);
	virtual ~LearnerDensityBased();

  /**
   * Create a grid for each class
   *
   * @param GridConfig grid config
   * @param nrClasses number of classes
   */
  virtual void InitializeGrid(const sg::base::RegularGridConfiguration& GridConfig, size_t nrClasses);

  	/**
	 * Learning a dataset with spatially adaptive sparse grids
	 *
	 * @param testDataset the training dataset
	 * @param classes classes corresponding to the training dataset
	 * @param GridConfig configuration of the regular start grid
	 * @param SolverConfigRefine configuration of the SLE solver during the adaptive refinements of the grid
	 * @param SolverConfigFinal configuration of the final SLE solving step on the refined grid
	 * @param AdaptConfig configuration of the adaptivity strategy
	 * @param testAccDuringAdapt set to true if the training accuracy should be determined in evert refinement step
	 * @param lambda regularization parameter lambda
	 */
	virtual LearnerTiming train(sg::base::DataMatrix& testDataset, sg::base::DataVector& classes,
			const sg::base::RegularGridConfiguration& GridConfig, const sg::solver::SLESolverConfiguration& SolverConfigRefine,
			const sg::solver::SLESolverConfiguration& SolverConfigFinal, const sg::base::AdpativityConfiguration& AdaptConfig,
				    bool testAccDuringAdapt, const double lambda);


	/**
	 * Learning a dataset with spatially adaptive sparse grids
	 *
	 * @param testDataset the training dataset
	 * @param classes classes corresponding to the training dataset
	 * @param GridConfig configuration of the regular start grid
	 * @param SolverConfigRefine configuration of the SLE solver during the adaptive refinements of the grid
	 * @param SolverConfigFinal configuration of the final SLE solving step on the refined grid
	 * @param AdaptConfig configuration of the adaptivity strategy
	 * @param testAccDuringAdapt set to true if the training accuracy should be determined in evert refinement step
	 * @param lambda regularization parameter lambda
	 * @param usePrior use prior information for prediction (Bayes)
	 */
	LearnerTiming train(sg::base::DataMatrix& testDataset, sg::base::DataVector& classes,
			const sg::base::RegularGridConfiguration& GridConfig, const sg::solver::SLESolverConfiguration& SolverConfigRefine,
			const sg::solver::SLESolverConfiguration& SolverConfigFinal, const sg::base::AdpativityConfiguration& AdaptConfig,
				    bool testAccDuringAdapt, const double lambda, bool usePrior = true);
	virtual sg::base::DataVector predict(sg::base::DataMatrix& testDataset);
	/// construct system matrix
        virtual sg::datadriven::DMSystemMatrixBase* createDMSystem(sg::base::DataMatrix& trainDataset, double lambda);

  /**
   * Returns the execution time
   */
  time_t getExecTime();

  /**
   * Returns number of grid points for the density
   * with the maximum number of grid points
   */
  size_t getNrGridPoints();
};

}

}

#endif /* LEARNERDENSITYBASED_HPP_ */
