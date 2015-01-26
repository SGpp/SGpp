/* ****************************************************************************
 * Copyright (C) 2012 Technische Universitaet Muenchen                         *
 * This file is part of the SG++ project. For conditions of distribution and   *
 * use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
 **************************************************************************** */
// @author Alexander Heinecke (Alexander.Heinecke@mytum.de)
#ifndef LEARNERBASE_HPP
#define LEARNERBASE_HPP

#include "base/grid/Grid.hpp"
#include "base/datatypes/DataVector.hpp"
#include "base/datatypes/DataMatrix.hpp"

#include "solver/SLESolver.hpp"

#include "datadriven/algorithm/DMSystemMatrixBase.hpp"
#include "datadriven/tools/TypesDatadriven.hpp"

namespace sg {

namespace datadriven {

/**
 * Abstract class that implements a regression/classification learner
 * based on spatial adaptive Sparse Grids.
 *
 * Furthermore this class is intended to provide a framework
 * for adavanded regression and classification methods
 * by allowing to override basic routines like train
 * and test.
 */
class LearnerBase {
protected:
	/// the grid's coefficients
	sg::base::DataVector* alpha_;
	/// sparse grid object
	sg::base::Grid* grid_;
	/// is verbose output enabled
	bool isVerbose_;
	/// is regression selected
	bool isRegression_;
	/// is the grid trained
	bool isTrained_;
	/// execution time
	double execTime_;
	/// execution time for current refinement
	/// to calculate the GFlops at the current timestep only
	/// otherwise accumulated GFlops (all refinement steps)
	/// are calculated
	double stepExecTime_;
	/// number of executed Floating Point operations
	double GFlop_;
	/// number of executed Floating Point operations
	/// in the current refinement step
	double stepGFlop_;
	/// number of transferred Gbytes
	double GByte_;
	/// number of transferred Gbytes
	/// in the current refinement step
	double stepGByte_;
	/// the current refinment step during training
	size_t currentRefinementStep;

	std::vector<std::pair<size_t, double> > ExecTimeOnStep;

	/**
	 * Hook-Method for pre-processing before
	 * starting learning.
	 *
	 * can be overwritten by derived classes
	 */
	virtual void preProcessing();

	/**
	 * Hook-Method for post-processing after each
	 * refinement learning.
	 *
	 * can be overwritten by derived classes
	 * @param trainDataset matrix with training data
	 * @param solver solver
	 * @param numNeededIterations number of required iterations
	 *
	 */
	virtual void postProcessing(const sg::base::DataMatrix& trainDataset, const sg::solver::SLESolverType& solver,
			const size_t numNeededIterations);

	/**
	 * Initialize the grid and its coefficients
	 *
	 * @param GridConfig structure which describes the regular start grid
	 */
	virtual void InitializeGrid(const sg::base::RegularGridConfiguration& GridConfig);

	/**
	 * abstract method that constructs the corresponding system of linear equations
	 * Derived classes MUST overwrite this functions!
	 *
	 * @param trainDataset training dataset
	 * @param lambda lambda regularization parameter
	 */
	virtual sg::datadriven::DMSystemMatrixBase* createDMSystem(sg::base::DataMatrix& trainDataset, double lambda) = 0;

public:
	/**
	 * Constructor
	 *
	 * @param isRegression flag for regression
	 * @param isVerbose flag for verbose output
	 */
	LearnerBase(const bool isRegression, const bool isVerbose = true);

	/**
	 * Constructor
	 *
	 * @param tGridFilename path to file that contains a serialized grid
	 * @param tAlphaFilename path to file that contains the grid's coefficients
	 * @param isRegression set to true if a regression task should be executed
	 * @param isVerbose set to true in order to allow console output
	 */
	LearnerBase(std::string tGridFilename, std::string tAlphaFilename, const bool isRegression, const bool isVerbose =
			true);

	/**
	 * Copy-Constructor
	 *
	 * @param copyMe LearnerBase that should be duplicated
	 */
	LearnerBase(const LearnerBase& copyMe);

	/**
	 * Destructor
	 */
	virtual ~LearnerBase();

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
			const sg::base::RegularGridConfiguration& GridConfig,
			const sg::solver::SLESolverConfiguration& SolverConfigRefine,
			const sg::solver::SLESolverConfiguration& SolverConfigFinal,
			const sg::base::AdpativityConfiguration& AdaptConfig, bool testAccDuringAdapt, const double lambda);

	/**
	 * Learning a dataset with regular sparse grids
	 *
	 * @param testDataset the training dataset
	 * @param classes classes corresponding to the training dataset
	 * @param GridConfig configuration of the regular grid
	 * @param SolverConfig configuration of the SLE solver
	 * @param lamda regularization parameter lambda
	 */
	LearnerTiming train(sg::base::DataMatrix& testDataset, sg::base::DataVector& classes,
			const sg::base::RegularGridConfiguration& GridConfig,
			const sg::solver::SLESolverConfiguration& SolverConfig, const double lamda);

	/**
	 * executes a Regression test for a given dataset and returns the result
	 *
	 * @param testDataset dataset that is evaluated with the current learner
	 *
	 * @return regression values of testDataset
	 */
	virtual sg::base::DataVector predict(sg::base::DataMatrix& testDataset);

	/**
	 * compute the accuracy for given testDataset. test is automatically called
	 * in order to determine the regression values of the current learner
	 *
	 * In case if classification (isRegression == false) this routine returns the learner's accuracy
	 * In case of regressions (isRegression == true) this routine returns the learner's MSE
	 *
	 * @param testDataset dataset to be tested
	 * @param classesReference reference labels of the test dataset
	 * @param threshold threshold used for classification, ignored when performing regressions
	 *
	 * @return accuracy, percent or MSE, depending on the execution mode
	 */
	virtual double getAccuracy(sg::base::DataMatrix& testDataset, const sg::base::DataVector& classesReference,
			const double threshold = 0.0);

	/**
	 * compute the accuracy for given testDataset.
	 *
	 * In case if classification (isRegression == false) this routine returns the learner's accuracy
	 * In case of regressions (isRegression == true) this routine returns the learner's MSE
	 *
	 * @param classesComputed regression results of the test dataset
	 * @param classesReference reference labels of the test dataset
	 * @param threshold threshold used for classification, ignored when performing regressions
	 *
	 * @return accuracy, percent or MSE, depending on the execution mode
	 */
	virtual double getAccuracy(const sg::base::DataVector& classesComputed,
			const sg::base::DataVector& classesReference, const double threshold = 0.0);

	/**
	 * compute the quality for given testDataset, classification ONLY!
	 * test is automatically called
	 * in order to determine the regression values of the current learner
	 *
	 * @param testDataset dataset to be tested
	 * @param classesReference reference labels of the test dataset
	 * @param threshold threshold used for classification, ignored when performing regressions
	 *
	 * @return quality structure containing tp, tn, fp, fn counts
	 */
	virtual ClassificatorQuality getCassificatorQuality(sg::base::DataMatrix& testDataset,
			const sg::base::DataVector& classesReference, const double threshold = 0.0);

	/**
	 * compute the quality for given testDataset, classification ONLY!
	 *
	 * @param classesComputed regression results of the test dataset
	 * @param classesReference reference labels of the test dataset
	 * @param threshold threshold used for classification, ignored when performing regressions
	 *
	 * @return quality structure containing tp, tn, fp, fn counts
	 */
	virtual ClassificatorQuality getCassificatorQuality(const sg::base::DataVector& classesComputed,
			const sg::base::DataVector& classesReference, const double threshold = 0.0);

	/**
	 * store the grid and its current coefficients into files for
	 * further usage.
	 *
	 * @param tGridFilename filename of grid file
	 * @param tAlphaFilename filename of coefficient file
	 */
	void store(std::string tGridFilename, std::string tAlphaFilename);

	/**
	 * simple dump of grid points into file, e.g. used to
	 * plot with gnuplot
	 *
	 * only executed if grid is trained
	 *
	 * @param tFilename filename to store the dump to
	 */
	void dumpGrid(std::string tFilename);

	/**
	 * simple dump of sparse grid function into file, e.g. used to
	 * plot with gnuplot.
	 *
	 * only executed if grid is trained and number of dimensions <= 2.
	 *
	 * @param tFilename filename to store the dump to
	 * @param resolution resolution of function plot
	 */
	void dumpFunction(std::string tFilename, size_t resolution);

	/**
	 * determines the current mode
	 *
	 * @return returns whether the current mode is regression or not
	 */
	bool getIsRegression() const;

	/**
	 * determines the current verbose mode of learner
	 *
	 * @return returns whether the current learner has verbose output
	 */
	bool getIsVerbose() const;

	/**
	 * sets the current verbose mode of learner
	 *
	 * @param isVerbose the current learner's verbose output
	 */
	void setIsVerbose(const bool isVerbose);

	std::vector<std::pair<size_t, double> > getRefinementExecTimes();
};

}

}

#endif /* LEARNERBASE_HPP */
