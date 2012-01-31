/******************************************************************************
* Copyright (C) 2012 Technische Universitaet Muenchen                         *
* This file is part of the SG++ project. For conditions of distribution and   *
* use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
******************************************************************************/
// @author Alexander Heinecke (Alexander.Heinecke@mytum.de)

#ifndef LEARNER_HPP
#define LEARNER_HPP

#include "base/grid/Grid.hpp"
#include "base/datatypes/DataVector.hpp"
#include "base/datatypes/DataMatrix.hpp"

#include "solver/SLESolver.hpp"

#include "datadriven/algorithm/BaseDMSystemMatrix.hpp"

namespace sg
{

namespace datadriven
{

struct ClassificatorQuality
{
    double truePositive;
    double trueNegative;
    double falsePositive;
    double falseNegative;
};

class Learner
{
protected:
	sg::base::DataVector* alpha_;

	sg::base::Grid* grid_;

	bool verbose_;

	bool isRegression_;

	virtual void postProcessing();

	void createInitialGrid(sg::base::RegularGridConfiguration& GridConfig);

	virtual void trainGrid(sg::base::DataMatrix& testDataset, sg::base::DataVector& classes, sg::base::DataMatrix& testDataset,
			sg::solver::SLESolverConfiguration& SolverConfigRefine, sg::solver::SLESolverConfiguration& SolverConfigFinal,
			sg::base::AdpativityConfiguration& AdaptConfig, sg::datadriven::BaseDMSystemMatrix& SLESystem,
			bool testAccDuringAdapt);

	virtual void trainGrid(sg::base::DataMatrix& testDataset, sg::base::DataVector& classes,
			sg::solver::SLESolverConfiguration& SolverConfig, sg::datadriven::BaseDMSystemMatrix& SLESystem);

public:
	Learner(bool isRegression, bool verbose);

	Learner(std::string tGridFilename, std::string tAlphaFilename, bool isRegression, bool verbose);

	virtual ~Learner();

	virtual void train(sg::base::DataMatrix& testDataset, sg::base::DataVector& classes,
			sg::base::RegularGridConfiguration& GridConfig, sg::base::AdpativityConfiguration& AdaptConfig,
			sg::solver::SLESolverConfiguration& SolverConfigRefine, sg::solver::SLESolverConfiguration& SolverConfigFinal,
			bool testAccDuringAdapt) = 0;

	virtual void train(sg::base::DataMatrix& testDataset, sg::base::RegularGridConfiguration& GridConfig,
			sg::solver::SLESolverConfiguration& SolverConfig) = 0;

	/**
	 * executes a Regression test for a given dataset and returns the result
	 *
	 * @param testDataset dataset that is evaluated with the current learner
	 *
	 * @return regression values of testDataset
	 */
	virtual sg::base::DataVector test(sg::base::DataMatrix& testDataset);

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
	virtual double getAccuracy(sg::base::DataMatrix& testDataset, sg::base::DataVector& classesReference, double threshold = 0.0)

	/**
	 * compute the accuracy for given testDataset.
	 *
	 * In case if classification (isRegression == false) this routine returns the learner's accuracy
	 * In case of regressions (isRegression == true) this routine returns the learner's MSE
	 *
	 * @param testDataset dataset to be tested
	 * @param classesComputed regression results of the test dataset
	 * @param classesReference reference labels of the test dataset
	 * @param threshold threshold used for classification, ignored when performing regressions
	 *
	 * @return accuracy, percent or MSE, depending on the execution mode
	 */
	virtual double getAccuracy(sg::base::DataMatrix& testDataset, sg::base::DataVector& classesComputed, sg::base::DataVector& classesReference, double threshold = 0.0)

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
	virtual ClassificatorQuality getCassificatorQuality(sg::base::DataMatrix& testDataset, sg::base::DataVector& classesReference, double threshold = 0.0)

	/**
	 * compute the quality for given testDataset, classification ONLY!
	 *
	 * @param testDataset dataset to be tested
	 * @param classesComputed regression results of the test dataset
	 * @param classesReference reference labels of the test dataset
	 * @param threshold threshold used for classification, ignored when performing regressions
	 *
	 * @return quality structure containing tp, tn, fp, fn counts
	 */
	virtual ClassificatorQuality getCassificatorQuality(sg::base::DataMatrix& testDataset, sg::base::DataVector& classesComputed, sg::base::DataVector& classesReference, double threshold = 0.0)

	/**
	 * store the grid and its current coefficients into files for
	 * further usage.
	 *
	 * @param tGridFilename filename of grid file
	 * @param tAlphaFilename filename of coefficient file
	 */
	void store(std::string tGridFilename, std::string tAlphaFilename);

	/**
	 * determines the current mode
	 *
	 * @return returns whether the current mode is regression or not
	 */
	bool getIsRegression();

	/**
	 * determines the current verbose mode of learner
	 *
	 * @return returns whether the current learner has verbose output
	 */
	bool getVerbose();

	/**
	 * sets the current verbose mode of learner
	 *
	 * @param verbose the current learner's verbose output
	 */
	bool setVerbose(bool verbose);
};

}

}

#endif /* LEARNER_HPP */
