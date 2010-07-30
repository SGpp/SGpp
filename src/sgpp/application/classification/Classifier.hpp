/******************************************************************************
* Copyright (C) 2009 Technische Universitaet Muenchen                         *
* This file is part of the SG++ project. For conditions of distribution and   *
* use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
******************************************************************************/
// @author Alexander Heinecke (Alexander.Heinecke@mytum.de)

#ifndef CLASSIFIER_HPP
#define CLASSIFIER_HPP

#include "application/classification/Classifier.hpp"
#include "sgpp.hpp"

namespace sg
{

/**
 * class that provides a Simple Classifier based on a Sparse Grid
 */
class Classifier
{
private:
	/// The Sparse Grid needed in this classificator
	Grid* myGrid;
	/// the number of levels used for an regular grid
	size_t levels;
	/// the type of grid used
	std::string GridType;
	/// the dimension of the grid
	size_t dim;
	/// the number of instances
	size_t instancesNo;
	/// the number of test instances
	size_t testinstancesNo;
	/// regularisation parameter lambda
	double lambda;
	/// kind of stiffness matrix
	std::string StiffnessMode;
	/// maximum number of iterations that are executed to train the grid
	size_t IterationMax;
	/// epsilon needed in the cg method
	double epsilon;

	/**
	 * Creates a regular Sparse Grid
	 */
	void createRegularGrid();

	/**
	 * Trains the grid with a given training data
	 *
	 * @param alpha DataVector that contains the Sparse Grid's coefficients
	 * @param tfileTrain file containing the train data
	 */
	void trainGrid(DataVector& alpha, std::string tfileTrain);

	/**
	 * Tests the grid with a given test data
	 *
	 * @param alpha DataVector that contains the Sparse Grid's coefficients
	 * @param tfileTest the filename of the file containing the test data
	 * @return the precentage of correct classified instances
	 */
	double applyTestdata(DataVector& alpha, std::string tfileTest);

public:
	/**
	 * STD-Constructor
	 */
	Classifier();

	/**
	 * STD-Destructor
	 */
	~Classifier();

	/**
	 * Learns a classification problem based on some training data, afterwards
	 * the learned function is applied to an unknown test dataset.
	 *
	 * A regular Sparse Grid is used for the classification
	 *
	 * @param tfileTrain the filename of the file containing the train data
	 * @param tfileTest the filename of the file containing the test data
	 * @param level the number of levels used for the regular Sparse Grid
	 * @param lambda the regularisation parameter
	 * @param GridType the kind of grid that is used in the classification problem: N = no boundary, B = Boundaries, E = expanded inner functions, U = trapezoid boundaries (trapeze boundaries)
	 * @param StiffMode the kind of Stiffnessmatrix used: L: Laplacianmatrix, I: Identitymatrix
	 * @param epsilon epsilon needed in the cg method
	 * @param imax maximum number of iterations that are executed to train the grid
	 */
	void trainNtestRegular(std::string tfileTrain, std::string tfileTest, size_t level, double lambda, std::string GridType, std::string StiffMode, double epsilon, size_t imax);
};

}

#endif /* CLASSIFIER_HPP */
