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

namespace sg
{

namespace datadriven
{

class Learner
{
protected:
	sg::base::DataVector* alpha_;

	sg::base::Grid* grid_;

	void trainGrid(sg::base::DataMatrix& testDataset, sg::base::DataVector& classes,
			sg::solver::SLESolver& solver, sg::base::AdpativityConfiguration& AdaptConfig, double finalEps, size_t maxIterations);

public:
	Learner();

	~Learner();

	void train(sg::base::DataMatrix& testDataset, sg::base::RegularGridConfiguration& GridConfig,
			sg::base::AdpativityConfiguration& AdaptConfig, sg::solver::SLESolverConfiguration& SolverRefine,
			sg::solver::SLESolverConfiguration& SolverFinal);

	void train(sg::base::DataMatrix& testDataset, sg::base::RegularGridConfiguration& GridConfig,
			sg::solver::SLESolverConfiguration& SolverConfig);

	void test(sg::base::DataMatrix& testDataset, sg::base::DataVector& classes, bool isRegression);

	void store(std::string tGridFilename, std::string tAlphaFilename);

	sg::base::DataVector getAlpha();
};

}

}

#endif /* LEARNER_HPP */
