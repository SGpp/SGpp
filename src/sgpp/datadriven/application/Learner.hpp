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

	bool verbose_;

	void createInitialGrid(sg::base::RegularGridConfiguration& GridConfig);

	void trainGrid(sg::base::DataMatrix& testDataset, sg::base::DataVector& classes, sg::base::DataMatrix& testDataset,
			sg::solver::SLESolverConfiguration& SolverConfigRefine, sg::solver::SLESolverConfiguration& SolverConfigFinal,
			sg::base::AdpativityConfiguration& AdaptConfig, sg::base::OperationMatrix& SLESystem);

	void trainGrid(sg::base::DataMatrix& testDataset, sg::base::DataVector& classes,
			sg::solver::SLESolverConfiguration& SolverConfig, sg::base::OperationMatrix& SLESystem);

public:
	Learner(bool verbose);

	virtual ~Learner();

	virtual void train(sg::base::DataMatrix& testDataset, sg::base::DataVector& classes,
			sg::base::RegularGridConfiguration& GridConfig, sg::base::AdpativityConfiguration& AdaptConfig,
			sg::solver::SLESolverConfiguration& SolverConfigRefine, sg::solver::SLESolverConfiguration& SolverConfigFinal);

	virtual void train(sg::base::DataMatrix& testDataset, sg::base::RegularGridConfiguration& GridConfig,
			sg::solver::SLESolverConfiguration& SolverConfig);

	virtual void test(sg::base::DataMatrix& testDataset, sg::base::DataVector& classes, bool isRegression);

	void store(std::string tGridFilename, std::string tAlphaFilename);

	sg::base::DataVector getAlpha();
};

}

}

#endif /* LEARNER_HPP */
