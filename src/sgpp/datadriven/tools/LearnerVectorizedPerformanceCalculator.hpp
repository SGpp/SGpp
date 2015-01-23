/* ****************************************************************************
 * Copyright (C) 2012 Technische Universitaet Muenchen                         *
 * This file is part of the SG++ project. For conditions of distribution and   *
 * use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
 **************************************************************************** */
// @author Alexander Heinecke (Alexander.Heinecke@mytum.de), David Pfander (David.Pfander@ipvs.uni-stuttgart.de)
#ifndef LEARNERVECTORIZEDPERFORMANCECALCULATOR_HPP
#define LEARNERVECTORIZEDPERFORMANCECALCULATOR_HPP

#include "base/grid/Grid.hpp"

#include "solver/SLESolver.hpp"

namespace sg {
namespace datadriven {

/**
 * struct that defines return
 * for calculation the performance
 * of a vectorized learner
 */
struct LearnerVectorizedPerformance {
	/// achieved GFLOP
	double GFlop_;
	/// achieved GByte
	double GByte_;
};

/**
 * Class that provides functionality
 * in order to determine a LearnerVectorized's
 * performance.
 */
class LearnerVectorizedPerformanceCalculator {
public:
	/**
	 * Calculate the performance of LearnerVectorized
	 *
	 * @param Grid reference to grid used bt the Learner
	 * @param numInstances number of training instances
	 * @param numIterations number of iterations the solver executed
	 * @param solver the selected solver
	 * @param sizeDatatype the size of the employed datatype in bytes
	 *
	 * @return a LearnerVectorizedPerformance struct containing the results
	 */
	static LearnerVectorizedPerformance getGFlopAndGByte(sg::base::Grid& Grid, size_t numInstances,
			sg::solver::SLESolverType solver, size_t numIterations, size_t sizeDatatype);
};

}
}

#endif /* LEARNERVECTORIZEDPERFORMANCECALCULATOR_HPP */

