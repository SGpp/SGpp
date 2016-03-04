// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#ifndef LEARNERVECTORIZEDPERFORMANCECALCULATOR_HPP
#define LEARNERVECTORIZEDPERFORMANCECALCULATOR_HPP

#include <sgpp/base/grid/Grid.hpp>

#include <sgpp/solver/SLESolver.hpp>

#include <sgpp/globaldef.hpp>


namespace sgpp {
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
  static LearnerVectorizedPerformance getGFlopAndGByte(sgpp::base::Grid& Grid,
      size_t numInstances,
      sgpp::solver::SLESolverType solver, size_t numIterations, size_t sizeDatatype);
};

}  // namespace datadriven
}  // namespace sgpp

#endif /* LEARNERVECTORIZEDPERFORMANCECALCULATOR_HPP */
