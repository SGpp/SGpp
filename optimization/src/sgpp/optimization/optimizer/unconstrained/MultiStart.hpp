// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#ifndef SGPP_OPTIMIZATION_OPTIMIZER_UNCONSTRAINED_MULTISTART_HPP
#define SGPP_OPTIMIZATION_OPTIMIZER_UNCONSTRAINED_MULTISTART_HPP

#include <sgpp/globaldef.hpp>

#include <sgpp/optimization/optimizer/unconstrained/UnconstrainedOptimizer.hpp>
#include <sgpp/optimization/optimizer/unconstrained/NelderMead.hpp>

#include <vector>

namespace SGPP {
namespace optimization {
namespace optimizer {

/**
 * Meta optimization algorithm calling local algorithm multiple times.
 * MultiStart generates a random population of a given number of
 * starting points, and then runs a local optimization algorithm
 * for each of the starting point.
 * The best point wins.
 */
class MultiStart : public UnconstrainedOptimizer {
 public:
  /// default maximal number of function evaluations
  static const size_t DEFAULT_MAX_FCN_EVAL_COUNT = 1000;

  /**
   * Constructor.
   * By default, Nelder-Mead is used as optimization algorithm.
   *
   * @param f               objective function
   * @param maxFcnEvalCount maximal number of function evaluations
   * @param populationSize  number of individual points
   *                        (default: \f$\min(10d, 100)\f$)
   */
  MultiStart(ScalarFunction& f, size_t maxFcnEvalCount = DEFAULT_MAX_FCN_EVAL_COUNT,
             size_t populationSize = 0);

  /**
   * Constructor with custom optimization algorithm.
   * The current values of the optimizer's N and starting point
   * properties will not be used.
   *
   * @param optimizer        optimization algorithm and
   *                         objective function
   * @param maxFcnEvalCount  maximal number of function evaluations
   * @param populationSize   number of individual points
   *                         (default: \f$\min(10d, 100)\f$)
   */
  MultiStart(UnconstrainedOptimizer& optimizer, size_t maxFcnEvalCount = DEFAULT_MAX_FCN_EVAL_COUNT,
             size_t populationSize = 0);

  /**
   * Destructor.
   */
  ~MultiStart() override;

  void optimize() override;

  /**
   * @return                  number of individual points
   */
  size_t getPopulationSize() const;

  /**
   * @param populationSize    number of individual points
   */
  void setPopulationSize(size_t populationSize);

  /**
   * @return vector in which the k-th entry indicates the number of
   *         inner iterations in the k-th (outer) iteration,
   *         empty vector on error
   */
  const std::vector<size_t>& getHistoryOfInnerIterations() const;

  /**
   * @param[out] clone pointer to cloned object
   */
  void clone(std::unique_ptr<UnconstrainedOptimizer>& clone) const override;

 protected:
  /// default optimization algorithm
  NelderMead defaultOptimizer;
  /// optimization algorithm
  UnconstrainedOptimizer& optimizer;
  /// number of individual points
  size_t populationSize;
  /// search history (inner iterations)
  std::vector<size_t> kHist;

  /**
   * Initializes populationSize.
   *
   * @param populationSize     number of individual points
   *                           (zero to use default value)
   */
  void initialize(size_t populationSize);
};
}  // namespace optimizer
}  // namespace optimization
}  // namespace SGPP

#endif /* SGPP_OPTIMIZATION_OPTIMIZER_UNCONSTRAINED_MULTISTART_HPP */
