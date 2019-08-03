// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#ifndef SGPP_OPTIMIZATION_OPTIMIZER_UNCONSTRAINED_DIFFERENTIALEVOLUTION_HPP
#define SGPP_OPTIMIZATION_OPTIMIZER_UNCONSTRAINED_DIFFERENTIALEVOLUTION_HPP

#include <sgpp/globaldef.hpp>

#include <sgpp/optimization/optimizer/unconstrained/UnconstrainedOptimizer.hpp>

namespace sgpp {
namespace optimization {
namespace optimizer {

/**
 * Gradient-free Differential Evolution method.
 */
class DifferentialEvolution : public UnconstrainedOptimizer {
 public:
  /// default crossover probability
  static constexpr double DEFAULT_CROSSOVER_PROBABILITY = 0.5;
  /// default crossover scaling factor
  static constexpr double DEFAULT_SCALING_FACTOR = 0.6;
  /// default stopping criterion parameter 1
  static const size_t DEFAULT_IDLE_GENERATIONS_COUNT = 20;
  /// default stopping criterion parameter 2
  static constexpr double DEFAULT_AVG_IMPROVEMENT_THRESHOLD = 1e-6;
  /// default stopping criterion parameter 3
  static constexpr double DEFAULT_MAX_DISTANCE_THRESHOLD = 1e-4;

  /**
   * Constructor.
   *
   * @param f                         objective function
   * @param maxFcnEvalCount           maximal number of
   *                                  function evaluations
   * @param populationSize            number of individuals
   *                                  (default: \f$10d\f$)
   * @param crossoverProbability      crossover probability
   * @param scalingFactor             crossover scaling factor
   * @param idleGenerationsCount      stopping criterion parameter 1
   * @param avgImprovementThreshold   stopping criterion parameter 2
   * @param maxDistanceThreshold      stopping criterion parameter 3
   */
  DifferentialEvolution(const base::ScalarFunction& f, size_t maxFcnEvalCount = DEFAULT_N,
                        size_t populationSize = 0,
                        double crossoverProbability = DEFAULT_CROSSOVER_PROBABILITY,
                        double scalingFactor = DEFAULT_SCALING_FACTOR,
                        size_t idleGenerationsCount = DEFAULT_IDLE_GENERATIONS_COUNT,
                        double avgImprovementThreshold = DEFAULT_AVG_IMPROVEMENT_THRESHOLD,
                        double maxDistanceThreshold = DEFAULT_MAX_DISTANCE_THRESHOLD);

  /**
   * Copy constructor.
   *
   * @param other optimizer to be copied
   */
  DifferentialEvolution(const DifferentialEvolution& other);

  /**
   * Destructor.
   */
  ~DifferentialEvolution() override;

  void optimize() override;

  /**
   * @return                  number of individuals
   */
  size_t getPopulationSize() const;

  /**
   * @param populationSize    number of individuals
   */
  void setPopulationSize(size_t populationSize);

  /**
   * @param[out] clone pointer to cloned object
   */
  void clone(std::unique_ptr<UnconstrainedOptimizer>& clone) const override;

 protected:
  /// number of individuals
  size_t populationSize;
  /// crossover probability
  double crossoverProbability;
  /// crossover scaling factor
  double scalingFactor;
  /// stopping criterion parameter 1
  size_t idleGenerationsCount;
  /// stopping criterion parameter 2
  double avgImprovementThreshold;
  /// stopping criterion parameter 3
  double maxDistanceThreshold;
};
}  // namespace optimizer
}  // namespace optimization
}  // namespace sgpp

#endif /* SGPP_OPTIMIZATION_OPTIMIZER_UNCONSTRAINED_DIFFERENTIALEVOLUTION_HPP */
