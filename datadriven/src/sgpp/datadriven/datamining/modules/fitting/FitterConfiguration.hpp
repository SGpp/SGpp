// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#pragma once

#include <sgpp/base/grid/Grid.hpp>
#include <sgpp/datadriven/application/RegularizationConfiguration.hpp>
#include <sgpp/datadriven/datamining/configuration/DataMiningConfigParser.hpp>
#include <sgpp/datadriven/operation/hash/DatadrivenOperationCommon.hpp>
#include <sgpp/solver/TypesSolver.hpp>

namespace sgpp {
namespace datadriven {
// forward declaration to break dependency cycle
class DataMiningConfigParser;

/**
 * Different fitter scenarios have different default values and support different operations
 */
enum class FitterType { RegressionLeastSquares };

/**
 * General configuration object for fitters. Bundles all structures needed to build a sparse grid,
 * fit a sparse grid based model, and perform adaptive refinement.
 */
class FitterConfiguration {
 public:
  /**
   * Sets up a Fitter configuration with its default values.
   */
  FitterConfiguration();

  /**
   * Copy constructor
   * @param rhs const reference to the scorer object to copy from.
   */
  FitterConfiguration(const FitterConfiguration& rhs) = default;

  /**
   * Move constructor
   * @param rhs R-value reference to a scorer object to moved from.
   */
  FitterConfiguration(FitterConfiguration&& rhs) = default;

  /**
   * Copy assign operator
   * @param rhs const reference to the scorer object to copy from.
   * @return rerefernce to this with updated values.
   */
  FitterConfiguration& operator=(const FitterConfiguration& rhs) = default;

  /**
   * Move assign operator
   * @param rhs R-value reference to an a scorer object to move from.
   * @return rerefernce to this with updated values.
   */
  FitterConfiguration& operator=(FitterConfiguration&& rhs) = default;

  /**
   * virtual destructor.
   */
  virtual ~FitterConfiguration() = default;

  /**
   * Polymorphic clone pattern
   * @return deep copy of this object. New object is owned by caller.
   */
  virtual FitterConfiguration* clone() const = 0;

  /**
   * Get initial conditions for the grid before adaptive refinement.
   * @return immutable RegularGridConfiguration
   */
  const base::RegularGridConfiguration& getGridConfig() const;

  /**
   * Get how the adaptivity algorithms for the grid should behave.
   * @return immutable AdpativityConfiguration
   */
  const base::AdpativityConfiguration& getRefinementConfig() const;

  /**
   * Get configuration for the linear system solver which should be used while building
   * adaptive grids
   * @return immutable SLESolverConfiguration
   */
  const solver::SLESolverConfiguration& getSolverRefineConfig() const;

  /**
   * Get configuration for the linear system solver when solving the final, refined system
   * @return immutable SLESolverConfiguration
   */
  const solver::SLESolverConfiguration& getSolverFinalConfig() const;

  /**
   * Get the type of regularization operation to use
   * @return immutable RegularizationConfiguration
   */
  const datadriven::RegularizationConfiguration& getRegularizationConfig() const;

  /**
   * Get implementation (openMP, MPI, GPU) that should be used for
   * #sgpp::base::operationMultipleEval.
   * @return immutable OperationMultipleEvalConfiguration
   */
  const datadriven::OperationMultipleEvalConfiguration& getMultipleEvalConfig() const;

  /**
   * Get or set initial conditions for the grid before adaptive refinement.
   * @return RegularGridConfiguration
   */
  base::RegularGridConfiguration& getGridConfig();

  /**
   * Get or set how the adaptivity algorithms for the grid should behave.
   * @return AdpativityConfiguration
   */
  base::AdpativityConfiguration& getRefinementConfig();

  /**
   * Get or set configuration for the linear system solver which should be used while building
   * adaptive grids
   * @return SLESolverConfiguration
   */
  solver::SLESolverConfiguration& getSolverRefineConfig();

  /**
   * Get or set configuration for the linear system solver when solving the final, refined system
   * @return SLESolverConfiguration
   */
  solver::SLESolverConfiguration& getSolverFinalConfig();

  /**
   * Get or set the type of regularization operation to use
   * @return RegularizationConfiguration
   */
  datadriven::RegularizationConfiguration& getRegularizationConfig();

  /**
   * Get or set implementation (openMP, MPI, GPU) that should be used for
   * #sgpp::base::operationMultipleEval.
   * @return current OperationMultipleEvalConfiguration
   */
  datadriven::OperationMultipleEvalConfiguration& getMultipleEvalConfig();

  /**
   * Get the regularization weight.
   * @return regularization weight.
   */
  double getLambda();

  /**
   * Set the regularization weight.
   * @param lambda new regularization weight.
   */
  void setLambda(double lambda);

  /**
   * set default values for all members based on the desired scenario.
   */
  virtual void setupDefaults() = 0;

  /**
   * obtain parameters from a parser
   * @params parser: the parser object to read from
   */
  virtual void readParams(const DataMiningConfigParser& parser) = 0;

 protected:
  /**
   * Initial conditions for the grid before adaptive refinement.
   */
  base::RegularGridConfiguration gridConfig;

  /**
   * Configure how the adaptivity algorithms for the grid should behave.
   */
  base::AdpativityConfiguration adaptivityConfig;

  /**
   * Configuration for the linear system solver which should be used while building adaptive grids
   */
  solver::SLESolverConfiguration solverRefineConfig;

  /**
   * Configuration for the linear system solver when solving the final, refined system
   */
  solver::SLESolverConfiguration solverFinalConfig;

  /**
   * Set the type of regularization operation to use
   */
  datadriven::RegularizationConfiguration regularizationConfig;

  /**
   * Influence of regularization term vs data term from 0 (no regularization) to 1 (no data term).
   */
  double lambda;

  /**
   * Determine implementation (openMP, MPI, GPU) that should be used for
   * #sgpp::base::operationMultipleEval
   */
  datadriven::OperationMultipleEvalConfiguration multipleEvalConfig;
};

} /* namespace datadriven */
} /* namespace sgpp */
