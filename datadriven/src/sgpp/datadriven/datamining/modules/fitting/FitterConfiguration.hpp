// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#pragma once

#include <sgpp/base/grid/Grid.hpp>
#include <sgpp/datadriven/configuration/CrossvalidationConfiguration.hpp>
#include <sgpp/datadriven/configuration/DatabaseConfiguration.hpp>
#include <sgpp/datadriven/configuration/DensityEstimationConfiguration.hpp>

#include <sgpp/datadriven/configuration/GeometryConfiguration.hpp>
#include <sgpp/datadriven/configuration/LearnerConfiguration.hpp>
#include <sgpp/datadriven/configuration/ParallelConfiguration.hpp>

#include <sgpp/datadriven/configuration/RegularizationConfiguration.hpp>
#include <sgpp/datadriven/datamining/configuration/DataMiningConfigParser.hpp>
#include <sgpp/datadriven/operation/hash/DatadrivenOperationCommon.hpp>
#include <sgpp/solver/TypesSolver.hpp>

#include <string>

namespace sgpp {
namespace datadriven {
// forward declaration to break dependency cycle
class DataMiningConfigParser;

/**
 * Different fitter scenarios have different default values and support different operations
 */
enum class FitterType {
  RegressionLeastSquares,
  DensityEstimation,
  DensityRatioEstimation,
  DensityDifferenceEstimation,
  Classification
};

/**
 * General configuration object for fitters. Bundles all structures needed to build a sparse grid,
 * fit a sparse grid based model, and perform adaptive refinement.
 */
class FitterConfiguration {
 public:
  /**
   * Sets up a Fitter configuration with its default values.
   */
  FitterConfiguration() = default;

  /**
   * Copy constructor
   * @param rhs const reference to the scorer object to copy from.
   */
  FitterConfiguration(const FitterConfiguration &rhs) = default;

  /**
   * Move constructor
   * @param rhs R-value reference to a scorer object to moved from.
   */
  FitterConfiguration(FitterConfiguration &&rhs) = default;

  /**
   * Copy assign operator
   * @param rhs const reference to the scorer object to copy from.
   * @return rerefernce to this with updated values.
   */
  FitterConfiguration &operator=(const FitterConfiguration &rhs) = default;

  /**
   * Move assign operator
   * @param rhs R-value reference to an a scorer object to move from.
   * @return rerefernce to this with updated values.
   */
  FitterConfiguration &operator=(FitterConfiguration &&rhs) = default;

  /**
   * virtual destructor.
   */
  virtual ~FitterConfiguration() = default;

  /**
   * Polymorphic clone pattern
   * @return deep copy of this object. New object is owned by caller.
   */
  virtual FitterConfiguration *clone() const = 0;

  /**
   * Get initial conditions for the grid before adaptive refinement.
   * @return immutable GeneralGridConfiguration
   */
  const base::GeneralGridConfiguration &getGridConfig() const;

  /**
   * Get how the adaptivity algorithms for the grid should behave.
   * @return immutable AdaptivityConfiguration
   */
  const base::AdaptivityConfiguration &getRefinementConfig() const;

  /**
   * Get how the crossvalidation should behave.
   * @return immutable CrossvalidationConfiguration
   */
  const datadriven::CrossvalidationConfiguration &getCrossvalidationConfig() const;

  /**
   * Get how the density estimation should behave.
   * @return immutable DensityEstimationConfiguration
   */
  const datadriven::DensityEstimationConfiguration &getDensityEstimationConfig() const;

  /**
   * Get configuration for the linear system solver which should be used while building adaptive
   * grids
   * @return immutable SLESolverConfiguration
   */
  const solver::SLESolverConfiguration &getSolverRefineConfig() const;

  /**
   * Get configuration for the linear system solver when solving the final, refined system
   * @return immutable SLESolverConfiguration
   */
  const solver::SLESolverConfiguration &getSolverFinalConfig() const;

  /**
   * Get the type of regularization operation to use
   * @return immutable RegularizationConfiguration
   */
  const datadriven::RegularizationConfiguration &getRegularizationConfig() const;

  /**
   * Get implementation (openMP, MPI, GPU) that should be used for
   * #sgpp::base::OperationMultipleEval.
   * @return immutable OperationMultipleEvalConfiguration
   */
  const datadriven::OperationMultipleEvalConfiguration &getMultipleEvalConfig() const;

  /**
   * Returns the database configuration, i.e. the filepath
   * @return immutable DatabaseConfiguration
   */
  const datadriven::DatabaseConfiguration &getDatabaseConfig() const;

  /**
   * Returns the configuration for the learner's behaviour
   * @return immutable LearnerConfiguration
   */
  const datadriven::LearnerConfiguration &getLearnerConfig() const;

  /**
   * Returns the configuration for parallelization with ScaLAPACK
   * @return immutable ParallelConfiguration
   */
  const datadriven::ParallelConfiguration &getParallelConfig() const;

  /*
   * Returns the configuration for the geometry parameters
   * @return immutable GeometryConfiguration
   */
  const datadriven::GeometryConfiguration &getGeometryConfig() const;

  /**
   * Get or set initial conditions for the grid before adaptive refinement.
   * @return GeneralGridConfiguration
   */
  base::GeneralGridConfiguration &getGridConfig();

  /**
   * Get or set how the adaptivity algorithms for the grid should behave.
   * @return AdaptivityConfiguration
   */
  base::AdaptivityConfiguration &getRefinementConfig();

  /**
   * Get or set how the crossvalidation should behave.
   * @return CrossvalidationConfiguration
   */
  datadriven::CrossvalidationConfiguration &getCrossvalidationConfig();

  /**
   * Get or set how the density estimation should behave.
   * @return DensityEstimationConfiguration
   */
  datadriven::DensityEstimationConfiguration &getDensityEstimationConfig();

  /**
   * Get or set configuration for the linear system solver which should be used while building
   * adaptive grids
   * @return SLESolverConfiguration
   */
  solver::SLESolverConfiguration &getSolverRefineConfig();

  /**
   * Get or set configuration for the linear system solver when solving the final, refined system
   * @return SLESolverConfiguration
   */
  solver::SLESolverConfiguration &getSolverFinalConfig();

  /**
   * Get or set the type of regularization operation to use
   * @return RegularizationConfiguration
   */
  datadriven::RegularizationConfiguration &getRegularizationConfig();

  /**
   * Get or set implementation (openMP, MPI, GPU) that should be used for
   * #sgpp::base::OperationMultipleEval.
   * @return current OperationMultipleEvalConfiguration
   */
  datadriven::OperationMultipleEvalConfiguration &getMultipleEvalConfig();

  /**
   * set default values for all members based on the desired scenario.
   */
  virtual void setupDefaults();

  /**
   * obtain parameters from a parser
   * @param parser: the parser object to read from
   */
  virtual void readParams(const DataMiningConfigParser &parser) = 0;

  /**
   * print out all parameters to stream
   */
  void dumpToStream(std::ostream &stream_out = std::cout) const;

 protected:
  /**
   * Initial conditions for the grid before adaptive refinement.
   */
  base::GeneralGridConfiguration gridConfig;

  /**
   * Configure how the adaptivity algorithms for the grid should behave.
   */
  base::AdaptivityConfiguration adaptivityConfig;

  /**
   * Configure how the crossvalidation should behave.
   */
  datadriven::CrossvalidationConfiguration crossvalidationConfig;

  /**
   * Configure how the density estimation should behave.
   */
  datadriven::DensityEstimationConfiguration densityEstimationConfig;

  /**
   * Configure where the lhs datamatrix decomposition database is stored
   */
  datadriven::DatabaseConfiguration databaseConfig;

  /**
   * Configuration for the linear system solver which should be used while building adaptive grids
   */
  solver::SLESolverConfiguration solverRefineConfig;

  /**
   * Configuration for the linear system solver when solving the final, refined system
   */
  solver::SLESolverConfiguration solverFinalConfig;

  /**
   * Set the type of regularization operation to use and specify the influence of the regularization
   * term vs data term from 0 (no regularization) to 1 (no data term).
   */
  datadriven::RegularizationConfiguration regularizationConfig;

  /**
   * Determine implementation (openMP, MPI, GPU) that should be used for
   * #sgpp::base::OperationMultipleEval
   */
  datadriven::OperationMultipleEvalConfiguration multipleEvalConfig;

  /**
   * Configuration for the learner's behaviour
   */
  datadriven::LearnerConfiguration learnerConfig;

  /*
   * Configuration of the geometry parameters
   */
  datadriven::GeometryConfiguration geometryConfig;

  /**
   *  Configuration for parallelization with ScaLAPACK
   */
  datadriven::ParallelConfiguration parallelConfig;
};
} /* namespace datadriven */
} /* namespace sgpp */
