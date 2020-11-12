// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#pragma once

#include <sgpp/globaldef.hpp>
#include <sgpp/base/datatypes/DataMatrix.hpp>
#include <sgpp/base/datatypes/DataVector.hpp>
#include <sgpp/base/exception/application_exception.hpp>
#include <sgpp/base/grid/Grid.hpp>
#include <sgpp/base/operation/hash/OperationMatrix.hpp>
#include <sgpp/datadriven/algorithm/DMSystemMatrixBase.hpp>
#include <sgpp/datadriven/configuration/RegularizationConfiguration.hpp>
#include <sgpp/solver/SLESolver.hpp>
#include <sgpp/solver/TypesSolver.hpp>
#include <sgpp/solver/sle/fista/FistaBase.hpp>

#include <algorithm>
#include <memory>
#include <set>
#include <vector>
#include <utility>

namespace sgpp {
namespace datadriven {
/**
 * @brief The RegressionLearner class
 * Solves a regression problem with continuous target vector.
 */
class RegressionLearner {
 public:
// Fista and (Bi)-CG do not share a common interace.
// This is why we need this slightly dirty solution.
// Excluded from SWIG.
#ifndef SWIG
  class Solver {
   public:
    enum class solverCategory { cg, fista, none } type = solverCategory::none;
    Solver() {}
    explicit Solver(std::unique_ptr<sgpp::solver::SLESolver>&& s) {  // NOLINT(build/c++11)
      new (&solverCG) std::unique_ptr<sgpp::solver::SLESolver>{std::move(s)};
      type = solverCategory::cg;
    }
    explicit Solver(std::unique_ptr<sgpp::solver::FistaBase>&& s) {  // NOLINT(build/c++11)
      new (&solverFista) std::unique_ptr<sgpp::solver::FistaBase>{std::move(s)};
      type = solverCategory::fista;
    }
    void solveCG(sgpp::base::OperationMatrix& systemMatrix, sgpp::base::DataVector& alpha,
                 sgpp::base::DataVector& b, bool reuse, bool verbose, double maxTreshold) {
      if (type != solverCategory::cg) {
        throw sgpp::base::application_exception("Tried to solve with incorrect solver!");
      }
      solverCG->solve(systemMatrix, alpha, b, reuse, verbose, maxTreshold);
    }
    void solveFista(sgpp::base::OperationMultipleEval& op, sgpp::base::DataVector& weights,
                    const sgpp::base::DataVector& classes, size_t maxIt, double treshold,
                    double L) {
      if (type != solverCategory::fista) {
        throw sgpp::base::application_exception("Tried to solve with incorrect solver!");
      }
      solverFista->solve(op, weights, classes, maxIt, treshold, L);
    }
    double getL() {
      if (type != solverCategory::fista) {
        throw sgpp::base::application_exception("Solver doesn't support L!");
      }
      return solverFista->getL();
    }

    friend void swap(Solver& first, Solver& second) {
      using std::swap;
      swap(first.type, second.type);
      switch (first.type) {
        case solverCategory::cg:
          first.solverCG.swap(second.solverCG);
          break;
        case solverCategory::fista:
          first.solverFista.swap(second.solverFista);
          break;
        case solverCategory::none:
          // Do nothing.
          break;
      }
    }

    Solver(const Solver&) = delete;
    Solver(Solver&& other) : Solver() { swap(*this, other); }  // NOLINT(build/c++11)
    Solver& operator=(const Solver&) = delete;
    Solver& operator=(Solver&& other) {  // NOLINT(build/c++11)
      swap(*this, other);
      return *this;
    }

    ~Solver() {
      switch (type) {
        case solverCategory::cg:
          solverCG.~unique_ptr<sgpp::solver::SLESolver>();
          break;
        case solverCategory::fista:
          solverFista.~unique_ptr<sgpp::solver::FistaBase>();
          break;
        case solverCategory::none:
          // do nothing
          break;
      }
    }

   private:
    union {
      std::unique_ptr<sgpp::solver::SLESolver> solverCG;
      std::unique_ptr<sgpp::solver::FistaBase> solverFista;
    };
  };
#endif  // end inner class

  /**
   * @brief RegressionLearner
   * @param gridConfig
   * @param adaptivityConfig
   * @param solverConfig is the solver used during each adaptivity step
   * @param finalSolverConfig is the solver used to build the final model
   * @param regularizationConfig
   * @param terms is a vector that contains all desired interaction terms.
   * For example, if we want to include grid points that model an
   * interaction between the first and the second predictor, we would
   * include the vector [1,2] in terms.
   */
  RegressionLearner(sgpp::base::RegularGridConfiguration gridConfig,
                    sgpp::base::AdaptivityConfiguration adaptivityConfig,
                    sgpp::solver::SLESolverConfiguration solverConfig,
                    sgpp::solver::SLESolverConfiguration finalSolverConfig,
                    datadriven::RegularizationConfiguration regularizationConfig,
                    std::set<std::set<size_t>> terms);

  /**
   * @brief RegressionLearner
   * @param gridConfig
   * @param adaptivityConfig
   * @param solverConfig is the solver used during each adaptivity step
   * @param finalSolverConfig is the solver used to build the final model
   * @param regularizationConfig
   */
  RegressionLearner(sgpp::base::RegularGridConfiguration gridConfig,
                    sgpp::base::AdaptivityConfiguration adaptivityConfig,
                    sgpp::solver::SLESolverConfiguration solverConfig,
                    sgpp::solver::SLESolverConfiguration finalSolverConfig,
                    datadriven::RegularizationConfiguration regularizationConfig);
  /**
   * @brief train fits a sparse grid regression model.
   * @param trainDataset is the design matrix
   * @param classes is the (continuous) target
   */
  void train(sgpp::base::DataMatrix& trainDataset, sgpp::base::DataVector& classes);
  /**
   * @brief predict
   * @param data are observations
   * @return the predicted target for matrix data
   */
  sgpp::base::DataVector predict(sgpp::base::DataMatrix& data);
  /**
   * @brief getGridSize
   * @return the size of the grid
   */
  size_t getGridSize() const;
  /**
   * @brief getGrid
   * @return the grid
   */
  sgpp::base::Grid& getGrid();
  /**
   * @brief getMSE
   * @param data is the design matrix
   * @param y is the target
   * @return the mean-squared-error of the prediction of the model for the matrix data
   */
  double getMSE(sgpp::base::DataMatrix& data, const sgpp::base::DataVector& y);
  /**
   * @brief getWeights
   * @return the weights
   */
  sgpp::base::DataVector getWeights() const;
  /**
   * @brief setWeights
   * @param weights are the new weights.
   */
  void setWeights(sgpp::base::DataVector weights);

 private:
  sgpp::base::RegularGridConfiguration gridConfig;
  sgpp::base::AdaptivityConfiguration adaptivityConfig;
  sgpp::solver::SLESolverConfiguration solverConfig;
  sgpp::solver::SLESolverConfiguration finalSolverConfig;
  datadriven::RegularizationConfiguration regularizationConfig;
  std::set<std::set<size_t>> terms;
  std::unique_ptr<sgpp::base::OperationMultipleEval> op;
  std::unique_ptr<datadriven::DMSystemMatrixBase> systemMatrix;

  /// sparse grid object
  std::unique_ptr<sgpp::base::Grid> grid;
  /// the grid's coefficients
  sgpp::base::DataVector weights;

  void initializeGrid(sgpp::base::RegularGridConfiguration gridConfig);
  std::unique_ptr<datadriven::DMSystemMatrixBase> createDMSystem(
      sgpp::base::DataMatrix& trainDataset);
  Solver createSolver(size_t n_rows);
  Solver createSolverFista(size_t n_rows);

  void fit(Solver& solver, sgpp::base::DataVector& classes);
  void refine(sgpp::base::DataMatrix& data, sgpp::base::DataVector& classes);

  double getMSE(const sgpp::base::DataVector& y, sgpp::base::DataVector yPrediction);
};

}  // namespace datadriven
}  // namespace sgpp
