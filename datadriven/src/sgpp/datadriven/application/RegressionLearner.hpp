// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#ifndef REGRESSIONLEARNER_H
#define REGRESSIONLEARNER_H

#include <sgpp/base/datatypes/DataVector.hpp>
#include <sgpp/base/datatypes/DataMatrix.hpp>
#include <sgpp/base/operation/hash/OperationMatrix.hpp>
#include <sgpp/datadriven/algorithm/DMSystemMatrixBase.hpp>
#include <sgpp/base/grid/Grid.hpp>
#include <sgpp/base/exception/application_exception.hpp>
#include <sgpp/datadriven/application/RegularizationConfiguration.hpp>
#include <sgpp/solver/TypesSolver.hpp>
#include <sgpp/solver/SLESolver.hpp>
#include <sgpp/solver/sle/fista/FistaBase.hpp>
#include <sgpp/globaldef.hpp>

#include <memory>

namespace sgpp {
namespace datadriven {

class RegressionLearner {
 public:
  // Fista and (Bi)-CG do not share a common interace.
  // This is why we need this slightly dirty solution.
  class Solver {
   public:
    enum class solverCategory { cg, fista, none } type = solverCategory::none;
    Solver() {}
    explicit Solver(std::unique_ptr<sgpp::solver::SLESolver>&& s) {
      new (&solverCG) std::unique_ptr<sgpp::solver::SLESolver>{std::move(s)};
      type = solverCategory::cg;
    }
    explicit Solver(std::unique_ptr<sgpp::solver::FistaBase>&& s) {
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
                    const sgpp::base::DataVector& b, size_t maxIt, double treshold) {
      if (type != solverCategory::fista) {
        throw sgpp::base::application_exception("Tried to solve with incorrect solver!");
      }
      solverFista->solve(op, weights, b, maxIt, treshold);
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
        default:
          // Do nothing.
          break;
      }
    }

    Solver(const Solver&) = delete;
    Solver(Solver&& other) : Solver() { swap(*this, other); }
    Solver& operator=(const Solver&) = delete;
    Solver& operator=(Solver&& other) {
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
        default:
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

  RegressionLearner(sgpp::base::RegularGridConfiguration gridConfig,
                    sgpp::base::AdpativityConfiguration adaptivityConfig,
                    sgpp::solver::SLESolverConfiguration solverConfig,
                    sgpp::solver::SLESolverConfiguration finalSolverConfig,
                    datadriven::RegularizationConfiguration regularizationConfig);
  void train(sgpp::base::DataMatrix& trainDataset, sgpp::base::DataVector& classes);
  sgpp::base::DataVector predict(sgpp::base::DataMatrix& data);
  size_t getGridSize() const;
  double getMSE(sgpp::base::DataMatrix& data, const sgpp::base::DataVector& y);
  void initializeWeights();
  sgpp::base::DataVector getWeights() const;

 private:
  sgpp::base::RegularGridConfiguration gridConfig;
  sgpp::base::AdpativityConfiguration adaptivityConfig;
  sgpp::solver::SLESolverConfiguration solverConfig;
  sgpp::solver::SLESolverConfiguration finalSolverConfig;
  datadriven::RegularizationConfiguration regularizationConfig;
  std::unique_ptr<sgpp::base::OperationMultipleEval> op;
  std::unique_ptr<datadriven::DMSystemMatrixBase> systemMatrix;

  /// sparse grid object
  std::unique_ptr<sgpp::base::Grid> grid;
  /// the grid's coefficients
  sgpp::base::DataVector weights;

  void initializeGrid(sgpp::base::RegularGridConfiguration GridConfig);
  std::unique_ptr<datadriven::DMSystemMatrixBase> createDMSystem(
      sgpp::base::DataMatrix& trainDataset);
  Solver createSolver();
  Solver createSolverFista();

  void fit(Solver& solver, sgpp::base::DataVector& classes);
  void refine(sgpp::base::DataMatrix& data, sgpp::base::DataVector& classes);

  double getMSE(const sgpp::base::DataVector& y, sgpp::base::DataVector yPrediction);
};

}  // namespace datadriven
}  // namespace sgpp

#endif  // REGRESSIONLEARNER_H
