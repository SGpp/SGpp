// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <algorithm>

#include <sgpp/globaldef.hpp>

#include <sgpp/optimization/operation/hash/OperationMultipleHierarchisationBsplineBoundary.hpp>
#include <sgpp/base/operation/hash/OperationNaiveEvalBsplineBoundary.hpp>
#include <sgpp/optimization/sle/solver/Auto.hpp>
#include <sgpp/optimization/sle/system/HierarchisationSLE.hpp>

namespace SGPP {
  namespace optimization {

    bool OperationMultipleHierarchisationBsplineBoundary::doHierarchisation(
      base::DataVector& nodeValues) {
      HierarchisationSLE system(grid);
      sle_solver::Auto solver;
      base::DataVector b(nodeValues);
      return solver.solve(system, b, nodeValues);
    }

    void OperationMultipleHierarchisationBsplineBoundary::doDehierarchisation(
      base::DataVector& alpha) {
      base::GridStorage& storage = *grid.getStorage();
      const size_t d = storage.dim();
      base::OperationNaiveEvalBsplineBoundary opNaiveEval(&storage, grid.getDegree());
      base::DataVector nodeValues(storage.size());
      base::DataVector x(d, 0.0);

      for (size_t j = 0; j < storage.size(); j++) {
        const base::GridIndex& gp = *storage.get(j);

        for (size_t t = 0; t < d; t++) {
          x[t] = gp.getCoord(t);
        }

        nodeValues[j] = opNaiveEval.eval(alpha, x);
      }

      alpha.resize(storage.size());
      alpha = nodeValues;
    }

    bool OperationMultipleHierarchisationBsplineBoundary::doHierarchisation(
      std::vector<base::DataVector>& nodeValues) {
      HierarchisationSLE system(grid);
      sle_solver::Auto solver;
      std::vector<base::DataVector> B(nodeValues);
      return solver.solve(system, B, nodeValues);
    }

    void OperationMultipleHierarchisationBsplineBoundary::doDehierarchisation(
      std::vector<base::DataVector>& alpha) {
      base::GridStorage& storage = *grid.getStorage();
      const size_t d = storage.dim();
      base::OperationNaiveEvalBsplineBoundary opNaiveEval(&storage, grid.getDegree());
      base::DataVector nodeValues(storage.size(), 0.0);
      base::DataVector x(d, 0.0);

      for (size_t i = 0; i < storage.size(); i++) {
        for (size_t j = 0; j < storage.size(); j++) {
          const base::GridIndex& gp = storage.get(j);

          for (size_t t = 0; t < d; t++) {
            x[t] = gp.getCoord(t);
          }

          nodeValues[j] = opNaiveEval.eval(alpha[i], x);
        }

        alpha[i].resize(storage.size());
        alpha[i] = nodeValues;
      }
    }

  }
}
