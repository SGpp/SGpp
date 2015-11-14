// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <algorithm>

#include <sgpp/globaldef.hpp>

#include <sgpp/optimization/operation/hash/OperationMultipleHierarchisationWaveletBoundary.hpp>
#include <sgpp/base/operation/hash/OperationNaiveEvalWaveletBoundary.hpp>
#include <sgpp/optimization/sle/solver/Auto.hpp>
#include <sgpp/optimization/sle/system/HierarchisationSLE.hpp>

namespace SGPP {
  namespace optimization {

    OperationMultipleHierarchisationWaveletBoundary::OperationMultipleHierarchisationWaveletBoundary(
      base::WaveletBoundaryGrid& grid) :
      grid(grid) {
    }

    OperationMultipleHierarchisationWaveletBoundary::~OperationMultipleHierarchisationWaveletBoundary() {
    }

    bool OperationMultipleHierarchisationWaveletBoundary::doHierarchisation(
      base::DataVector& nodeValues) {
      HierarchisationSLE system(grid);
      sle_solver::Auto solver;
      base::DataVector b(nodeValues);
      return solver.solve(system, b, nodeValues);
    }

    void OperationMultipleHierarchisationWaveletBoundary::doDehierarchisation(
      base::DataVector& alpha) {
      base::GridStorage& storage = *grid.getStorage();
      const size_t d = storage.dim();
      base::OperationNaiveEvalWaveletBoundary opNaiveEval(&storage);
      base::DataVector nodeValues(storage.size());
      base::DataVector x(d, 0.0);

      for (size_t j = 0; j < storage.size(); j++) {
        const base::GridIndex& gp = *storage[j];

        for (size_t t = 0; t < d; t++) {
          x[t] = gp.getCoord(t);
        }

        nodeValues[j] = opNaiveEval.eval(alpha, x);
      }

      alpha.resize(storage.size());
      alpha = nodeValues;
    }

    bool OperationMultipleHierarchisationWaveletBoundary::doHierarchisation(
      base::DataMatrix& nodeValues) {
      HierarchisationSLE system(grid);
      sle_solver::Auto solver;
      base::DataMatrix B(nodeValues);
      return solver.solve(system, B, nodeValues);
    }

    void OperationMultipleHierarchisationWaveletBoundary::doDehierarchisation(
      base::DataMatrix& alpha) {
      base::GridStorage& storage = *grid.getStorage();
      const size_t d = storage.dim();
      base::OperationNaiveEvalWaveletBoundary opNaiveEval(&storage);
      base::DataVector nodeValues(storage.size(), 0.0);
      base::DataVector x(d, 0.0);
      base::DataVector alpha1(storage.size(), 0.0);

      for (size_t i = 0; i < alpha.getNcols(); i++) {
        alpha.getColumn(i, alpha1);

        for (size_t j = 0; j < storage.size(); j++) {
          const base::GridIndex& gp = *storage[j];

          for (size_t t = 0; t < d; t++) {
            x[t] = gp.getCoord(t);
          }

          nodeValues[j] = opNaiveEval.eval(alpha1, x);
        }

        alpha.setColumn(i, nodeValues);
      }
    }

  }
}
