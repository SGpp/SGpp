// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

//#include <algorithm>
#include <queue>
#include <stack>

#include <sgpp/globaldef.hpp>

#include <sgpp/optimization/operation/hash/OperationMultipleHierarchisationModFundamentalSpline.hpp>
#include <sgpp/base/operation/hash/OperationNaiveEvalModFundamentalSpline.hpp>
//#include <sgpp/optimization/sle/solver/Auto.hpp>
//#include <sgpp/optimization/sle/system/HierarchisationSLE.hpp>
#include <sgpp/base/operation/hash/common/basis/FundamentalSplineModifiedBasis.hpp>

namespace SGPP {
  namespace optimization {

    bool OperationMultipleHierarchisationModFundamentalSpline::doBFS(
      const std::vector<base::DataVector>& nodeValues,
      std::vector<base::DataVector>& alpha) {
      base::GridStorage& gridStorage = *grid.getStorage();
      const size_t n = grid.getSize();
      const size_t d = gridStorage.dim();
      alpha.clear();

      for (size_t m = 0; m < nodeValues.size(); m++) {
        alpha.push_back(nodeValues[m]);
      }

      std::vector<bool> bfsVisited(n, false);
      base::GridIndex bfsStart(d);
      //base::HashGridIterator iterator(grid.getStorage());

      for (size_t t = 0; t < d; t++) {
        bfsStart.push(t, 1, 1);
      }

      bfsStart.rehash();
      std::queue<size_t> queue;
      const size_t bfsStartIndex = gridStorage.seq(&bfsStart);

      if (bfsStartIndex >= n) {
        return false;
      }

      queue.push(bfsStartIndex);
      bfsVisited[bfsStartIndex] = true;

      base::SFundamentalSplineModifiedBase base(grid.getDegree());

      while (!queue.empty()) {
        size_t bfsPointIndex = queue.front();
        queue.pop();
        base::GridIndex& bfsPoint = *gridStorage.get(bfsPointIndex);

        {
          std::vector<bool> dfsVisited(n, false);
          std::stack<size_t> stack;
          stack.push(bfsPointIndex);

          while (!stack.empty()) {
            size_t dfsPointIndex = stack.top();
            stack.pop();
            base::GridIndex& dfsPoint = *gridStorage.get(dfsPointIndex);

            if (!dfsVisited[dfsPointIndex]) {
              dfsVisited[dfsPointIndex] = true;

              if (dfsPointIndex != bfsPointIndex) {
                float_t value = 1.0;

                for (size_t t = 0; t < d; t++) {
                  const float_t val1d = base.eval(bfsPoint.getLevel(t),
                                                  bfsPoint.getIndex(t),
                                                  dfsPoint.getCoord(t));

                  if (val1d == 0.0) {
                    value = 0.0;
                    break;
                  }

                  value *= val1d;
                }

                for (size_t m = 0; m < alpha.size(); m++) {
                  alpha[m][dfsPointIndex] -= alpha[m][bfsPointIndex] * value;
                }
              }

              for (size_t t = 0; t < d; t++) {
                base::GridIndex::level_type l;
                base::GridIndex::index_type i;
                dfsPoint.get(t, l, i);

                {
                  dfsPoint.set(t, l + 1, 2 * i - 1);
                  dfsPointIndex = gridStorage.seq(&dfsPoint);

                  if (dfsPointIndex < n) {
                    stack.push(dfsPointIndex);
                  }
                }

                {
                  dfsPoint.set(t, l + 1, 2 * i + 1);
                  dfsPointIndex = gridStorage.seq(&dfsPoint);

                  if (dfsPointIndex < n) {
                    stack.push(dfsPointIndex);
                  }
                }

                dfsPoint.push(t, l, i);
              }
            }
          }
        }

        for (size_t t = 0; t < d; t++) {
          base::GridIndex::level_type l;
          base::GridIndex::index_type i;
          bfsPoint.get(t, l, i);

          {
            bfsPoint.set(t, l + 1, 2 * i - 1);
            bfsPointIndex = gridStorage.seq(&bfsPoint);

            if ((bfsPointIndex < n) && (!bfsVisited[bfsPointIndex])) {
              queue.push(bfsPointIndex);
              bfsVisited[bfsPointIndex] = true;
            }
          }

          {
            bfsPoint.set(t, l + 1, 2 * i + 1);
            bfsPointIndex = gridStorage.seq(&bfsPoint);

            if ((bfsPointIndex < n) && (!bfsVisited[bfsPointIndex])) {
              queue.push(bfsPointIndex);
              bfsVisited[bfsPointIndex] = true;
            }
          }

          bfsPoint.push(t, l, i);
        }
      }

      for (size_t k = 0; k < n; k++) {
        if (!bfsVisited[k]) {
          return false;
        }
      }

      return true;
    }

    bool OperationMultipleHierarchisationModFundamentalSpline::doHierarchisation(
      base::DataVector& nodeValues) {
      /*HierarchisationSLE system(grid);
      sle_solver::Auto solver;
      base::DataVector b(nodeValues);
      return solver.solve(system, b, nodeValues);*/

      std::vector<base::DataVector> nodeValuesVec{nodeValues};
      std::vector<base::DataVector> resultVec;

      if (!doBFS(nodeValuesVec, resultVec)) {
        return false;
      }

      nodeValues = resultVec[0];
      return true;
    }

    void OperationMultipleHierarchisationModFundamentalSpline::doDehierarchisation(
      base::DataVector& alpha) {
      base::GridStorage& storage = *grid.getStorage();
      const size_t d = storage.dim();
      base::OperationNaiveEvalModFundamentalSpline opNaiveEval(&storage, grid.getDegree());
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

    bool OperationMultipleHierarchisationModFundamentalSpline::doHierarchisation(
      std::vector<base::DataVector>& nodeValues) {
      /*HierarchisationSLE system(grid);
      sle_solver::Auto solver;
      std::vector<base::DataVector> B(nodeValues);
      return solver.solve(system, B, nodeValues);*/

      std::vector<base::DataVector> result;

      if (!doBFS(nodeValues, result)) {
        return false;
      }

      nodeValues = result;
      return true;
    }

    void OperationMultipleHierarchisationModFundamentalSpline::doDehierarchisation(
      std::vector<base::DataVector>& alpha) {
      base::GridStorage& storage = *grid.getStorage();
      const size_t d = storage.dim();
      base::OperationNaiveEvalModFundamentalSpline opNaiveEval(&storage, grid.getDegree());
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
