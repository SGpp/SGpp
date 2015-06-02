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

      std::vector<bool> visited(n, false);
      base::GridIndex startPoint(d);

      for (size_t t = 0; t < d; t++) {
        startPoint.push(t, 1, 1);
      }

      startPoint.rehash();
      std::queue<size_t> queue;
      size_t pointIndex = gridStorage.seq(&startPoint);

      if (pointIndex >= n) {
        return false;
      }

      queue.push(pointIndex);
      visited[pointIndex] = true;

      base::SFundamentalSplineModifiedBase base(grid.getDegree());

      /*size_t counter1 = 0, counter2 = 0;
      size_t counter3 = 0, counter4 = 0;*/

      while (!queue.empty()) {
        pointIndex = queue.front();
        queue.pop();

        base::GridIndex& point = *gridStorage[pointIndex];

        for (size_t q = 0; q < n; q++) {
          base::GridIndex& childPoint = *gridStorage[q];
          bool skipChild = false;

          if (q == pointIndex) {
            continue;
          }

          for (size_t t = 0; t < d; t++) {
            const base::GridIndex::level_type l = point.getLevel(t);
            const base::GridIndex::level_type k = childPoint.getLevel(t);

            if ((k <= l) && ((k != l) ||
                             (point.getIndex(t) != childPoint.getIndex(t)))) {
              skipChild = true;
              break;
            }
          }

          /*TODO
          if (skipChild) {
            counter4++;
          } else {
            counter3++;
          }*/

          if (!skipChild) {
            float_t value = 1.0;

            for (size_t t = 0; t < d; t++) {
              const float_t val1d = base.eval(point.getLevel(t),
                                              point.getIndex(t),
                                              childPoint.getCoord(t));

              if (val1d == 0.0) {
                value = 0.0;
                break;
              }

              value *= val1d;
            }

            if (value != 0.0) {
              //TODO
              //counter1++;

              for (size_t m = 0; m < alpha.size(); m++) {
                alpha[m][q] -= alpha[m][pointIndex] * value;
              }
            } /* else { //TODO

              counter2++;
            }*/
          }
        }

        for (size_t t = 0; t < d; t++) {
          const base::GridIndex::level_type l = point.getLevel(t);
          const base::GridIndex::index_type i = point.getIndex(t);

          {
            point.set(t, l + 1, 2 * i - 1);
            pointIndex = gridStorage.seq(&point);

            if ((pointIndex < n) && (!visited[pointIndex])) {
              queue.push(pointIndex);
              visited[pointIndex] = true;
            }
          }

          {
            point.set(t, l + 1, 2 * i + 1);
            pointIndex = gridStorage.seq(&point);

            if ((pointIndex < n) && (!visited[pointIndex])) {
              queue.push(pointIndex);
              visited[pointIndex] = true;
            }
          }

          point.push(t, l, i);
        }
      }

      for (size_t q = 0; q < n; q++) {
        if (!visited[q]) {
          return false;
        }
      }

      /*TODO
      std::cout << "counter1 = " << counter1 << "\n";
      std::cout << "counter2 = " << counter2 << "\n";
      std::cout << "counter3 = " << counter3 << "\n";
      std::cout << "counter4 = " << counter4 << "\n";*/

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
