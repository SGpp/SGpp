// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <algorithm>

#include <sgpp/globaldef.hpp>

#include <sgpp/optimization/operation/hash/OperationMultipleHierarchisationLinearClenshawCurtis.hpp>
#include <sgpp/base/operation/hash/OperationNaiveEvalLinearClenshawCurtis.hpp>
#include <sgpp/optimization/sle/system/Hierarchisation.hpp>
#include <sgpp/optimization/sle/solver/Auto.hpp>

namespace SGPP {
namespace optimization {

void OperationMultipleHierarchisationLinearClenshawCurtis::doHierarchisation(
    std::vector<base::DataVector*> nodeValues) {
  sle::system::Hierarchisation system(grid);
  sle::solver::Auto solver;
  std::vector<std::vector<float_t> > B;
  std::vector<std::vector<float_t> > X;

  for (size_t i = 0; i < nodeValues.size(); i++) {
    B.push_back(
        std::vector<float_t>(
            nodeValues[i]->getPointer(),
            nodeValues[i]->getPointer() + nodeValues[i]->getSize()));
  }

  if (solver.solve(system, B, X)) {
    for (size_t i = 0; i < nodeValues.size(); i++) {
      std::copy(X[i].begin(), X[i].begin() + nodeValues[i]->getSize(),
                nodeValues[i]->getPointer());
    }
  }
}

void OperationMultipleHierarchisationLinearClenshawCurtis::doDehierarchisation(
    std::vector<base::DataVector*> alpha) {
  base::GridStorage* storage = grid.getStorage();
  const size_t d = storage->dim();
  base::OperationNaiveEvalLinearClenshawCurtis op_eval(storage);
  std::vector<float_t> nodeValuesVector(storage->size(), 0.0);
  std::vector<float_t> x(d, 0.0);

  for (size_t i = 0; i < storage->size(); i++) {
    for (size_t j = 0; j < storage->size(); j++) {
      base::GridIndex* gp = storage->get(j);

      for (size_t t = 0; t < d; t++) {
        x[t] = gp->getCoord(t);
      }

      nodeValuesVector[j] = op_eval.eval(*alpha[i], x);
    }

    std::copy(nodeValuesVector.begin(),
              nodeValuesVector.begin() + alpha[i]->getSize(),
              alpha[i]->getPointer());
  }
}

}
}
