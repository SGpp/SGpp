// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#ifdef _OPENMP
#include <omp.h>
#endif

#include <sgpp/globaldef.hpp>
#include <sgpp/optimization/gridgen/IterativeGridGenerator.hpp>

#include <list>
#include <numeric>

namespace sgpp {
namespace optimization {

IterativeGridGenerator::IterativeGridGenerator(base::ScalarFunction& f, base::Grid& grid, size_t N)
    : f(f), grid(grid), N(N), functionValues(0) {}

IterativeGridGenerator::~IterativeGridGenerator() {}

base::Grid& IterativeGridGenerator::getGrid() const { return grid; }

const base::DataVector& IterativeGridGenerator::getFunctionValues() const { return functionValues; }

void IterativeGridGenerator::printIterativeGridGenerator() const {
  base::GridStorage& gridStorage = this->getGrid().getStorage();
  const base::DataVector& functionValues = this->getFunctionValues();

  for (size_t i = 0; i < gridStorage.getSize(); i++) {
    if (i > 0) {
      std::cout << "\n";
    }

    // print grid point and function value
    std::cout << gridStorage[i].toString() << ", " << functionValues[i];
  }

  std::cout << "\n";
}

void IterativeGridGenerator::undoRefinement(size_t oldGridSize) {
  base::GridStorage& gridStorage = grid.getStorage();
  std::list<size_t> indicesToRemove(gridStorage.getSize() - oldGridSize);
  std::iota(indicesToRemove.begin(), indicesToRemove.end(), oldGridSize);
  gridStorage.deletePoints(indicesToRemove);
}

void IterativeGridGenerator::evalFunction(size_t oldGridSize) {
  const size_t d = f.getNumberOfParameters();
  base::GridStorage& gridStorage = grid.getStorage();
  const size_t curGridSize = gridStorage.getSize();
  base::DataVector& fX = functionValues;

#pragma omp parallel shared(fX, oldGridSize, gridStorage)
  {
    base::DataVector x(d);
    base::ScalarFunction* curFPtr = &f;
#ifdef _OPENMP
    std::unique_ptr<base::ScalarFunction> curF;

    if (omp_get_max_threads() > 1) {
      f.clone(curF);
      curFPtr = curF.get();
    }

#endif /* _OPENMP */

#pragma omp for

    for (size_t i = oldGridSize; i < curGridSize; i++) {
      // convert grid point to coordinate vector
      const base::GridPoint& gp = gridStorage[i];

      for (size_t t = 0; t < d; t++) {
        x[t] = gridStorage.getCoordinate(gp, t);
      }

      const double fx = curFPtr->eval(x);
      fX[i] = fx;
    }
  }
}
}  // namespace optimization
}  // namespace sgpp
