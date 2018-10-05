// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <sgpp/base/datatypes/DataVector.hpp>
#include <sgpp/base/grid/type/NakBsplineGrid.hpp>
#include <sgpp/base/operation/hash/common/basis/BsplineBasis.hpp>
#include <sgpp/base/operation/hash/common/basis/NakBsplineBasis.hpp>
#include <sgpp/base/tools/SGppStopwatch.hpp>
#include <sgpp/optimization/function/scalar/InterpolantScalarFunction.hpp>
#include <sgpp/optimization/sle/solver/Armadillo.hpp>
#include <sgpp/optimization/sle/system/HierarchisationSLE.hpp>

double f(sgpp::base::DataVector v) {
  double res = 1;
  for (auto& a : v) res *= a;
  return res;
}

int main() {
  size_t numDim = 4;
  size_t degree = 3;
  sgpp::base::DataVector evalPoint(numDim, 0.33762);
  for (size_t level = 1; level < 12; level++) {
    sgpp::base::SNakBsplineBase basis(degree);
    auto grid = std::make_shared<sgpp::base::NakBsplineGrid>(numDim, degree);
    sgpp::base::GridStorage& gridStorage = grid->getStorage();
    grid->getGenerator().regular(level);

    sgpp::base::SGppStopwatch watch;
    watch.start();
    for (size_t j = 0; j < gridStorage.getSize(); j++) {
      sgpp::base::GridPoint& gpBasis = gridStorage.getPoint(j);
      double reducedBasisEval = 1;
      for (size_t t = 0; t < numDim; t++) {
        double basisEval1D = basis.eval(gpBasis.getLevel(t), gpBasis.getIndex(t), evalPoint[t]);
        reducedBasisEval *= basisEval1D;
      }
    }
    double runtime = watch.stop();
    std::cout << level << "  total: " << runtime << " numPoints: " << gridStorage.getSize()
              << " per point: " << runtime / static_cast<double>(gridStorage.getSize())
              << std::endl;
  }

  for (size_t level = 1; level < 7; level++) {
    auto grid = std::make_shared<sgpp::base::NakBsplineGrid>(numDim, degree);
    sgpp::base::GridStorage& gridStorage = grid->getStorage();
    grid->getGenerator().regular(level);
    sgpp::base::DataVector functionValues(gridStorage.getSize());
    for (size_t i = 0; i < gridStorage.getSize(); i++) {
      sgpp::base::GridPoint& gp = gridStorage.getPoint(i);
      sgpp::base::DataVector gridPointVector(gridStorage.getDimension());
      gp.getStandardCoordinates(gridPointVector);
      functionValues[i] = f(gridPointVector);
    }

    // solve linear system
    sgpp::base::DataVector alpha(functionValues.getSize());
    sgpp::optimization::HierarchisationSLE hierSLE(*grid);
    sgpp::optimization::sle_solver::Armadillo sleSolver;
    if (!sleSolver.solve(hierSLE, functionValues, alpha)) {
      std::cout << "ASMatrixNakBspline: Solving failed.\n";
      return 0;
    }

    sgpp::optimization::InterpolantScalarFunction I(*grid, alpha);
    sgpp::base::SGppStopwatch watch;
    watch.start();
    I.eval(evalPoint);
    double runtime = watch.stop();
    std::cout << level << "  total: " << runtime << " numPoints: " << gridStorage.getSize()
              << " per point: " << runtime / static_cast<double>(gridStorage.getSize())
              << std::endl;
  }
  return 0;
}
