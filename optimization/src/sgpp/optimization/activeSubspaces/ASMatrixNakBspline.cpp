// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

//#ifdef USE_EIGEN

#include <sgpp/optimization/activeSubspaces/ASMatrixNakBspline.hpp>

namespace sgpp {
namespace optimization {
void ASMatrixNakBspline::buildRegularInterpolant(size_t level) {
  size_t numDim = objectiveFunc.getNumberOfParameters();
  auto regularGrid = std::make_shared<sgpp::base::NakBsplineGrid>(numDim, degree);
  sgpp::base::GridStorage& gridStorage = regularGrid->getStorage();
  regularGrid->getGenerator().regular(level);
  sgpp::base::DataVector functionValues(gridStorage.getSize());
  for (size_t i = 0; i < gridStorage.getSize(); i++) {
    sgpp::base::GridPoint& gp = gridStorage.getPoint(i);
    sgpp::base::DataVector gridPointVector(gridStorage.getDimension());
    gp.getStandardCoordinates(gridPointVector);
    functionValues[i] = objectiveFunc.eval(gridPointVector);
  }

  // solve linear system
  sgpp::base::DataVector alpha(functionValues.getSize());
  sgpp::optimization::HierarchisationSLE hierSLE(*regularGrid);
  sgpp::optimization::sle_solver::Armadillo sleSolver;
  if (!sleSolver.solve(hierSLE, functionValues, alpha)) {
    std::cout << "ASMatrixNakBspline: Solving failed.\n";
    return;
  }
  this->coefficients = alpha;
  this->grid = regularGrid;
  interpolantFlag = 1;
}

void ASMatrixNakBspline::createMatrix(size_t numPoints) { this->createMatrixMonteCarlo(numPoints); }

void ASMatrixNakBspline::createMatrixMonteCarlo(size_t numPoints) {
  if (interpolantFlag == 0) {
    std::cout << "ASMatrixNakBspline: cannot create Matrix without interpolant!\n";
    return;
  }
  size_t numDim = objectiveFunc.getNumberOfParameters();
  sgpp::optimization::InterpolantScalarFunctionGradient interpolantGradient(*(this->grid),
                                                                            this->coefficients);

  sgpp::base::DataVector randomVector(numDim, 1);
  this->C.resize(numDim, numDim);
  this->C.setZero();
  for (size_t i = 0; i < numPoints; ++i) {
    for (size_t d = 0; d < numDim; ++d) {
      randomVector[d] = RandomNumberGenerator::getInstance().getUniformRN(0.0, 1.0);
    }
    sgpp::base::DataVector gradient(numDim);
    interpolantGradient.eval(randomVector, gradient);
    Eigen::VectorXd e = DataVectorToEigen(gradient);
    this->C += e * e.transpose();
  }
  this->C /= static_cast<double>(numPoints);
}

}  // namespace optimization
}  // namespace sgpp

//#endif /* USE_EIGEN */
