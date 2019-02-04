// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <sgpp/datadriven/activeSubspaces/ASResponseSurfaceNakBspline.hpp>

namespace sgpp {
namespace datadriven {

void ASResponseSurfaceNakBspline::initialize() {
  activeDim = W1.cols();
  if (gridType == sgpp::base::GridType::NakBspline) {
    grid = std::make_shared<sgpp::base::NakBsplineGrid>(activeDim, degree);
    basis = std::make_shared<sgpp::base::SNakBsplineBase>(degree);
  } else if (gridType == sgpp::base::GridType::NakBsplineBoundary) {
    grid = std::make_shared<sgpp::base::NakBsplineBoundaryGrid>(activeDim, degree);
    basis = std::make_shared<sgpp::base::SNakBsplineBoundaryBase>(degree);
  } else if (gridType == sgpp::base::GridType::NakBsplineModified) {
    grid = std::make_shared<sgpp::base::NakBsplineModifiedGrid>(activeDim, degree);
    basis = std::make_shared<sgpp::base::SNakBsplineModifiedBase>(degree);
  } else if (gridType == sgpp::base::GridType::NakBsplineExtended) {
    grid = std::make_shared<sgpp::base::NakBsplineExtendedGrid>(activeDim, degree);
    basis = std::make_shared<sgpp::base::SNakBsplineExtendedBase>(degree);
  } else {
    throw sgpp::base::generation_exception("ASMatrixNakBspline: gridType not supported.");
  }
}

void ASResponseSurfaceNakBspline::createRegularReducedSurfaceFromData(
    sgpp::base::DataMatrix evaluationPoints, sgpp::base::DataVector functionValues, size_t level,
    double lambda) {
  Eigen::MatrixXd transformedPoints_Eigen = prepareDataForEigenRegression(evaluationPoints);

  grid->getGenerator().regular(level);
  coefficients = EigenRegression(grid, degree, transformedPoints_Eigen, functionValues, mse,
                                 errorPerBasis, lambda);

  interpolant =
      std::make_shared<sgpp::optimization::ASInterpolantScalarFunction>(*grid, coefficients);
  interpolantGradient = std::make_shared<sgpp::optimization::ASInterpolantScalarFunctionGradient>(
      *grid, coefficients);
}

void ASResponseSurfaceNakBspline::createRegularReducedSurfaceFromData_DataDriven(
    sgpp::base::DataMatrix evaluationPoints, sgpp::base::DataVector functionValues, size_t level,
    double lambda, double exponentBase) {
  Eigen::MatrixXd transformedPoints_Eigen = prepareDataForEigenRegression(evaluationPoints);
  sgpp::base::DataMatrix transformedPoints = EigenToDataMatrix(transformedPoints_Eigen);
  transformedPoints.transpose();

  //  const auto regularizationType = sgpp::datadriven::RegularizationType::Diagonal;
  const auto regularizationType = sgpp::datadriven::RegularizationType::Identity;
  auto regularizationConfig = sgpp::datadriven::RegularizationConfiguration();
  regularizationConfig.type_ = regularizationType;
  regularizationConfig.lambda_ = lambda;
  regularizationConfig.exponentBase_ = exponentBase;  // what is this for? does it do anything?

  auto gridConfig = sgpp::base::RegularGridConfiguration();
  gridConfig.dim_ = W1.cols();
  gridConfig.level_ = static_cast<int>(level);
  gridConfig.type_ = gridType;
  gridConfig.maxDegree_ = 3;

  auto adaptivityConfig = sgpp::base::AdpativityConfiguration();  // no adaptivity
  adaptivityConfig.noPoints_ = 0;
  adaptivityConfig.numRefinements_ = 0;

  auto solverConfig = sgpp::solver::SLESolverConfiguration();
  solverConfig.type_ = sgpp::solver::SLESolverType::CG;
  solverConfig.maxIterations_ = 1000;
  solverConfig.eps_ = 1e-8;
  solverConfig.threshold_ = 1e-5;
  auto learner = sgpp::datadriven::RegressionLearner(gridConfig, adaptivityConfig, solverConfig,
                                                     solverConfig, regularizationConfig);
  learner.train(transformedPoints, functionValues);
  grid = learner.getGridPtr();
  coefficients = learner.getWeights();

  //---------------
  //  sgpp::base::GridStorage& gridStroage = grid->getStorage();
  //  Eigen::MatrixXd regressionMatrix(evaluationPoints.getNcols(), gridStorage.getSize());
  //  for (size_t i = 0; i < evaluationPoints.getNcols(); i++) {
  //    sgpp::base::DataVector point(evaluationPoints.getNrows());
  //    evaluationPoints.getColumn(i, point);
  //
  //    Eigen::VectorXd y_i = W1.transpose() * DataVectorToEigen(point);
  //    // transformation in 1D case
  //    if (W1.cols() == 1) {
  //      y_i(0) = (y_i(0) - leftBound1D) / (rightBound1D - leftBound1D);
  //    }
  //
  //    for (size_t j = 0; j < gridStorage.getSize(); j++) {
  //      double basisEval = 1;
  //      for (size_t t = 0; t < gridStorage.getDimension(); t++) {
  //        double basisEval1D =
  //            basis->eval(gridStorage.getPointLevel(j, t), gridStorage.getPointIndex(j, t),
  //            y_i(t));
  //        if (basisEval1D == 0) {
  //          basisEval = 0;
  //          break;
  //        } else {
  //          basisEval *= basisEval1D;
  //        }
  //      }
  //      regressionMatrix(i, j) = basisEval;
  //    }
  //  }
  //  Eigen::VectorXd functionValues_Eigen = DataVectorToEigen(functionValues);
  //---------------

  interpolant =
      std::make_shared<sgpp::optimization::ASInterpolantScalarFunction>(*grid, coefficients);
  interpolantGradient = std::make_shared<sgpp::optimization::ASInterpolantScalarFunctionGradient>(
      *grid, coefficients);
}

void ASResponseSurfaceNakBspline::createAdaptiveReducedSurfaceFromData(
    size_t maxNumGridPoints, sgpp::base::DataMatrix evaluationPoints,
    sgpp::base::DataVector functionValues, size_t initialLevel, size_t refinementsNum,
    double lambda) {
  // number of points to be refined in each step
  grid->getGenerator().regular(initialLevel);
  std::shared_ptr<sgpp::optimization::ScalarFunction> transformedObjectiveFunc;
  // Special case one dimensional active subspace allows for reasonable response grid structure.
  // (Transform the one dimensional grid to an interval of according size). For active subspace
  // dimensions >1 we do not yet know how to transform the grid and recommend using regression on
  // the detection points
  Eigen::MatrixXd transformedPoints_Eigen = prepareDataForEigenRegression(evaluationPoints);
  coefficients = EigenRegression(grid, degree, transformedPoints_Eigen, functionValues, mse,
                                 errorPerBasis, lambda);
  while (grid->getSize() < maxNumGridPoints) {
    refineRegressionSurplusAdaptive(refinementsNum, transformedPoints_Eigen, functionValues,
                                    lambda);
    //    refineRegressionErrorAdaptive(refinementsNum, transformedPoints_Eigen, functionValues,
    //                                    lambda);
  }
  interpolant =
      std::make_shared<sgpp::optimization::ASInterpolantScalarFunction>(*grid, coefficients);
  interpolantGradient = std::make_shared<sgpp::optimization::ASInterpolantScalarFunctionGradient>(
      *grid, coefficients);
}

void ASResponseSurfaceNakBspline::createRegularReducedSurfaceWithPseudoInverse(
    size_t level, std::shared_ptr<sgpp::optimization::ScalarFunction> objectiveFunc) {
  grid->getGenerator().regular(level);
  // Special case one dimensional active subspace allows for reasonable response grid structure.
  // (Transform the one dimensional grid to an interval of according size). For active subspace
  // dimensions >1 we do not yet know how to transform the grid and recommend using regression on
  // the detection points
  if (W1.cols() == 1) {
    transformationfor1DActiveSubspace();
  }
  calculateInterpolationCoefficientsWithPseudoInverse(objectiveFunc);
  interpolant =
      std::make_shared<sgpp::optimization::ASInterpolantScalarFunction>(*grid, coefficients);
  interpolantGradient = std::make_shared<sgpp::optimization::ASInterpolantScalarFunctionGradient>(
      *grid, coefficients);

  // Check interpolation property g(p) = f(pinvW1' * p) for all grid points p
  //  Eigen::MatrixXd pinvW1 = W1.transpose().completeOrthogonalDecomposition().pseudoInverse();
  //  double diff = 0;
  //  for (size_t i = 0; i < grid->getSize(); i++) {
  //    auto gp = grid->getStorage().getPoint(i);
  //    Eigen::VectorXd p(grid->getDimension());
  //    for (size_t d = 0; d < grid->getDimension(); d++) {
  //      p[d] = gp.getStandardCoordinate(d);
  //    }
  //    sgpp::base::DataVector pinv = EigenToDataVector(pinvW1 * p);
  //    diff += fabs(objectiveFunc.eval(pinv) - interpolant->eval(EigenToDataVector(p)));
  //  }
  //  std::cout << "total diff: " << diff << std::endl;
}

void ASResponseSurfaceNakBspline::createAdaptiveReducedSurfaceWithPseudoInverse(
    size_t maxNumGridPoints, std::shared_ptr<sgpp::optimization::ScalarFunction> objectiveFunc,
    size_t initialLevel, size_t refinementsNum) {
  // number of points to be refined in each step
  grid->getGenerator().regular(initialLevel);
  std::shared_ptr<sgpp::optimization::ScalarFunction> transformedObjectiveFunc;
  // Special case one dimensional active subspace allows for reasonable response grid structure.
  // (Transform the one dimensional grid to an interval of according size). For active subspace
  // dimensions >1 we do not yet know how to transform the grid and recommend using regression on
  // the detection points
  if (W1.cols() == 1) {
    transformationfor1DActiveSubspace();
  }
  calculateInterpolationCoefficientsWithPseudoInverse(objectiveFunc);
  while (grid->getSize() < maxNumGridPoints) {
    this->refineInterpolationSurplusAdaptive(refinementsNum, objectiveFunc);
  }
  interpolant =
      std::make_shared<sgpp::optimization::ASInterpolantScalarFunction>(*grid, coefficients);
  interpolantGradient = std::make_shared<sgpp::optimization::ASInterpolantScalarFunctionGradient>(
      *grid, coefficients);
}

double ASResponseSurfaceNakBspline::eval(sgpp::base::DataVector v) {
  Eigen::VectorXd v_Eigen = sgpp::datadriven::DataVectorToEigen(v);
  Eigen::VectorXd trans_v_Eigen;
  // Special case one dimensional active subspace allows for reasonable response grid structure.
  // (Transform the one dimensional grid to an interval of according size).
  if (W1.cols() == 1) {
    trans_v_Eigen = W1.transpose() * v_Eigen;
    trans_v_Eigen(0) = (trans_v_Eigen(0) - leftBound1D) / (rightBound1D - leftBound1D);
  } else {
    trans_v_Eigen = W1.transpose() * v_Eigen;
  }
  sgpp::base::DataVector trans_v_DataVector = EigenToDataVector(trans_v_Eigen);
  return interpolant->eval(trans_v_DataVector);
}

double ASResponseSurfaceNakBspline::evalGradient(sgpp::base::DataVector v,
                                                 sgpp::base::DataVector& gradient) {
  Eigen::VectorXd v_Eigen = sgpp::datadriven::DataVectorToEigen(v);
  Eigen::VectorXd trans_v_Eigen;
  // Special case one dimensional active subspace allows for reasonable response grid structure.
  // (Transform the one dimensional grid to an interval of according size).
  if (W1.cols() == 1) {
    trans_v_Eigen = W1.transpose() * v_Eigen;
    trans_v_Eigen(0) = (trans_v_Eigen(0) - leftBound1D) / (rightBound1D - leftBound1D);
  } else {
    trans_v_Eigen = W1.transpose() * v_Eigen;
  }
  sgpp::base::DataVector trans_v_DataVector = EigenToDataVector(trans_v_Eigen);
  return interpolantGradient->eval(trans_v_DataVector, gradient);
}

double ASResponseSurfaceNakBspline::eval1D(double x) {
  sgpp::base::DataVector v(1, x);
  return interpolant->eval(v);
}

double ASResponseSurfaceNakBspline::getMCIntegral(size_t numMCPoints, size_t numHistogramMCPoints,
                                                  std::string pointStrategy) {
  if (W1.cols() != 1) {
    std::cerr << "ASResponseSurface::getMCIntegral currently supports only 1D active subspaces\n";
    return -1;
  }
  // uniform points spread over the active subspace (without boundary points)
  sgpp::base::DataVector points(numMCPoints);
  double delta = (rightBound1D - leftBound1D) / static_cast<double>(numMCPoints);
  for (size_t i = 0; i < numMCPoints; i++) {
    points[i] = leftBound1D + (static_cast<double>(2 * i + 1) / 2.0) * delta;
  }
  sgpp::base::DataVector weights =
      uniformIntervalHistogram(numHistogramMCPoints, points, delta, pointStrategy);
  sgpp::base::DataVector numVec(points.getSize(), static_cast<double>(numHistogramMCPoints));
  weights.componentwise_div(numVec);
  double integral = 0.0;
  sgpp::base::DataVector point(1, 0);
  for (size_t j = 0; j < numMCPoints; j++) {
    point[0] = (points[j] - leftBound1D) / (rightBound1D - leftBound1D);
    integral += weights[j] * interpolant->eval(point);
  }
  return integral;
}

double ASResponseSurfaceNakBspline::getHistogramBasedIntegral(size_t level,
                                                              size_t numHistogramMCPoints,
                                                              std::string pointStrategy) {
  if (W1.cols() != 1) {
    std::cerr << "ASResponseSurface::getHistogramBasedIntegral currently supports only 1D active "
                 "subspaces\n";
    return -1;
  }
  std::shared_ptr<sgpp::base::Grid> volGrid;
  sgpp::base::DataVector volCoefficients;
  size_t volDegree = 1;
  sgpp::base::GridType volGridType = sgpp::base::GridType::NakBsplineBoundary;
  histogramIntervalQuadrature(level, numHistogramMCPoints, volDegree, volGridType, volGrid,
                              volCoefficients, pointStrategy);

  double integral = getIntegralFromVolumeInterpolant(volGrid, volCoefficients, volDegree);
  return integral;
}

double ASResponseSurfaceNakBspline::getSplineBasedIntegral(size_t quadOrder) {
  if (W1.cols() != 1) {
    std::cerr << "ASResponseSurface::getSplineBasedIntegral supports only 1D active subspaces\n";
    return -1;
  }
  //  dd_set_global_constants();  // For cdd this must be called in the beginning.
  size_t dim = W1.rows();
  // matrix containing the orthogonally projected simplex corners (columnwise)
  Eigen::MatrixXd projectedCorners(dim + 1, factorial(dim));
  double simplexVolume = simplexDecomposition(projectedCorners);

  double integral = 0.0;
  // ToDo (rehmemk) local instances of grid, coefficients, gridType, degree for each thread?
  // #pragma omp parallel for reduction(+ : integral)
  for (unsigned int i = 0; i < projectedCorners.cols(); i++) {
    // the iterative M-spline definition needs sorted input knots
    sgpp::base::DataVector xi = EigenToDataVector(projectedCorners.col(i));
    std::sort(xi.begin(), xi.end());
    sgpp::datadriven::MSplineNakBsplineScalarProducts scalarProducts(gridType, degree, xi,
                                                                     quadOrder);

    double simplexIntegral = scalarProducts.calculateScalarProduct(grid, coefficients);

    integral += simplexIntegral;
  }
  integral *= simplexVolume;

  return integral;
}

double ASResponseSurfaceNakBspline::getApproximateSplineBasedIntegral(size_t approxLevel,
                                                                      size_t approxDegree) {
  if (W1.cols() != 1) {
    std::cerr << "ASResponseSurface::getSplineBasedIntegral  supports only 1D active subspaces\n";
    return -1;
  }
  //  dd_set_global_constants();  // For cdd this must be called in the beginning.
  size_t dim = W1.rows();
  // matrix containing the orthogonally projected simplex corners (columnwise)
  Eigen::MatrixXd projectedCorners(dim + 1, factorial(dim));
  sgpp::base::SGppStopwatch watch;
  watch.start();
  double simplexVolume = simplexDecomposition(projectedCorners);
  //  std::cout << "simplexDecomposition " << watch.stop() << "s\n";

  auto volGrid = std::make_shared<sgpp::base::NakBsplineBoundaryGrid>(1, approxDegree);
  volGrid->getGenerator().regular(approxLevel);
  sgpp::base::GridStorage& volGridStorage = volGrid->getStorage();
  sgpp::base::DataVector interpolPoints(volGridStorage.getSize());
  double width = rightBound1D - leftBound1D;
  for (size_t i = 0; i < volGridStorage.getSize(); i++) {
    interpolPoints[i] = leftBound1D + width * volGridStorage.getPointCoordinate(i, 0);
  }
  watch.start();
  sgpp::base::DataVector volumes =
      caclculateVolumeSimplexWise(interpolPoints, simplexVolume, projectedCorners);
  //  std::cout << "caclculateVolumeSimplexWise " << watch.stop() << "s\n";

  sgpp::optimization::HierarchisationSLE hierSLE(*volGrid);
  sgpp::optimization::sle_solver::Armadillo sleSolver;
  sgpp::base::DataVector volCoefficients;
  if (!sleSolver.solve(hierSLE, volumes, volCoefficients)) {
    std::cerr << "ASMatrixNakBspline: Solving failed.\n";
  }
  double integral = getIntegralFromVolumeInterpolant(volGrid, volCoefficients, approxDegree);

  return integral;
}

double ASResponseSurfaceNakBspline::getIntegralFromVolumeInterpolant(
    std::shared_ptr<sgpp::base::Grid> volGrid, sgpp::base::DataVector volCoefficients,
    size_t volDegree) {
  if (W1.cols() != 1) {
    std::cerr << "ASResponseSurface::getIntegralFromVolumeInterpolant currently supports only 1D "
                 "active subspaces\n";
    return -1;
  }
  sgpp::base::GridType volGridType = volGrid->getType();
  size_t quadOrder = degree + volDegree;

  sgpp::datadriven::NakBsplineScalarProducts scalarProducts(gridType, volGridType, degree,
                                                            volDegree, quadOrder);
  double integral =
      scalarProducts.calculateScalarProduct(grid, coefficients, volGrid, volCoefficients);
  return integral * (rightBound1D - leftBound1D);
}

// ----------------- auxiliary routines -----------
int ASResponseSurfaceNakBspline::factorial(size_t n) {
  return (n == 1 || n == 0) ? 1 : factorial(n - 1) * n;
}

// dd_MatrixPtr ASResponseSurfaceNakBspline::createHPolytope(std::vector<int> permutations) {
//  size_t dim = permutations.size();
//  dd_MatrixPtr hMatrix = dd_CreateMatrix(dim + 3, dim + 1);
//
//  // x_p[0] >= 0
//  hMatrix->matrix[0][permutations[0] + 1][0] = 1.0;
//  // x_p[end] <= 1
//  hMatrix->matrix[1][0][0] = 1.0;
//  hMatrix->matrix[1][permutations.back() + 1][0] = -1.0;
//  // x_p[i] <= x_p[i+1]
//  for (size_t d = 0; d < dim - 1; d++) {
//    hMatrix->matrix[d + 2][permutations[d] + 1][0] = -1.0;
//    hMatrix->matrix[d + 2][permutations[d + 1] + 1][0] = 1.0;
//  }
//  hMatrix->numbtype = dd_Real;
//  hMatrix->representation = dd_Inequality;
//  return hMatrix;
//}

// ToDo(rehmemk) isn't it somehow possible to directly generate the vertices of the simplices
// instead of creating a H(yperplane)representation and transforming that to a
// V(ertex)-representation as done here?
double ASResponseSurfaceNakBspline::simplexDecomposition(Eigen::MatrixXd& projectedCorners) {
  size_t dim = W1.rows();
  double simplexVolume = 0.0;
  projectedCorners.resize(dim + 1, factorial(dim));
  std::vector<int> permutations(dim);
  std::iota(permutations.begin(), permutations.end(), 0);
  size_t i = 0;
  do {
    //    dd_MatrixPtr hMatrix = createHPolytope(permutations);
    //    dd_ErrorType err;
    //    dd_PolyhedraPtr poly = dd_DDMatrix2Poly(hMatrix, &err);
    //    dd_FreeMatrix(hMatrix);
    //
    //    dd_MatrixPtr vRep = dd_CopyGenerators(poly);
    //    dd_FreePolyhedra(poly);
    //    dd_Amatrix vRepMatrix = vRep->matrix;
    //
    //    Eigen::MatrixXd V(vRep->colsize - 1, vRep->rowsize);  // simplex points (columnwise)
    //    for (unsigned int i = 0; i < V.rows(); i++) {
    //      for (unsigned int j = 0; j < V.cols(); j++) {
    //        V(i, j) = vRepMatrix[j][i + 1][0];
    //      }
    //    }

    // simplex points (column wise)
    Eigen::MatrixXd V = Eigen::MatrixXd::Zero(dim, dim + 1);
    for (size_t k = 0; k < dim; k++) {
      for (size_t l = k; l < dim; l++) {
        V(permutations[l], k) = 1.0;
      }
    }

    //    std::cout << "\n";
    //    for (auto& p : permutations) {
    //      std::cout << p << " ";
    //    }
    //    std::cout << "\n";
    //    std::cout << V << "\n";

    projectedCorners.col(i) = (W1.transpose() * V).transpose();

    // calculate only once as it is equal for all simplices
    if (i == 0) {
      Eigen::MatrixXd D(dim, dim);
      for (unsigned int k = 1; k <= dim; k++) {
        D.col(k - 1) = V.col(k) - V.col(0);
      }
      simplexVolume = abs(D.determinant() / static_cast<double>(factorial(dim)));
    }
    i++;
  } while (std::next_permutation(permutations.begin(), permutations.end()));
  return simplexVolume;
}  // namespace datadriven

sgpp::base::DataVector ASResponseSurfaceNakBspline::caclculateVolumeSimplexWise(
    sgpp::base::DataVector points, double simplexVolume, Eigen::MatrixXd projectedCorners) {
  sgpp::base::DataVector volumes(points.getSize(), 0.0);
  MSplineBasis mSplineBasis;
#pragma omp parallel for private(mSplineBasis)
  for (size_t p = 0; p < points.getSize(); p++) {
    double res = 0.0;
    for (unsigned int i = 0; i < projectedCorners.cols(); i++) {
      // the iterative M-spline definition needs sorted input knots
      sgpp::base::DataVector xi = EigenToDataVector(projectedCorners.col(i));
      std::sort(xi.begin(), xi.end());
      mSplineBasis.setXi(xi);
      res += mSplineBasis.eval(xi.getSize() - 1, 0, points[p]);
    }
    volumes[p] = res;
  }
  volumes.mult(simplexVolume);
  return volumes;
}

sgpp::base::DataVector ASResponseSurfaceNakBspline::uniformIntervalHistogram(
    size_t numHistogramMCPoints, sgpp::base::DataVector points, double delta,
    std::string pointStrategy) {
  // the weight for each point p_i is the proportion of all random points in point p_i's bucket
  //  when transformed  to the active subspace (by multiplication with W1)
  // p_i's bucket is [(p_{i-1}+p_i)/2, (p_i,p_{i+1})/2]
  sgpp::base::DataVector weights(points.getSize());
  Eigen::VectorXd randomUnitPoint(W1.rows());
  Eigen::VectorXd randomPoint(1);
  double hdelta = delta / 2.0;

  if (pointStrategy == "MC") {
    // random Monte Carlo points in the unit hypercube
    Eigen::MatrixXd randomUnitPoints = (Eigen::MatrixXd::Random(W1.rows(), numHistogramMCPoints) +
                                        Eigen::MatrixXd::Ones(W1.rows(), numHistogramMCPoints)) /
                                       2;
    Eigen::MatrixXd randomPoints = W1.transpose() * randomUnitPoints;

    for (size_t i = 0; i < numHistogramMCPoints; i++) {
      for (size_t j = 0; j < points.getSize(); j++) {
        if ((randomPoints(i) >= (points[j] - hdelta)) && (randomPoints(i) < (points[j] + hdelta))) {
          weights[j]++;
          break;
        }
      }
    }
  } else if (pointStrategy == "Halton") {
    // Quasi Monte Carlo with Halton Sequence in the unit hypercube
    // ( Wikipedia: "Halton sequence works best up to ~ 6 dimensions, then Sobol sequence
    // performs better")
    int* haltonsequence = new int[W1.rows()];
    for (int i = 0; i < W1.rows(); i++) {
      haltonsequence[i] = sgpp::datadriven::prime(i + 1);
    }
    double* haltonPoint;
    for (size_t j = 0; j < numHistogramMCPoints; j++) {
      haltonPoint = sgpp::datadriven::halton_base(static_cast<int>(j), static_cast<int>(W1.rows()),
                                                  haltonsequence);
      for (int i = 0; i < W1.rows(); i++) {
        randomUnitPoint(i) = haltonPoint[i];
      }
      randomPoint = W1.transpose() * randomUnitPoint;
      for (size_t j = 0; j < points.getSize(); j++) {
        if ((randomPoint(0) >= (points[j] - hdelta)) && (randomPoint(0) < (points[j] + hdelta))) {
          weights[j]++;
          break;
        }
      }
    }
    delete[] haltonsequence;
  } else if (pointStrategy == "Sobol") {
    // Quasi Monte Carlo with Sobol sequence in the unit hypercube
    double* sobolPoint = new double[W1.cols()];
    // ToDo (rehmemk) Is this seed really good? For example simple2D seed 12345 was much better
    // !?
    long long int seed = sgpp::datadriven::tau_sobol(static_cast<int>(W1.cols()));
    for (size_t j = 0; j < numHistogramMCPoints; j++) {
      sgpp::datadriven::i8_sobol(static_cast<int>(W1.rows()), &seed, sobolPoint);
      for (int i = 0; i < W1.rows(); i++) {
        randomUnitPoint(i) = sobolPoint[i];
      }
      randomPoint = W1.transpose() * randomUnitPoint;
      for (size_t j = 0; j < points.getSize(); j++) {
        if ((randomPoint(0) >= (points[j] - hdelta)) && (randomPoint(0) < (points[j] + hdelta))) {
          weights[j]++;
          break;
        }
      }
    }
    delete[] sobolPoint;
  } else {
    std::cerr << "ASResponseSurfaceNakBspline: pointStrategy not supported!\n";
  }

  return weights;
}

void ASResponseSurfaceNakBspline::histogramIntervalQuadrature(
    size_t level, size_t numHistogramMCPoints, size_t quadDegree, sgpp::base::GridType quadGridType,
    std::shared_ptr<sgpp::base::Grid>& quadGrid, sgpp::base::DataVector& quadCoefficients,
    std::string pointStrategy) {
  if (W1.cols() != 1) {
    std::cerr << "ASResponseSurface::continuousIntervalQuadrature currently supports only 1D "
                 "active subspaces\n";
    return;
  }
  size_t dim = 1;

  if (quadGridType == sgpp::base::GridType::NakBsplineBoundary) {
    // this is only for 1D. So we don't care about the two extra boundary points
    quadGrid.reset(new sgpp::base::NakBsplineBoundaryGrid(dim, quadDegree));
  } else {
    throw sgpp::base::generation_exception("ASMatrixNakBspline: gridType not supported.");
  }
  quadGrid->getGenerator().regular(level);
  sgpp::base::GridStorage& quadGridStorage = quadGrid->getStorage();
  sgpp::base::DataVector points(quadGridStorage.getSize());
  for (size_t i = 0; i < quadGridStorage.getSize(); i++) {
    points[i] =
        leftBound1D + (rightBound1D - leftBound1D) * quadGridStorage.getPointCoordinate(i, 0);
  }
  double delta = abs((rightBound1D - leftBound1D) / std::pow(2, level));
  sgpp::base::DataVector weights =
      uniformIntervalHistogram(numHistogramMCPoints, points, delta, pointStrategy);
  sgpp::base::DataVector numVec(points.getSize(),
                                delta * static_cast<double>(numHistogramMCPoints));
  weights.componentwise_div(numVec);

  sgpp::optimization::HierarchisationSLE hierSLE(*quadGrid);
  sgpp::optimization::sle_solver::Armadillo sleSolver;
  if (!sleSolver.solve(hierSLE, weights, quadCoefficients)) {
    std::cout << "ASMatrixNakBspline: Solving failed.\n";
    return;
  }
}

void ASResponseSurfaceNakBspline::refineInterpolationSurplusAdaptive(
    size_t refinementsNum, std::shared_ptr<sgpp::optimization::ScalarFunction> objectiveFunc) {
  sgpp::base::SurplusRefinementFunctor functor(coefficients, refinementsNum);
  grid->getGenerator().refine(functor);
  calculateInterpolationCoefficientsWithPseudoInverse(objectiveFunc);
}

void ASResponseSurfaceNakBspline::calculateInterpolationCoefficientsWithPseudoInverse(
    std::shared_ptr<sgpp::optimization::ScalarFunction> objectiveFunc) {
  sgpp::base::GridStorage& gridStorage = grid->getStorage();
  sgpp::base::DataVector functionValues(grid->getSize());
  sgpp::optimization::HierarchisationSLE hierSLE(*grid);

  for (size_t i = 0; i < gridStorage.getSize(); i++) {
    Eigen::VectorXd p(static_cast<int>(gridStorage.getDimension()));
    p = DataVectorToEigen(gridStorage.getPointCoordinates(i));

    // transformation in 1D case
    if (W1.cols() == 1) {
      p(0) = p(0) * (rightBound1D - leftBound1D) + leftBound1D;
    }

    // computes approximate solution "pseudo inverse style".
    sgpp::base::DataVector pinv = EigenToDataVector(
        W1.transpose().jacobiSvd(Eigen::ComputeThinU | Eigen::ComputeThinV).solve(p));

    functionValues[i] = objectiveFunc->eval(pinv);
  }

  sgpp::optimization::sle_solver::Auto sleSolver;
  if (!sleSolver.solve(hierSLE, functionValues, coefficients)) {
    std::cout << "Solving failed, exiting.\n";
    return;
  }
}

void ASResponseSurfaceNakBspline::refineRegressionSurplusAdaptive(
    size_t refinementsNum, Eigen::MatrixXd transformedPoints_Eigen,
    sgpp::base::DataVector functionValues, double lambda) {
  sgpp::base::SurplusRefinementFunctor functor(coefficients, refinementsNum);
  grid->getGenerator().refine(functor);
  coefficients = EigenRegression(grid, degree, transformedPoints_Eigen, functionValues, mse,
                                 errorPerBasis, lambda);
}

void ASResponseSurfaceNakBspline::refineRegressionErrorAdaptive(
    size_t refinementsNum, Eigen::MatrixXd transformedPoints_Eigen,
    sgpp::base::DataVector functionValues, double lambda) {
  sgpp::base::SurplusRefinementFunctor functor(errorPerBasis, refinementsNum);
  grid->getGenerator().refine(functor);
  coefficients = EigenRegression(grid, degree, transformedPoints_Eigen, functionValues, mse,
                                 errorPerBasis, lambda);
}

Eigen::MatrixXd ASResponseSurfaceNakBspline::hypercubeVertices(size_t dimension) {
  size_t twoDim = static_cast<size_t>(std::pow(2, dimension));
  Eigen::MatrixXd corners = Eigen::MatrixXd::Zero(
      static_cast<unsigned int>(dimension), static_cast<unsigned int>(std::pow(2, dimension)));
  if (dimension <= 0) {
    throw sgpp::base::algorithm_exception("ASResponseSurfaceNakBspline: dimension must be > 0!");
  } else if (dimension == 1) {
    corners(0, 0) = 1;
    corners(0, 1) = 0;
    return corners;
  }
  size_t jump = static_cast<size_t>(std::pow(2, dimension - 1));
  for (size_t i = 0; i < dimension; i++) {
    size_t j = 0;
    while (j < twoDim) {
      for (size_t n = 0; n < jump; n++) {
        if (j + n >= twoDim) {
          break;
        }
        corners(i, j + n) = 1;
      }
      j = j + 2 * jump;
    }
    jump = jump / 2;
  }
  return corners;
}

void ASResponseSurfaceNakBspline::transformationfor1DActiveSubspace() {
  // find maximum value "rightBound" of W1^T * x for all corners x of [0,1]^d,
  // then set up an interpolant on [leftBound,rightBound]
  Eigen::MatrixXd corners = hypercubeVertices(W1.rows());
  rightBound1D = std::numeric_limits<double>::lowest();
  leftBound1D = std::numeric_limits<double>::max();
  for (unsigned int j = 0; j < corners.cols(); j++) {
    double tempBound = (W1.transpose() * corners.col(j))(0, 0);
    rightBound1D = tempBound > rightBound1D ? tempBound : rightBound1D;
    leftBound1D = tempBound < leftBound1D ? tempBound : leftBound1D;
  }
}

Eigen::MatrixXd ASResponseSurfaceNakBspline::prepareDataForEigenRegression(
    sgpp::base::DataMatrix evaluationPoints) {
  Eigen::MatrixXd evaluationPoints_Eigen = DataMatrixToEigen(evaluationPoints);
  Eigen::MatrixXd transformedPoints_Eigen = W1.transpose() * evaluationPoints_Eigen;
  if (W1.cols() == 1) {
    transformationfor1DActiveSubspace();
    Eigen::MatrixXd leftMatrix =
        leftBound1D *
        Eigen::MatrixXd::Ones(transformedPoints_Eigen.rows(), transformedPoints_Eigen.cols());
    transformedPoints_Eigen = (transformedPoints_Eigen - leftMatrix) / (rightBound1D - leftBound1D);
  }
  return transformedPoints_Eigen;
}

}  // namespace datadriven
}  // namespace sgpp
