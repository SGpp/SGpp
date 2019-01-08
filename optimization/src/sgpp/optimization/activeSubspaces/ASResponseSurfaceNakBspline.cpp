// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <sgpp/optimization/activeSubspaces/ASResponseSurfaceNakBspline.hpp>

namespace sgpp {
namespace optimization {

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

void ASResponseSurfaceNakBspline::createRegularReducedSurfaceFromDetectionPoints(
    sgpp::base::DataMatrix evaluationPoints, sgpp::base::DataVector functionValues, size_t level) {
  sgpp::base::GridStorage& gridStorage = grid->getStorage();
  grid->getGenerator().regular(level);

  // interpolationMatrix(i,j) = b_j (y_i) = b_j (W1T * evaluationPoint_i)
  Eigen::MatrixXd interpolationMatrix(evaluationPoints.getNrows(), gridStorage.getSize());
  for (size_t i = 0; i < evaluationPoints.getNrows(); i++) {
    sgpp::base::DataVector point(evaluationPoints.getNcols());
    evaluationPoints.getRow(i, point);
    Eigen::VectorXd y_i = W1.transpose() * DataVectorToEigen(point);
    for (size_t j = 0; j < gridStorage.getSize(); j++) {
      double basisEval = 1;
      for (size_t t = 0; t < gridStorage.getDimension(); t++) {
        double basisEval1D =
            basis->eval(gridStorage.getPointLevel(j, t), gridStorage.getPointIndex(j, t), y_i(t));
        if (basisEval1D == 0) {
          basisEval = 0;
          break;
        } else {
          basisEval *= basisEval1D;
        }
      }
      interpolationMatrix(i, j) = basisEval;
    }
  }

  Eigen::VectorXd functionValues_Eigen = DataVectorToEigen(functionValues);
  // Least Squares Fit
  Eigen::VectorXd alpha_Eigen =
      interpolationMatrix.bdcSvd(Eigen::ComputeThinU | Eigen::ComputeThinV)
          .solve(functionValues_Eigen);

  coefficients = EigenToDataVector(alpha_Eigen);
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
    this->refineSurplusAdaptive(refinementsNum, objectiveFunc);
  }
  interpolant =
      std::make_shared<sgpp::optimization::ASInterpolantScalarFunction>(*grid, coefficients);
  interpolantGradient = std::make_shared<sgpp::optimization::ASInterpolantScalarFunctionGradient>(
      *grid, coefficients);
}

double ASResponseSurfaceNakBspline::eval(sgpp::base::DataVector v) {
  Eigen::VectorXd v_Eigen = sgpp::optimization::DataVectorToEigen(v);
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
  Eigen::VectorXd v_Eigen = sgpp::optimization::DataVectorToEigen(v);
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
  /*
  double integral = 0.0;
    sgpp::base::DataVector points(numMCPoints + 1);
    double delta = (rightBound1D - leftBound1D) / static_cast<double>(numMCPoints);
    for (size_t i = 0; i < numMCPoints + 1; i++) {
      points[i] = leftBound1D + static_cast<double>(i) * delta;
    }
    sgpp::base::DataVector weights(numMCPoints, 0.0);

    Eigen::MatrixXd randomUnitPoints = Eigen::MatrixXd::Random(W1.rows(), numHistogramMCPoints);
    randomUnitPoints += Eigen::MatrixXd::Ones(W1.rows(), numHistogramMCPoints);
    randomUnitPoints /= 2;
    Eigen::MatrixXd randomPoints = W1.transpose() * randomUnitPoints;

    for (size_t i = 0; i < numHistogramMCPoints; i++) {
      for (size_t j = 0; j < numMCPoints; j++) {
        if (randomPoints(i) >= points[j] && randomPoints(i) < points[j + 1]) {
          weights[j]++;
          break;
        }
      }
    }
    sgpp::base::DataVector numVec(numMCPoints, static_cast<double>(numHistogramMCPoints));
    weights.componentwise_div(numVec);

    sgpp::base::DataVector point(1, 0);
    for (size_t j = 0; j < numMCPoints; j++) {
      point[0] = (points[j] + points[j + 1]) / 2;
      point[0] = (point[0] - leftBound1D) / (rightBound1D - leftBound1D);
      //    std::cout << point[0] << " " << interpolant->eval(point) << " " << weights[j] << "\n";
      integral += weights[j] * interpolant->eval(point);
    }*/

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
    //    std::cout << point[0] << " " << interpolant->eval(point) << " " << weights[j] << "\n";
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
  sgpp::base::GridStorage& gridStorage = grid->getStorage();
  std::shared_ptr<sgpp::base::Grid> volGrid;
  sgpp::base::DataVector volCoefficients;
  size_t volDegree = 1;
  sgpp::base::GridType volGridType = sgpp::base::GridType::NakBsplineBoundary;
  histogramIntervalQuadrature(level, numHistogramMCPoints, volDegree, volGridType, volGrid,
                              volCoefficients, pointStrategy);

  // ---------DEBUG-------
  //  sgpp::optimization::InterpolantScalarFunction quadInt(*quadGrid, quadCoefficients);
  //  sgpp::optimization::InterpolantScalarFunction Int(*grid, coefficients);
  //  sgpp::base::DataVector pointi(1, 0.0);
  //  size_t npi = 100;
  //  sgpp::base::DataVector quadRes(npi + 1);
  //  sgpp::base::DataVector Res(npi + 1);
  //  for (size_t i = 0; i < npi + 1; i++) {
  //    pointi[0] = static_cast<double>(i) / static_cast<double>(npi);
  //    quadRes[i] = quadInt.eval(pointi);
  //    Res[i] = Int.eval(pointi);
  //  }
  //  std::cout << "qRes:\n" << quadRes.toString() << "\n";
  //  std::cout << "Res:\n" << Res.toString() << "\n";
  // ------------

  double integral = getIntegralFromVolumeInterpolant(volGrid, volCoefficients, volDegree);
  return integral;
}

// ToDo (rehmemk) nterpolating the ;-Spline with B-splines is of course nonsense. Rewrite the scalar
// product routine s.t. it allows directly using the M-spline
double ASResponseSurfaceNakBspline::getSplineBasedIntegral(size_t level, size_t volDegree) {
  if (W1.cols() != 1) {
    std::cerr << "ASResponseSurface::getSplineBasedIntegral currently supports only 1D active "
                 "subspaces\n";
    return -1;
  }
  dd_set_global_constants();  // For cdd this must be called in the beginning.

  // matrix containing the orthogonally projected simplex corners (columnwise)
  size_t dim = W1.rows();
  Eigen::MatrixXd projectedCorners(dim + 1, factorial(dim));
  double simplexVolume = simplexWiseVolume(projectedCorners);
  //  std::cout << "projected corners:\n" << projectedCorners << "\n";
  //  std::cout << "simplex volumes:\n" << simplexVolumes << "\n---\n";

  //  size_t nPoints = 100;
  //  Eigen::VectorXd B(nPoints);
  //  Eigen::VectorXd X(nPoints);
  //  for (size_t i = 0; i < nPoints; i++) {
  //    double x = static_cast<double>(i) / static_cast<double>(nPoints) * 1.31306;
  //    X(i) = x;
  //    B(i) = evalSimplexWiseVolume(x, projectedCorners, simplexVolumes);
  //  }
  //  std::cout << "X:\n" << X.transpose() << "\nB:\n" << B.transpose() << "\n";

  auto volGrid = std::make_shared<sgpp::base::NakBsplineBoundaryGrid>(1, volDegree);
  volGrid->getGenerator().regular(level);
  sgpp::base::GridStorage& volGridStorage = volGrid->getStorage();
  sgpp::base::DataVector volumes(volGridStorage.getSize());
  for (size_t i = 0; i < volGridStorage.getSize(); i++) {
    //  size_t i = 7;
    double point =
        leftBound1D + (rightBound1D - leftBound1D) * volGridStorage.getPointCoordinate(i, 0);
    //  std::cout << "point: " << point << "\n";
    volumes[i] = evalSimplexWiseVolume(point, simplexVolume, projectedCorners);
    //  std::cout << "vol[i]: " << volumes[i] << "\n";
  }

  sgpp::optimization::HierarchisationSLE hierSLE(*volGrid);
  sgpp::optimization::sle_solver::Armadillo sleSolver;
  sgpp::base::DataVector volCoefficients;
  if (!sleSolver.solve(hierSLE, volumes, volCoefficients)) {
    std::cerr << "ASMatrixNakBspline: Solving failed.\n";
  }
  double integral = getIntegralFromVolumeInterpolant(volGrid, volCoefficients, volDegree);
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
  double integral = 0.0;
  sgpp::base::GridStorage& gridStorage = grid->getStorage();

  sgpp::base::GridStorage& volGridStorage = volGrid->getStorage();
  sgpp::base::GridType volGridType = volGrid->getType();
  size_t quadOrder = static_cast<size_t>(std::ceil(static_cast<double>(degree) + 1.0 / 2.0)) * 2;

  sgpp::optimization::NakBsplineScalarProducts scalarProducts(gridType, volGridType, degree,
                                                              volDegree, quadOrder);
  for (size_t k = 0; k < coefficients.getSize(); k++) {
    for (size_t l = 0; l < volCoefficients.getSize(); l++) {
      integral += coefficients[k] * volCoefficients[l] *
                  scalarProducts.univariateScalarProduct(gridStorage.getPointLevel(k, 0),
                                                         gridStorage.getPointIndex(k, 0), false,
                                                         volGridStorage.getPointLevel(l, 0),
                                                         volGridStorage.getPointIndex(l, 0), false);
    }
  }
  return integral * (rightBound1D - leftBound1D);
}

// ----------------- auxiliary routines -----------
int ASResponseSurfaceNakBspline::factorial(size_t n) {
  return (n == 1 || n == 0) ? 1 : factorial(n - 1) * n;
}

double ASResponseSurfaceNakBspline::xpowplus(double x, size_t n) {
  if (x >= 0)
    return std::pow(x, static_cast<double>(n));
  else
    return 0.0;
}

double ASResponseSurfaceNakBspline::w(size_t v, Eigen::VectorXd xi) {
  double res = 1.0;
  for (unsigned int i = 0; i < xi.size(); i++) {
    if (i != v) {
      res *= (xi(v) - xi(i));
    }
  }
  return res;
}

double ASResponseSurfaceNakBspline::Mspline(double x, Eigen::VectorXd xi) {
  double res = 0.0;
  size_t n = xi.size() - 1;
  for (size_t v = 0; v < n + 1; v++) {
    //    std::cout << "v " << v << " x+ " << xpowplus(xi[v] - x, n - 1) << " w " << w(v, xi) <<
    //    "\n";
    res += (static_cast<double>(n) * xpowplus(xi[v] - x, n - 1)) / w(v, xi);
  }
  return res;
}

dd_MatrixPtr ASResponseSurfaceNakBspline::createHPolytope(std::vector<int> permutations) {
  size_t dim = permutations.size();
  dd_MatrixPtr hMatrix = dd_CreateMatrix(dim + 3, dim + 1);

  // x_p[0] >= 0
  hMatrix->matrix[0][permutations[0] + 1][0] = 1.0;
  // x_p[end] <= 1
  hMatrix->matrix[1][0][0] = 1.0;
  hMatrix->matrix[1][permutations.back() + 1][0] = -1.0;
  // x_p[i] <= x_p[i+1]
  for (size_t d = 0; d < dim - 1; d++) {
    hMatrix->matrix[d + 2][permutations[d] + 1][0] = -1.0;
    hMatrix->matrix[d + 2][permutations[d + 1] + 1][0] = 1.0;
  }
  hMatrix->numbtype = dd_Real;
  hMatrix->representation = dd_Inequality;
  return hMatrix;
}

double ASResponseSurfaceNakBspline::simplexWiseVolume(Eigen::MatrixXd& projectedCorners) {
  size_t dim = W1.rows();
  double simplexVolume = 0.0;
  projectedCorners.resize(dim + 1, factorial(dim));
  std::vector<int> permutations(dim);
  std::iota(permutations.begin(), permutations.end(), 0);
  size_t i = 0;
  do {
    //    for (size_t d = 0; d < dim; d++) {
    //      std::cout << permutations[d] << " ";
    //    }
    //    std::cout << ":\n";

    dd_MatrixPtr hMatrix = createHPolytope(permutations);
    dd_ErrorType err;
    dd_PolyhedraPtr poly = dd_DDMatrix2Poly(hMatrix, &err);
    dd_FreeMatrix(hMatrix);

    dd_MatrixPtr vRep = dd_CopyGenerators(poly);
    dd_FreePolyhedra(poly);
    dd_Amatrix vRepMatrix = vRep->matrix;

    Eigen::MatrixXd V(vRep->colsize - 1, vRep->rowsize);  // simplex points (columnwise)
    for (unsigned int i = 0; i < V.rows(); i++) {
      for (unsigned int j = 0; j < V.cols(); j++) {
        V(i, j) = vRepMatrix[j][i + 1][0];
      }
    }

    projectedCorners.col(i) = (W1.transpose() * V).transpose();

    // calculate only one as it is equal for all simplices
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
}

double ASResponseSurfaceNakBspline::evalSimplexWiseVolume(double x, double simplexVolume,
                                                          Eigen::MatrixXd projectedCorners) {
  double res = 0.0;
  for (unsigned int i = 0; i < projectedCorners.cols(); i++) {
    res += Mspline(x, projectedCorners.col(i));
  }
  return simplexVolume * res;
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
      haltonsequence[i] = sgpp::optimization::prime(i + 1);
    }
    double* haltonPoint;
    for (size_t j = 0; j < numHistogramMCPoints; j++) {
      haltonPoint = sgpp::optimization::halton_base(static_cast<int>(j),
                                                    static_cast<int>(W1.rows()), haltonsequence);
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
    long long int seed = sgpp::optimization::tau_sobol(static_cast<int>(W1.cols()));
    for (size_t j = 0; j < numHistogramMCPoints; j++) {
      sgpp::optimization::i8_sobol(static_cast<int>(W1.rows()), &seed, sobolPoint);
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

void ASResponseSurfaceNakBspline::refineSurplusAdaptive(
    size_t refinementsNum, std::shared_ptr<sgpp::optimization::ScalarFunction> objectiveFunc) {
  sgpp::base::SurplusRefinementFunctor functor(coefficients, refinementsNum);
  grid->getGenerator().refine(functor);
  calculateInterpolationCoefficientsWithPseudoInverse(objectiveFunc);
}

void ASResponseSurfaceNakBspline::calculateInterpolationCoefficientsWithPseudoInverse(
    std::shared_ptr<sgpp::optimization::ScalarFunction> objectiveFunc) {
  Eigen::MatrixXd pinvW1 = W1.transpose().completeOrthogonalDecomposition().pseudoInverse();
  sgpp::base::GridStorage& gridStorage = grid->getStorage();
  sgpp::base::DataVector functionValues(grid->getSize());
  sgpp::optimization::HierarchisationSLE hierSLE(*grid);

  for (size_t i = 0; i < gridStorage.getSize(); i++) {
    Eigen::VectorXd p(static_cast<int>(gridStorage.getDimension()));
    for (size_t d = 0; d < gridStorage.getDimension(); d++) {
      p[d] = gridStorage.getPointCoordinate(i, d);

      // transformation in 1D case
      if (W1.cols() == 1) {
        p(0) = p(0) * (rightBound1D - leftBound1D) + leftBound1D;
      }
      sgpp::base::DataVector pinv =
          EigenToDataVector(pinvW1 * p);  // introduce a wrapper for eigen functions so
                                          // we don't have to transform here every time?
      functionValues[i] = objectiveFunc->eval(pinv);
    }
  }
  sgpp::optimization::sle_solver::Auto sleSolver;
  if (!sleSolver.solve(hierSLE, functionValues, coefficients)) {
    std::cout << "Solving failed, exiting.\n";
    return;
  }
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
  std::cout << "left Bound: " << leftBound1D << " right Bound: " << rightBound1D << "\n";
}

}  // namespace optimization
}  // namespace sgpp
