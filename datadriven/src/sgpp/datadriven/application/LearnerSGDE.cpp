// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <sgpp/datadriven/application/LearnerSGDE.hpp>

#include <sgpp/base/exception/application_exception.hpp>
#include <sgpp/base/exception/data_exception.hpp>
#include <sgpp/base/grid/GridStorage.hpp>
#include <sgpp/base/grid/generation/GridGenerator.hpp>
#include <sgpp/base/grid/generation/functors/SurplusRefinementFunctor.hpp>
#include <sgpp/base/grid/storage/hashmap/HashGridStorage.hpp>
#include <sgpp/base/operation/BaseOpFactory.hpp>
#include <sgpp/base/operation/hash/OperationEval.hpp>
#include <sgpp/base/operation/hash/OperationFirstMoment.hpp>
#include <sgpp/base/operation/hash/OperationMultipleEval.hpp>
#include <sgpp/base/tools/json/json_exception.hpp>
#include <sgpp/datadriven/algorithm/DensitySystemMatrix.hpp>
#include <sgpp/datadriven/algorithm/RefinementMonitor.hpp>
#include <sgpp/datadriven/algorithm/RefinementMonitorConvergence.hpp>
#include <sgpp/datadriven/algorithm/RefinementMonitorPeriodic.hpp>
#include <sgpp/datadriven/functors/MultiGridRefinementFunctor.hpp>
#include <sgpp/datadriven/functors/classification/DataBasedRefinementFunctor.hpp>
#include <sgpp/datadriven/functors/classification/ZeroCrossingRefinementFunctor.hpp>
#include <sgpp/datadriven/operation/hash/simple/OperationCovariance.hpp>
#include <sgpp/pde/operation/PdeOpFactory.hpp>
#include <sgpp/solver/TypesSolver.hpp>
#include <sgpp/solver/sle/ConjugateGradients.hpp>

#include <sgpp/datadriven/DatadrivenOpFactory.hpp>
#include <sgpp/datadriven/operation/hash/simple/OperationDensityMargTo1D.hpp>
#include <sgpp/datadriven/operation/hash/simple/OperationDensityMarginalize.hpp>

#include <stddef.h>
#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <ctime>
#include <iostream>
#include <limits>
#include <string>
#include <utility>
#include <vector>

namespace sgpp {
namespace datadriven {

// --------------------------------------------------------------------------------------------
LearnerSGDEConfiguration::LearnerSGDEConfiguration() : json::JSON() { initConfig(); }

LearnerSGDEConfiguration::LearnerSGDEConfiguration(const std::string& fileName)
    : json::JSON(fileName) {
  initConfig();
  // initialize structs from file
  // configure grid
  try {
    if (this->contains("grid_filename")) gridConfig.filename_ = (*this)["grid_filename"].get();
    if (this->contains("grid_dim")) gridConfig.dim_ = (*this)["grid_level"].getUInt();
    if (this->contains("grid_level"))
      gridConfig.level_ = static_cast<int>((*this)["grid_level"].getInt());
    if (this->contains("grid_type"))
      gridConfig.type_ = stringToGridType((*this)["grid_type"].get());

    // configure adaptive refinement
    if (this->contains("refinement_numSteps"))
      adaptivityConfig.numRefinements_ = (*this)["refinement_numSteps"].getUInt();
    if (this->contains("refinement_numPoints"))
      adaptivityConfig.numRefinementPoints_ = (*this)["refinement_numPoints"].getUInt();

    // configure solver
    if (this->contains("solver_type"))
      solverConfig.type_ = stringToSolverType((*this)["solver_type"].get());
    if (this->contains("solver_maxIterations"))
      solverConfig.maxIterations_ = (*this)["solver_maxIterations"].getUInt();
    if (this->contains("solver_eps")) solverConfig.eps_ = (*this)["solver_eps"].getDouble();
    if (this->contains("solver_threshold"))
      solverConfig.threshold_ = (*this)["solver_threshold"].getDouble();

    // configure regularization
    if (this->contains("regularization_type"))
      regularizationConfig.type_ = stringToRegularizationType((*this)["regularization_type"].get());

    // configure learner
    if (this->contains("crossValidation_lambda"))
      crossValidationConfig.lambda_ = (*this)["crossValidation_lambda"].getDouble();
    if (this->contains("crossValidation_enable"))
      crossValidationConfig.enable_ = (*this)["crossValidation_enable"].getBool();
    if (this->contains("crossValidation_kfold"))
      crossValidationConfig.kfold_ = (*this)["crossValidation_kfold"].getUInt();
    if (this->contains("crossValidation_lambdaStart"))
      crossValidationConfig.lambdaStart_ = (*this)["crossValidation_lambdaStart"].getDouble();
    if (this->contains("crossValidation_lambdaEnd"))
      crossValidationConfig.lambdaEnd_ = (*this)["crossValidation_lambdaEnd"].getDouble();
    if (this->contains("crossValidation_lambdaSteps"))
      crossValidationConfig.lambdaSteps_ = (*this)["crossValidation_lambdaSteps"].getUInt();
    if (this->contains("crossValidation_logScale"))
      crossValidationConfig.logScale_ = (*this)["crossValidation_logScale"].getBool();
    if (this->contains("crossValidation_shuffle"))
      crossValidationConfig.shuffle_ = (*this)["crossValidation_shuffle"].getBool();
    if (this->contains("crossValidation_seed"))
      crossValidationConfig.seed_ = static_cast<int>((*this)["crossValidation_seed"].getInt());
    if (this->contains("crossValidation_silent"))
      crossValidationConfig.silent_ = (*this)["crossValidation_silent"].getBool();
  } catch (json::json_exception& e) {
    std::cout << e.what() << std::endl;
  }
}

void LearnerSGDEConfiguration::initConfig() {
  // set default config
  gridConfig.dim_ = 0;
  gridConfig.level_ = 6;
  gridConfig.type_ = base::GridType::Linear;
  gridConfig.maxDegree_ = 1;
  gridConfig.boundaryLevel_ = 0;

  // configure adaptive refinement
  adaptivityConfig.numRefinements_ = 0;
  adaptivityConfig.numRefinementPoints_ = 5;

  // configure solver
  solverConfig.type_ = solver::SLESolverType::CG;
  solverConfig.maxIterations_ = 1000;
  solverConfig.eps_ = 1e-10;
  solverConfig.threshold_ = 1e-14;

  // configure regularization
  regularizationConfig.type_ = datadriven::RegularizationType::Laplace;

  // configure learner
  crossValidationConfig.enable_ = true;
  crossValidationConfig.kfold_ = 5;
  crossValidationConfig.lambda_ = 1e-5;
  crossValidationConfig.lambdaStart_ = 1e-1;
  crossValidationConfig.lambdaEnd_ = 1e-10;
  crossValidationConfig.lambdaSteps_ = 5;
  crossValidationConfig.logScale_ = true;
  crossValidationConfig.shuffle_ = false;
  crossValidationConfig.seed_ = 1234567;
  crossValidationConfig.silent_ = true;
}

LearnerSGDEConfiguration* LearnerSGDEConfiguration::clone() {
  LearnerSGDEConfiguration* clone = new LearnerSGDEConfiguration(*this);
  return clone;
}

sgpp::base::GridType LearnerSGDEConfiguration::stringToGridType(std::string& gridType) {
  if (gridType.compare("Linear") == 0) {
    return sgpp::base::GridType::Linear;
  } else if (gridType.compare("LinearStretched") == 0) {
    return sgpp::base::GridType::LinearStretched;
  } else if (gridType.compare("LinearL0Boundary") == 0) {
    return sgpp::base::GridType::LinearL0Boundary;
  } else if (gridType.compare("LinearBoundary") == 0) {
    return sgpp::base::GridType::LinearBoundary;
  } else if (gridType.compare("LinearStretchedBoundary") == 0) {
    return sgpp::base::GridType::LinearStretchedBoundary;
  } else if (gridType.compare("LinearTruncatedBoundary") == 0) {
    return sgpp::base::GridType::LinearTruncatedBoundary;
  } else if (gridType.compare("ModLinear") == 0) {
    return sgpp::base::GridType::ModLinear;
  } else if (gridType.compare("Poly") == 0) {
    return sgpp::base::GridType::Poly;
  } else if (gridType.compare("PolyBoundary") == 0) {
    return sgpp::base::GridType::PolyBoundary;
  } else if (gridType.compare("ModPoly") == 0) {
    return sgpp::base::GridType::ModPoly;
  } else if (gridType.compare("ModWavelet") == 0) {
    return sgpp::base::GridType::ModWavelet;
  } else if (gridType.compare("ModBspline") == 0) {
    return sgpp::base::GridType::ModBspline;
  } else if (gridType.compare("Prewavelet") == 0) {
    return sgpp::base::GridType::Prewavelet;
  } else if (gridType.compare("SquareRoot") == 0) {
    return sgpp::base::GridType::SquareRoot;
  } else if (gridType.compare("Periodic") == 0) {
    return sgpp::base::GridType::Periodic;
  } else if (gridType.compare("LinearClenshawCurtis") == 0) {
    return sgpp::base::GridType::LinearClenshawCurtis;
  } else if (gridType.compare("Bspline") == 0) {
    return sgpp::base::GridType::Bspline;
  } else if (gridType.compare("BsplineBoundary") == 0) {
    return sgpp::base::GridType::BsplineBoundary;
  } else if (gridType.compare("BsplineClenshawCurtis") == 0) {
    return sgpp::base::GridType::BsplineClenshawCurtis;
  } else if (gridType.compare("Wavelet") == 0) {
    return sgpp::base::GridType::Wavelet;
  } else if (gridType.compare("WaveletBoundary") == 0) {
    return sgpp::base::GridType::WaveletBoundary;
  } else if (gridType.compare("FundamentalNakSplineBoundary") == 0) {
    return sgpp::base::GridType::FundamentalNakSplineBoundary;
  } else if (gridType.compare("FundamentalSpline") == 0) {
    return sgpp::base::GridType::FundamentalSpline;
  } else if (gridType.compare("FundamentalSplineBoundary") == 0) {
    return sgpp::base::GridType::FundamentalSplineBoundary;
  } else if (gridType.compare("ModFundamentalSpline") == 0) {
    return sgpp::base::GridType::ModFundamentalSpline;
  } else if (gridType.compare("ModBsplineClenshawCurtis") == 0) {
    return sgpp::base::GridType::ModBsplineClenshawCurtis;
  } else if (gridType.compare("LinearStencil") == 0) {
    return sgpp::base::GridType::LinearStencil;
  } else if (gridType.compare("ModLinearStencil") == 0) {
    return sgpp::base::GridType::ModLinearStencil;
  } else {
    throw sgpp::base::application_exception("grid type is unknown");
  }
}

sgpp::datadriven::RegularizationType LearnerSGDEConfiguration::stringToRegularizationType(
    std::string& regularizationType) {
  if (regularizationType.compare("Identity") == 0) {
    return sgpp::datadriven::RegularizationType::Identity;
  } else if (regularizationType.compare("Laplace") == 0) {
    return sgpp::datadriven::RegularizationType::Laplace;
  } else {
    throw sgpp::base::application_exception("regularization type is unknown");
  }
}

sgpp::solver::SLESolverType LearnerSGDEConfiguration::stringToSolverType(std::string& solverType) {
  if (solverType.compare("CG")) {
    return sgpp::solver::SLESolverType::CG;
  } else if (solverType.compare("BiCGSTAB")) {
    return sgpp::solver::SLESolverType::BiCGSTAB;
  } else {
    throw sgpp::base::application_exception("solver type is unknown");
  }
}

// --------------------------------------------------------------------------------------------
LearnerSGDE::LearnerSGDE(sgpp::base::RegularGridConfiguration& gridConfig,
                         sgpp::base::AdaptivityConfiguration& adaptivityConfig,
                         sgpp::solver::SLESolverConfiguration& solverConfig,
                         sgpp::datadriven::RegularizationConfiguration& regularizationConfig,
                         CrossvalidationConfiguration& crossValidationConfig)
    : error(0.0),
      grid(nullptr),
      alpha(nullptr),
      trainData(nullptr),
      trainLabels(nullptr),
      usePrior(false),
      lambdaReg(1e-6),
      gridConfig(gridConfig),
      adaptivityConfig(adaptivityConfig),
      solverConfig(solverConfig),
      regularizationConfig(regularizationConfig),
      crossValidationConfig(crossValidationConfig) {}

LearnerSGDE::LearnerSGDE(LearnerSGDEConfiguration& learnerSGDEConfig)
    : LearnerSGDE(learnerSGDEConfig.gridConfig, learnerSGDEConfig.adaptivityConfig,
                  learnerSGDEConfig.solverConfig, learnerSGDEConfig.regularizationConfig,
                  learnerSGDEConfig.crossValidationConfig) {}

LearnerSGDE::LearnerSGDE(const LearnerSGDE& learnerSGDE) {
  error = 0.0;
  grid = learnerSGDE.grid;
  alpha = learnerSGDE.alpha;
  trainData = learnerSGDE.trainData;
  usePrior = false;
  lambdaReg = 1e-6;
  gridConfig = learnerSGDE.gridConfig;
  adaptivityConfig = learnerSGDE.adaptivityConfig;
  solverConfig = learnerSGDE.solverConfig;
  regularizationConfig = learnerSGDE.regularizationConfig;
  crossValidationConfig = learnerSGDE.crossValidationConfig;
}

LearnerSGDE::~LearnerSGDE() {}

// -----------------------------------------------------------------------------------------------

void LearnerSGDE::initialize(base::DataMatrix& samples) {
  trainData = std::make_shared<base::DataMatrix>(samples);
  gridConfig.dim_ = samples.getNcols();

  grid = createRegularGrid();
  alpha = std::make_shared<base::DataVector>(grid->getSize());

  // optimize the regularization parameter if cv enabled
  if (crossValidationConfig.enable_) {
    lambdaReg = optimizeLambdaCV();
  } else {
    lambdaReg = crossValidationConfig.lambda_;
  }

  std::cout << "lambda: " << lambdaReg << std::endl;
}

// ---------------------------------------------------------------------------

double LearnerSGDE::pdf(base::DataVector& x) {
  return op_factory::createOperationEval(*grid)->eval(*alpha, x);
}

void LearnerSGDE::pdf(base::DataMatrix& points, base::DataVector& res) {
  op_factory::createOperationMultipleEval(*grid, points)->eval(*alpha, res);
}

double LearnerSGDE::mean(base::Grid& grid, base::DataVector& alpha) {
  return op_factory::createOperationFirstMoment(grid)->doQuadrature(alpha);
}

double LearnerSGDE::mean() { return mean(*grid, *alpha); }

double LearnerSGDE::variance(base::Grid& grid, base::DataVector& alpha) {
  double secondMoment = op_factory::createOperationSecondMoment(grid)->doQuadrature(alpha);

  // use Steiners translation theorem to compute the variance
  double firstMoment = mean();
  double res = secondMoment - firstMoment * firstMoment;
  return res;
}

double LearnerSGDE::variance() { return variance(*grid, *alpha); }

void LearnerSGDE::cov(base::DataMatrix& cov, base::DataMatrix* bounds) {
  std::unique_ptr<datadriven::OperationCovariance> opCov(
      op_factory::createOperationCovariance(*grid));
  opCov->doQuadrature(*alpha, cov, bounds);
}

std::shared_ptr<base::DataVector> LearnerSGDE::getSamples(size_t dim) {
  std::shared_ptr<base::DataVector> isamples = std::make_shared<base::DataVector>(getNsamples());
  trainData->getColumn(dim, *isamples);
  return isamples;
}

std::shared_ptr<base::DataMatrix> LearnerSGDE::getSamples() { return trainData; }

size_t LearnerSGDE::getDim() { return gridConfig.dim_; }

size_t LearnerSGDE::getNsamples() { return trainData->getNrows(); }

base::DataVector* LearnerSGDE::getSurpluses() { return alpha.get(); }

base::Grid* LearnerSGDE::getGrid() { return grid.get(); }

std::shared_ptr<base::DataVector> LearnerSGDE::getSharedSurpluses() { return alpha; }
std::shared_ptr<base::Grid> LearnerSGDE::getSharedGrid() { return grid; }

// ---------------------------------------------------------------------------

std::shared_ptr<base::Grid> LearnerSGDE::createRegularGrid() {
  // load grid
  std::unique_ptr<base::Grid> uGrid;
  if (gridConfig.filename_.length() > 0) {
    std::ifstream ifs(gridConfig.filename_);
    std::string content((std::istreambuf_iterator<char>(ifs)), (std::istreambuf_iterator<char>()));
    uGrid.reset(base::Grid::unserialize(content));
  } else {
    if (gridConfig.type_ == base::GridType::Linear) {
      uGrid.reset(base::Grid::createLinearGrid(gridConfig.dim_));
    } else if (gridConfig.type_ == base::GridType::LinearL0Boundary) {
      uGrid.reset(base::Grid::createLinearBoundaryGrid(gridConfig.dim_, 0));
    } else if (gridConfig.type_ == base::GridType::LinearBoundary) {
      uGrid.reset(base::Grid::createLinearBoundaryGrid(gridConfig.dim_, 1));
    } else {
      throw base::application_exception("LeanerSGDE::initialize : grid type is not supported");
    }
    uGrid->getGenerator().regular(gridConfig.level_);
  }

  // move the grid to be shared
  std::shared_ptr<base::Grid> sGrid{std::move(uGrid)};

  return sGrid;
}

double LearnerSGDE::optimizeLambdaCV() {
  double curLambda;
  double bestLambda = 0;
  double curMean = 0;
  double curMeanAcc = 0;
  double bestMeanAcc = 0;

  size_t kfold = crossValidationConfig.kfold_;

  std::vector<std::shared_ptr<base::DataMatrix>> kfold_train(kfold);
  std::vector<std::shared_ptr<base::DataMatrix>> kfold_test(kfold);
  splitset(kfold_train, kfold_test);

  double lambdaStart = crossValidationConfig.lambdaStart_;
  double lambdaEnd = crossValidationConfig.lambdaEnd_;

  if (crossValidationConfig.logScale_) {
    lambdaStart = std::log(lambdaStart);
    lambdaEnd = std::log(lambdaEnd);
  }

  for (size_t i = 0; i < crossValidationConfig.lambdaSteps_; i++) {
    // compute current lambda
    curLambda = lambdaStart +
                static_cast<double>(i) * (lambdaEnd - lambdaStart) /
                    static_cast<double>(crossValidationConfig.lambdaSteps_ - 1);

    if (crossValidationConfig.logScale_) curLambda = exp(curLambda);

    if (i % static_cast<size_t>(
                std::max(static_cast<double>(crossValidationConfig.lambdaSteps_) / 10.0,
                         static_cast<double>(1.0))) ==
        0) {
      if (!crossValidationConfig.silent_) {
        std::cout << i + 1 << "/" << crossValidationConfig.lambdaSteps_
                  << " (lambda = " << curLambda << ") " << std::endl;
        std::cout.flush();
      }
    }

    // cross-validation
    curMeanAcc = 0.0;
    curMean = 0.0;
    std::shared_ptr<base::Grid> grid;
    base::DataVector alpha(getNsamples());

    for (size_t j = 0; j < kfold; j++) {
      // initialize standard grid and alpha vector
      grid = createRegularGrid();
      alpha.resizeZero(grid->getSize());

      // compute density
      train(*grid, alpha, *(kfold_train[j]), curLambda);
      // get L2 norm of residual for test set
      curMean = computeResidual(*grid, alpha, *(kfold_test[j]), 0.0);
      curMeanAcc += curMean;

      if (!crossValidationConfig.silent_) {
        std::cout << "# " << curLambda << " " << i << " " << j << " " << curMeanAcc << " "
                  << curMean << "; alpha in [" << alpha.min() << ", " << alpha.max() << "]"
                  << "; data in " << kfold_test[j]->getNrows() << " x " << kfold_test[j]->getNcols()
                  << std::endl;
      }
    }

    curMeanAcc /= static_cast<double>(kfold);

    if (i == 0 || curMeanAcc < bestMeanAcc) {
      bestMeanAcc = curMeanAcc;
      bestLambda = curLambda;
    }

    if (!crossValidationConfig.silent_) {
      std::cout << "# " << curLambda << " " << bestLambda << " " << i << " " << curMeanAcc
                << std::endl;
    }
  }

  if (!crossValidationConfig.silent_) {
    std::cout << "# -> best lambda = " << bestLambda << std::endl;
  }

  return bestLambda;
}

void LearnerSGDE::train() {
  // learn the data -> do the density estimation
  train(*grid, *alpha, *trainData, lambdaReg);
}

void LearnerSGDE::train(base::Grid& grid, base::DataVector& alpha, base::DataMatrix& trainData,
                        double lambdaReg) {
  size_t dim = trainData.getNcols();

  base::GridStorage& gridStorage = grid.getStorage();
  base::GridGenerator& gridGen = grid.getGenerator();
  base::DataVector rhs(grid.getSize());
  alpha.resize(grid.getSize());
  alpha.setAll(0.0);

  if (!crossValidationConfig.silent_) {
    std::cout << "# LearnerSGDE: grid points " << grid.getSize() << std::endl;
  }

  for (size_t ref = 0; ref <= adaptivityConfig.numRefinements_; ref++) {
    auto C = computeRegularizationMatrix(grid);

    datadriven::DensitySystemMatrix SMatrix(grid, trainData, C, lambdaReg);
    SMatrix.generateb(rhs);

    if (!crossValidationConfig.silent_) {
      std::cout << "# LearnerSGDE: Solving " << std::endl;
    }

    solver::ConjugateGradients myCG(solverConfig.maxIterations_, solverConfig.eps_);
    myCG.solve(SMatrix, alpha, rhs, false, false, solverConfig.threshold_);

    if (myCG.getResiduum() > solverConfig.threshold_) {
      throw base::operation_exception("LearnerSGDE - train: conjugate gradients is not converged");
    }

    if (ref < adaptivityConfig.numRefinements_) {
      if (!crossValidationConfig.silent_) {
        std::cout << "# LearnerSGDE: Refine grid ... ";
      }

      // Weight surplus with function evaluation at grid points
      std::unique_ptr<base::OperationEval> opEval(op_factory::createOperationEval(grid));
      base::DataVector p(dim);
      base::DataVector alphaWeight(alpha.getSize());

      for (size_t i = 0; i < grid.getSize(); i++) {
        gridStorage.getPoint(i).getStandardCoordinates(p);
        alphaWeight[i] = alpha.get(i) * opEval->eval(alpha, p);
      }

      base::SurplusRefinementFunctor srf(alphaWeight, adaptivityConfig.numRefinementPoints_,
                                         adaptivityConfig.refinementThreshold_);
      gridGen.refine(srf);

      if (!crossValidationConfig.silent_) {
        std::cout << "# LearnerSGDE: ref " << ref << "/" << adaptivityConfig.numRefinements_ - 1
                  << ": " << grid.getSize() << std::endl;
      }

      alpha.resize(grid.getSize());
      rhs.resize(grid.getSize());
      alpha.setAll(0.0);
      rhs.setAll(0.0);
    }
  }

  return;
}

void LearnerSGDE::trainOnline(base::DataVector& labels, base::DataMatrix& testData,
                              base::DataVector& testLabels, base::DataMatrix* validData,
                              base::DataVector* validLabels, base::DataVector& classLabels,
                              size_t maxDataPasses, std::string refType, std::string refMonitor,
                              size_t refPeriod, double accDeclineThreshold,
                              size_t accDeclineBufferSize, size_t minRefInterval, bool usePrior) {
  this->trainLabels = std::make_shared<base::DataVector>(labels);
  this->usePrior = usePrior;

  size_t dim = trainData->getNcols();

  // initialize counter for dataset passes
  size_t cntDataPasses = 0;
  // counter for number of processed data points
  size_t processedPoints = 0;

  // initialize refinement variables
  double currentValidError = 0.0;
  double currentTrainError = 0.0;

  // create convergence monitor object
  RefinementMonitor* monitor = nullptr;
  if (refMonitor == "periodic") {
    monitor = new RefinementMonitorPeriodic(refPeriod);
  } else if (refMonitor == "convergence") {
    monitor =
        new RefinementMonitorConvergence(accDeclineThreshold, accDeclineBufferSize, minRefInterval);
  }

  // counts number of performed refinement steps
  size_t refCnt = 0;

  // create grid for each class
  for (size_t i = 0; i < classLabels.getSize(); i++) {
    int label = static_cast<int>(classLabels[i]);
    std::unique_ptr<base::Grid> uGrid;
    if (gridConfig.type_ == base::GridType::Linear) {
      uGrid.reset(base::Grid::createLinearGrid(gridConfig.dim_));
    } else if (gridConfig.type_ == base::GridType::ModLinear) {
      uGrid.reset(base::Grid::createModLinearGrid(gridConfig.dim_));
    } else {
      throw base::application_exception("LearnerSGDE::trainOnline : grid type is not supported");
    }

    uGrid->getGenerator().regular(gridConfig.level_);
    // move the grid to be shared
    std::shared_ptr<base::Grid> sGrid{std::move(uGrid)};
    // insert grid into grid collection
    grids.insert(std::pair<int, std::shared_ptr<base::Grid>>(label, sGrid));

    // create alpha vector for new label
    alphas.insert(std::pair<int, std::shared_ptr<base::DataVector>>(
        label, std::make_shared<base::DataVector>(sGrid->getSize())));
    alphas.at(label)->setAll(0.0);

    appearances.insert(std::pair<int, size_t>(label, 0));
    priors.insert(std::pair<int, double>(label, 0.0));
  }

  // auxiliary variable for accuracy (error) measurement
  double acc = 0.0;
  acc = getAccuracy(testData, testLabels, 0.0);
  avgErrors.append(1.0 - acc);

  // main loop which performs the training process
  while (cntDataPasses < maxDataPasses) {
    for (size_t i = 0; i < trainData->getNrows(); i++) {
      // get next training sample x and its label y
      sgpp::base::DataVector x(dim);
      trainData->getRow(static_cast<size_t>(i), x);
      double y = trainLabels->get(i);
      int label = static_cast<int>(y);

      appearances.at(label) += 1;

      grid = grids.at(label);
      alpha = alphas.at(label);

      // solve the system using CG to obtain new surplus vector
      base::DataVector newAlpha(*alpha);
      base::DataVector rhs(grid->getSize());
      base::DataMatrix dataSample(0, dim);
      dataSample.appendRow(x);
      auto C = computeRegularizationMatrix(*grid);
      datadriven::DensitySystemMatrix SMatrix(*grid, dataSample, C, lambdaReg);
      SMatrix.generateb(rhs);
      solver::ConjugateGradients myCG(solverConfig.maxIterations_, solverConfig.eps_);
      myCG.solve(SMatrix, newAlpha, rhs, false, false, solverConfig.threshold_);

      /*if (myCG.getResiduum() > solverConfig.threshold_) {
        throw base::operation_exception(
          "LearnerSGDE - train: conjugate gradients is not converged");
      }*/

      // apply weighting -> determine new surplus vector
      // corresponding to class label using new alpha vector and
      // the previously computed one -> use equal weighting depending
      // on the number of processed samples per class so far
      size_t alphaCnt = appearances.at(label);
      // reset
      alpha->mult(static_cast<double>(alphaCnt - 1));
      // determine new surplus vector
      alpha->add(newAlpha);
      alpha->mult(1.0 / static_cast<double>(alphaCnt));

      // check if refinement should be performed
      size_t refinementsNecessary = 0;
      if (refCnt < adaptivityConfig.numRefinements_ && processedPoints > 0 && monitor) {
        currentValidError = getError(*validData, *validLabels, 0.0, "Acc");
        // if train dataset is large use a subset for error evaluation
        currentTrainError = getError(*trainData, *trainLabels, 0.0, "Acc");
        monitor->pushToBuffer(1, currentValidError, currentTrainError);
        refinementsNecessary = monitor->refinementsNecessary();
      }

      // refinement
      while (refinementsNecessary > 0) {
        // acc = getAccuracy(testData, testLabels, 0.0);
        // avgErrors.append(1.0 - acc);
        std::cout << "Refinement at iteration: " << processedPoints + 1 << std::endl;
        // bundle grids and surplus vector pointer needed for refinement
        // (for zero-crossings refinement, data-based refinement)
        std::vector<sgpp::base::Grid*> refGrids;
        std::vector<sgpp::base::DataVector*> refAlphas;
        std::vector<double> refPriors;
        for (auto& g : grids) {
          refGrids.push_back(&*(g.second));
          refAlphas.push_back(&*(alphas.at(g.first)));
          if (usePrior) {
            refPriors.push_back(priors.at(g.first));
          } else {
            refPriors.push_back(1.0);
          }
        }
        bool levelPenalize = false;  // multiplies penalzing term for fine levels
        bool preCompute = true;      // precomputes and caches evals for zrcr
        sgpp::datadriven::MultiGridRefinementFunctor* func = nullptr;
        // Zero-crossing-based refinement
        sgpp::datadriven::ZeroCrossingRefinementFunctor funcZrcr(
            refGrids, refAlphas, refPriors, adaptivityConfig.numRefinementPoints_, levelPenalize,
            preCompute);
        // Data-based refinement. Needs a problem dependent coeffA. The values
        // can be determined by testing (aim at ~10 % of the training data is
        // to be marked relevant). Cross-validation or similar can/should be
        // employed
        // to determine this value.
        std::vector<double> coeffA;
        coeffA.push_back(1.2);  // ripley 1.2
        coeffA.push_back(1.2);  // ripley 1.2
        base::DataMatrix* refTrainData = trainData.get();
        base::DataVector* refTrainLabels = trainLabels.get();
        sgpp::datadriven::DataBasedRefinementFunctor funcData(
            refGrids, refAlphas, refPriors, refTrainData, refTrainLabels,
            adaptivityConfig.numRefinementPoints_, levelPenalize, coeffA);
        if (refType == "zero") {
          func = &funcZrcr;
        } else if (refType == "data") {
          func = &funcData;
        }

        // refine each grid
        int gIdx = 0;
        for (auto& g : grids) {
          grid = g.second;
          alpha = alphas.at(g.first);
          if (refType == "surplus") {
            // surplus refinement
            // weight surplus with function evaluation at grid points
            base::GridStorage& gridStorage = grid->getStorage();
            std::unique_ptr<base::OperationEval> opEval(op_factory::createOperationEval(*grid));
            base::DataVector p(dim);
            base::DataVector alphaWeight(alpha->getSize());
            for (size_t j = 0; j < grid->getSize(); j++) {
              gridStorage.getPoint(j).getStandardCoordinates(p);
              alphaWeight[j] = alpha->get(j) * opEval->eval(*alpha, p);
            }

            base::SurplusRefinementFunctor srf(alphaWeight, adaptivityConfig.numRefinementPoints_,
                                               adaptivityConfig.refinementThreshold_);
            // base::SurplusRefinementFunctor srf(
            //  *alpha, adaptivityConfig.numRefinementPoints_,
            //  adaptivityConfig.threshold_);
            // refine grid
            grid->getGenerator().refine(srf);
          } else if ((refType == "data") || (refType == "zero")) {
            if (preCompute) {
              // precompute the evals; needs to be done once per step, before
              // any refinement is done
              func->preComputeEvaluations();
            }

            func->setGridIndex(gIdx);

            grid->getGenerator().refine(*func);
          }
          std::cout << "# LearnerSGDE (class " << gIdx << "): ref " << refCnt + 1 << "/"
                    << adaptivityConfig.numRefinements_ << " new grid size: " << grid->getSize()
                    << std::endl;

          // keep computed alpha values and append zeros only for new points
          alpha->resizeZero(grid->getSize());
          gIdx++;
        }

        refCnt++;
        refinementsNecessary--;
      }

      // update prior probabilities
      priors[label] = ((priors[label] * static_cast<double>(i)) + 1) / (1 + static_cast<double>(i));

      // save current error
      if ((processedPoints > 0) && ((processedPoints + 1) % 10 == 0)) {
        acc = getAccuracy(testData, testLabels, 0.0);
        avgErrors.append(1.0 - acc);
      }

      processedPoints++;
    }

    cntDataPasses++;
  }
  std::cout << "# Training finished" << std::endl;
  error = 1.0 - getAccuracy(testData, testLabels, 0.0);
}

void LearnerSGDE::storeResults(base::DataMatrix& testDataset) {
  base::DataVector predictedLabels(testDataset.getNrows());
  predict(testDataset, predictedLabels);

  // write predicted classes to csv file
  std::ofstream output;
  output.open("SGDE_predicted_classes.csv");
  if (output.fail()) {
    std::cout << "failed to create csv file!" << std::endl;
  } else {
    for (size_t i = 0; i < predictedLabels.getSize(); i++) {
      base::DataVector x(2);
      testDataset.getRow(static_cast<size_t>(i), x);
      output << x[0] << ";" << x[1] << ";" << predictedLabels[i] << std::endl;
    }
    output.close();
  }
  // write grids to csv file
  output.open("SGDE_grid_1.csv");
  if (output.fail()) {
    std::cout << "failed to create csv file!" << std::endl;
  } else {
    grid = grids.at(-1);
    base::GridStorage& storage = grid->getStorage();
    base::GridStorage::grid_map_iterator end_iter = storage.end();
    for (base::GridStorage::grid_map_iterator iter = storage.begin(); iter != end_iter; iter++) {
      base::DataVector gpCoord(testDataset.getNcols());
      storage.getCoordinates(*(iter->first), gpCoord);
      for (size_t d = 0; d < gpCoord.getSize(); d++) {
        if (d < gpCoord.getSize() - 1) {
          output << gpCoord[d] << ";";
        } else {
          output << gpCoord[d] << std::endl;
        }
      }
    }
    output.close();
  }
  output.open("SGDE_grid_-1.csv");
  if (output.fail()) {
    std::cout << "failed to create csv file!" << std::endl;
  } else {
    grid = grids.at(1);
    base::GridStorage& storage = grid->getStorage();
    base::GridStorage::grid_map_iterator end_iter = storage.end();
    for (base::GridStorage::grid_map_iterator iter = storage.begin(); iter != end_iter; iter++) {
      base::DataVector gpCoord(testDataset.getNcols());
      storage.getCoordinates(*(iter->first), gpCoord);
      for (size_t d = 0; d < gpCoord.getSize(); d++) {
        if (d < gpCoord.getSize() - 1) {
          output << gpCoord[d] << ";";
        } else {
          output << gpCoord[d] << std::endl;
        }
      }
    }
    output.close();
  }

  // write density function evaluations to csv file
  double stepSize = 0.01;
  base::DataMatrix values(0, 2);
  base::DataVector range(101);
  for (size_t i = 0; i < 101; i++) {
    range.set(i, stepSize * (static_cast<double>(i)));
  }
  for (size_t i = 0; i < range.getSize(); i++) {
    for (size_t j = 0; j < range.getSize(); j++) {
      base::DataVector row(2);
      row.set(1, range.get(i));
      row.set(0, range.get(j));
      values.appendRow(row);
    }
  }
  // evaluate each density function at all points from 'values'
  // and write result to csv file
  for (auto const& g : grids) {
    output.open("SGDE_density_fun_" + std::to_string(g.first) + "_evals.csv");
    std::unique_ptr<base::OperationEval> opEval(op_factory::createOperationEval(*g.second));
    for (size_t i = 0; i < values.getNrows(); i++) {
      // Get next test sample x
      base::DataVector x(2);
      values.getRow(i, x);
      double res = opEval->eval(*alphas.at(g.first), x);
      output << res << ";" << std::endl;
    }
    output.close();
  }
}

double LearnerSGDE::getAccuracy(base::DataMatrix& testDataset,
                                const base::DataVector& referenceLabels, const double threshold) {
  // evaluate test dataset
  base::DataVector predictedLabels(testDataset.getNrows());
  predict(testDataset, predictedLabels);

  return getAccuracy(referenceLabels, threshold, predictedLabels);
}

double LearnerSGDE::getAccuracy(const base::DataVector& referenceLabels, const double threshold,
                                const base::DataVector& predictedLabels) {
  double result = -1.0;

  if (predictedLabels.getSize() != referenceLabels.getSize()) {
    throw base::application_exception(
        "LearnerSGDE::getAccuracy: lengths of classes vectors do not match!");
  }

  size_t correct = 0;

  for (size_t i = 0; i < predictedLabels.getSize(); i++) {
    if (predictedLabels.get(i) == referenceLabels.get(i)) correct++;
  }

  result = static_cast<double>(correct) / static_cast<double>(predictedLabels.getSize());

  return result;
}

void LearnerSGDE::predict(base::DataMatrix& testData, base::DataVector& predictedLabels) {
  predictedLabels.resize(testData.getNrows());
  size_t dim = testData.getNcols();

  double prior;

  for (size_t i = 0; i < testData.getNrows(); i++) {
    // get next test sample x
    base::DataVector x(dim);
    testData.getRow(static_cast<size_t>(i), x);
    // predict label using Bayes' Theorem
    double max = std::numeric_limits<double>::max() * (-1);
    int predLabel = 0;
    // compute each density function for current test sample x
    for (auto const& g : grids) {
      std::unique_ptr<base::OperationEval> opEval(op_factory::createOperationEval(*g.second));
      double res = opEval->eval(*alphas.at(g.first), x);
      // determine prior
      if (usePrior) {
        prior = priors.at(g.first);
      } else {
        prior = 1.0;
      }
      // multiply by prior
      res *= prior;
      // check if this class probability is larger than current maximum
      if (res > max) {
        max = res;
        predLabel = g.first;
      }
    }
    predictedLabels.set(i, predLabel);
  }
}

double LearnerSGDE::getError(base::DataMatrix& data, const base::DataVector& labels,
                             const double threshold, std::string errorType) {
  double res = -1.0;

  if (errorType == "Acc") {
    res = 1.0 - getAccuracy(data, labels, threshold);
  }

  return res;
}

double LearnerSGDE::computeResidual(base::Grid& grid, base::DataVector& alpha,
                                    base::DataMatrix& test, double lambdaReg) {
  auto C = computeRegularizationMatrix(grid);

  base::DataVector rhs(grid.getSize());
  base::DataVector res(grid.getSize());
  datadriven::DensitySystemMatrix SMatrix(grid, test, C, lambdaReg);
  SMatrix.generateb(rhs);

  SMatrix.mult(alpha, res);

  for (size_t i = 0; i < res.getSize(); i++) {
    res[i] = res[i] - rhs[i];
  }
  return res.l2Norm();
}

base::OperationMatrix* LearnerSGDE::computeRegularizationMatrix(base::Grid& grid) {
  base::OperationMatrix* C;

  if (regularizationConfig.type_ == datadriven::RegularizationType::Identity) {
    C = op_factory::createOperationIdentity(grid);
  } else if (regularizationConfig.type_ == datadriven::RegularizationType::Laplace) {
    C = op_factory::createOperationLaplace(grid);
  } else {
    throw base::application_exception("LearnerSGDE::train : unknown regularization type");
  }

  return C;
}

void LearnerSGDE::splitset(std::vector<std::shared_ptr<base::DataMatrix>>& strain,
                           std::vector<std::shared_ptr<base::DataMatrix>>& stest) {
  std::shared_ptr<base::DataMatrix> mydata = std::make_shared<base::DataMatrix>(*trainData);
  base::DataVector p(trainData->getNcols());
  base::DataVector tmp(trainData->getNcols());

  size_t kfold = crossValidationConfig.kfold_;

  std::vector<size_t> s(kfold);        // size of partition
  std::vector<size_t> ind(kfold + 1);  // index of partition
  size_t n = mydata->getNrows();       // size of data

  if (crossValidationConfig.shuffle_) {
    if (crossValidationConfig.seed_ == -1)
      srand(static_cast<unsigned int>(time(nullptr)));
    else
      srand(crossValidationConfig.seed_);

    for (size_t i = 0; i < mydata->getNrows(); i++) {
      size_t r = i + (static_cast<size_t>(rand()) % (mydata->getNrows() - i));
      mydata->getRow(i, p);
      mydata->getRow(r, tmp);
      mydata->setRow(r, p);
      mydata->setRow(i, tmp);
    }
  }

  // set size of partitions
  if (!crossValidationConfig.silent_) std::cout << "# kfold: ";

  ind[0] = 0;

  for (size_t i = 0; i < kfold - 1; i++) {
    s[i] = n / kfold;
    ind[i + 1] = ind[i] + s[i];

    if (!crossValidationConfig.silent_) std::cout << s[i] << " ";
  }

  ind[kfold] = n;
  s[kfold - 1] = n - (kfold - 1) * (n / kfold);

  if (!crossValidationConfig.silent_) std::cout << s[kfold - 1] << std::endl;

  if (!crossValidationConfig.silent_) {
    std::cout << "# kfold ind: ";

    for (size_t i = 0; i <= kfold; i++) std::cout << ind[i] << " ";

    std::cout << std::endl;
  }

  // fill data
  for (size_t i = 0; i < kfold; i++) {
    // allocate memory
    strain[i] = std::make_shared<base::DataMatrix>(mydata->getNrows() - s[i], mydata->getNcols());
    stest[i] = std::make_shared<base::DataMatrix>(s[i], mydata->getNcols());

    size_t local_test = 0;
    size_t local_train = 0;

    for (size_t j = 0; j < mydata->getNrows(); j++) {
      mydata->getRow(j, p);

      if (ind[i] <= j && j < ind[i + 1]) {
        stest[i]->setRow(local_test, p);
        local_test++;
      } else {
        strain[i]->setRow(local_train, p);
        local_train++;
      }
    }
  }
}

}  // namespace datadriven
}  // namespace sgpp
