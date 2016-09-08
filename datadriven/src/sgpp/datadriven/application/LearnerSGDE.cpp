// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <sgpp/datadriven/application/LearnerSGDE.hpp>

#include <sgpp/base/exception/application_exception.hpp>
#include <sgpp/base/grid/generation/functors/SurplusRefinementFunctor.hpp>
#include <sgpp/datadriven/functors/MultiGridRefinementFunctor.hpp>
#include <sgpp/datadriven/functors/classification/DataBasedRefinementFunctor.hpp>
#include <sgpp/datadriven/functors/classification/ZeroCrossingRefinementFunctor.hpp>
#include <sgpp/base/grid/generation/GridGenerator.hpp>
#include <sgpp/base/grid/GridStorage.hpp>
#include <sgpp/base/grid/storage/hashmap/HashGridStorage.hpp>
#include <sgpp/base/operation/BaseOpFactory.hpp>
#include <sgpp/base/operation/hash/OperationEval.hpp>
#include <sgpp/base/operation/hash/OperationFirstMoment.hpp>
#include <sgpp/base/operation/hash/OperationMultipleEval.hpp>
#include <sgpp/pde/operation/PdeOpFactory.hpp>
#include <sgpp/solver/sle/ConjugateGradients.hpp>
#include <sgpp/datadriven/algorithm/DensitySystemMatrix.hpp>
#include <sgpp/solver/TypesSolver.hpp>
#include <sgpp/base/tools/json/json_exception.hpp>
#include <sgpp/base/exception/data_exception.hpp>

#include <sgpp/datadriven/DatadrivenOpFactory.hpp>
#include <sgpp/datadriven/operation/hash/simple/OperationDensityMargTo1D.hpp>
#include <sgpp/datadriven/operation/hash/simple/OperationDensityMarginalize.hpp>

#include <stddef.h>
#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <ctime>
#include <iostream>
#include <vector>
#include <string>

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
      adaptivityConfig.noPoints_ = (*this)["refinement_numPoints"].getUInt();

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
      regularizationConfig.regType_ =
          stringToRegularizationType((*this)["regularization_type"].get());

    // configure learner
    if (this->contains("crossValidation_lambda"))
      crossvalidationConfig.lambda_ = (*this)["crossValidation_lambda"].getDouble();
    if (this->contains("crossValidation_enable"))
      crossvalidationConfig.enable_ = (*this)["crossValidation_enable"].getBool();
    if (this->contains("crossValidation_kfold"))
      crossvalidationConfig.kfold_ = (*this)["crossValidation_kfold"].getUInt();
    if (this->contains("crossValidation_lambdaStart"))
      crossvalidationConfig.lambdaStart_ = (*this)["crossValidation_lambdaStart"].getDouble();
    if (this->contains("crossValidation_lambdaEnd"))
      crossvalidationConfig.lambdaEnd_ = (*this)["crossValidation_lambdaEnd"].getDouble();
    if (this->contains("crossValidation_lambdaSteps"))
      crossvalidationConfig.lambdaSteps_ = (*this)["crossValidation_lambdaSteps"].getUInt();
    if (this->contains("crossValidation_logScale"))
      crossvalidationConfig.logScale_ = (*this)["crossValidation_logScale"].getBool();
    if (this->contains("crossValidation_shuffle"))
      crossvalidationConfig.shuffle_ = (*this)["crossValidation_shuffle"].getBool();
    if (this->contains("crossValidation_seed"))
      crossvalidationConfig.seed_ = static_cast<int>((*this)["crossValidation_seed"].getInt());
    if (this->contains("crossValidation_silent"))
      crossvalidationConfig.silent_ = (*this)["crossValidation_silent"].getBool();
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
  adaptivityConfig.noPoints_ = 5;

  // configure solver
  solverConfig.type_ = solver::SLESolverType::CG;
  solverConfig.maxIterations_ = 1000;
  solverConfig.eps_ = 1e-10;
  solverConfig.threshold_ = 1e-14;

  // configure regularization
  regularizationConfig.regType_ = datadriven::RegularizationType::Laplace;

  // configure learner
  crossvalidationConfig.enable_ = true;
  crossvalidationConfig.kfold_ = 5;
  crossvalidationConfig.lambda_ = 1e-5;
  crossvalidationConfig.lambdaStart_ = 1e-1;
  crossvalidationConfig.lambdaEnd_ = 1e-10;
  crossvalidationConfig.lambdaSteps_ = 5;
  crossvalidationConfig.logScale_ = true;
  crossvalidationConfig.shuffle_ = false;
  crossvalidationConfig.seed_ = 1234567;
  crossvalidationConfig.silent_ = true;
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
  } else if (gridType.compare("FundamentalSpline") == 0) {
    return sgpp::base::GridType::FundamentalSpline;
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
                         sgpp::base::AdpativityConfiguration& adaptivityConfig,
                         sgpp::solver::SLESolverConfiguration& solverConfig,
                         sgpp::datadriven::RegularizationConfiguration& regularizationConfig,
                         CrossvalidationForRegularizationConfiguration& crossvalidationConfig)
    : grid(nullptr),
      alpha(nullptr),
      trainData(nullptr),
      labels(nullptr),
      gridConfig(gridConfig),
      adaptivityConfig(adaptivityConfig),
      solverConfig(solverConfig),
      regularizationConfig(regularizationConfig),
      crossvalidationConfig(crossvalidationConfig) {}

LearnerSGDE::LearnerSGDE(LearnerSGDEConfiguration& learnerSGDEConfig)
    : LearnerSGDE(learnerSGDEConfig.gridConfig, learnerSGDEConfig.adaptivityConfig,
                  learnerSGDEConfig.solverConfig, learnerSGDEConfig.regularizationConfig,
                  learnerSGDEConfig.crossvalidationConfig) {}

LearnerSGDE::LearnerSGDE(const LearnerSGDE& learnerSGDE) {
  grid = learnerSGDE.grid;
  alpha = learnerSGDE.alpha;
  trainData = learnerSGDE.trainData;
  gridConfig = learnerSGDE.gridConfig;
  adaptivityConfig = learnerSGDE.adaptivityConfig;
  solverConfig = learnerSGDE.solverConfig;
  regularizationConfig = learnerSGDE.regularizationConfig;
  crossvalidationConfig = learnerSGDE.crossvalidationConfig;
}

LearnerSGDE::~LearnerSGDE() {}

// -----------------------------------------------------------------------------------------------

void LearnerSGDE::initialize(base::DataMatrix& samples) {
  trainData = std::make_shared<base::DataMatrix>(samples);
  gridConfig.dim_ = samples.getNcols();

  grid = createRegularGrid();
  alpha = std::make_shared<base::DataVector>(grid->getSize());

  // optimize the regularization parameter
  //double lambdaReg = 0.0; // define as member of LearnerSGDE

  if (crossvalidationConfig.enable_) {
    lambdaReg = optimizeLambdaCV();
  } else {
    lambdaReg = crossvalidationConfig.lambda_;
  }
  
  //lambdaReg *= 5000.0;
  std::cout << "lambda: " << lambdaReg << std::endl;

  // batch mode
  // learn the data -> do the density estimation
  //train(*grid, *alpha, *trainData, lambdaReg);
  
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

void LearnerSGDE::cov(base::DataMatrix& cov) {
  size_t ndim = grid->getStorage().getDimension();

  if ((cov.getNrows() != ndim) || (cov.getNcols() != ndim)) {
    // covariance matrix has wrong size -> resize
    cov.resize(ndim, ndim);
  }

  // prepare covariance marix
  cov.setAll(0.0);

  // generate 1d densities and compute means and variances
  base::DataVector means(ndim);
  base::DataVector variances(ndim);

  std::unique_ptr<datadriven::OperationDensityMargTo1D> opMarg =
      op_factory::createOperationDensityMargTo1D(*grid);

  base::Grid* marginalizedGrid = NULL;
  base::DataVector* marginalizedAlpha = new base::DataVector(0);

  for (size_t idim = 0; idim < ndim; idim++) {
    opMarg->margToDimX(&*alpha, marginalizedGrid, marginalizedAlpha, idim);
    // store moments
    means[idim] = mean(*marginalizedGrid, *marginalizedAlpha);
    variances[idim] = variance(*marginalizedGrid, *marginalizedAlpha);

    delete marginalizedGrid;
  }

  // helper variables
  std::vector<size_t> mdims(2);
  double covij = 0.0;

  for (size_t idim = 0; idim < ndim; idim++) {
    // diagonal is equal to the variance of the marginalized densities
    cov.set(idim, idim, variances[idim]);

    for (size_t jdim = idim + 1; jdim < ndim; jdim++) {
      // marginalize the density
      mdims[0] = idim;
      mdims[1] = jdim;
      opMarg->margToDimXs(&*alpha, marginalizedGrid, marginalizedAlpha, mdims);
      // -----------------------------------------------------
      // compute the covariance of Cov(X_i, X_j)
      covij = mean(*marginalizedGrid, *marginalizedAlpha) - means[idim] * means[jdim];
      cov.set(idim, jdim, covij);
      cov.set(jdim, idim, covij);
      // -----------------------------------------------------
      delete marginalizedGrid;
    }
  }

  delete marginalizedAlpha;
}

std::shared_ptr<base::DataVector> LearnerSGDE::getSamples(size_t dim) {
  std::shared_ptr<base::DataVector> isamples = std::make_shared<base::DataVector>(getNsamples());
  trainData->getColumn(dim, *isamples);
  return isamples;
}

std::shared_ptr<base::DataMatrix> LearnerSGDE::getSamples() { return trainData; }

size_t LearnerSGDE::getDim() { return gridConfig.dim_; }

size_t LearnerSGDE::getNsamples() { return trainData->getNrows(); }

std::shared_ptr<base::DataVector> LearnerSGDE::getSurpluses() { return alpha; }

std::shared_ptr<base::Grid> LearnerSGDE::getGrid() { return grid; }

// ---------------------------------------------------------------------------

std::shared_ptr<base::Grid> LearnerSGDE::createRegularGrid() {
  // load grid
  std::unique_ptr<base::Grid> uGrid;
  if (gridConfig.filename_.length() > 0) {
    std::ifstream ifs(gridConfig.filename_);
    std::string content((std::istreambuf_iterator<char>(ifs)), (std::istreambuf_iterator<char>()));
    uGrid = base::Grid::unserialize(content);
  } else {
    if (gridConfig.type_ == base::GridType::Linear) {
      uGrid = base::Grid::createLinearGrid(gridConfig.dim_);
    } else if (gridConfig.type_ == base::GridType::LinearL0Boundary) {
      uGrid = base::Grid::createLinearBoundaryGrid(gridConfig.dim_, 0);
    } else if (gridConfig.type_ == base::GridType::LinearBoundary) {
      uGrid = base::Grid::createLinearBoundaryGrid(gridConfig.dim_, 1);
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

  size_t kfold = crossvalidationConfig.kfold_;

  std::vector<std::shared_ptr<base::DataMatrix> > kfold_train(kfold);
  std::vector<std::shared_ptr<base::DataMatrix> > kfold_test(kfold);
  splitset(kfold_train, kfold_test);

  double lambdaStart = crossvalidationConfig.lambdaStart_;
  double lambdaEnd = crossvalidationConfig.lambdaEnd_;

  if (crossvalidationConfig.logScale_) {
    lambdaStart = std::log(lambdaStart);
    lambdaEnd = std::log(lambdaEnd);
  }

  for (size_t i = 0; i < crossvalidationConfig.lambdaSteps_; i++) {
    // compute current lambda
    curLambda = lambdaStart +
                static_cast<double>(i) * (lambdaEnd - lambdaStart) /
                    static_cast<double>(crossvalidationConfig.lambdaSteps_ - 1);

    if (crossvalidationConfig.logScale_) curLambda = exp(curLambda);

    if (i % static_cast<size_t>(
                std::max(static_cast<double>(crossvalidationConfig.lambdaSteps_) / 10.0f,
                         static_cast<double>(1.0f))) ==
        0) {
      if (!crossvalidationConfig.silent_) {
        std::cout << i + 1 << "/" << crossvalidationConfig.lambdaSteps_
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

      if (!crossvalidationConfig.silent_) {
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

    if (!crossvalidationConfig.silent_) {
      std::cout << "# " << curLambda << " " << bestLambda << " " << i << " " << curMeanAcc
                << std::endl;
    }
  }

  if (!crossvalidationConfig.silent_) {
    std::cout << "# -> best lambda = " << bestLambda << std::endl;
  }

  return bestLambda;
}

void LearnerSGDE::train() {
  // batch mode
  // learn the data -> do the density estimation
  train(*grid, *alpha, *trainData, lambdaReg);
}

void LearnerSGDE::train(base::Grid& grid, base::DataVector& alpha, base::DataMatrix& trainData,
                        double lambdaReg) {
//void LearnerSGDE::train() {
  size_t dim = trainData.getNcols();

  base::GridStorage& gridStorage = grid.getStorage();
  base::GridGenerator& gridGen = grid.getGenerator();
  base::DataVector rhs(grid.getSize());
  alpha.resize(grid.getSize());
  alpha.setAll(0.0);

  if (!crossvalidationConfig.silent_) {
    std::cout << "# LearnerSGDE: grid points " << grid.getSize() << std::endl;
  }

  for (size_t ref = 0; ref <= adaptivityConfig.numRefinements_; ref++) {
    std::unique_ptr<base::OperationMatrix> C = computeRegularizationMatrix(grid);

    datadriven::DensitySystemMatrix SMatrix(grid, trainData, *C, lambdaReg);
    SMatrix.generateb(rhs);

    if (!crossvalidationConfig.silent_) {
      std::cout << "# LearnerSGDE: Solving " << std::endl;
    }

    solver::ConjugateGradients myCG(solverConfig.maxIterations_, solverConfig.eps_);
    myCG.solve(SMatrix, alpha, rhs, false, false, solverConfig.threshold_);

    if (myCG.getResiduum() > solverConfig.threshold_) {
      throw base::operation_exception("LearnerSGDE - train: conjugate gradients is not converged");
    }

    if (ref < adaptivityConfig.numRefinements_) {
      if (!crossvalidationConfig.silent_) {
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

      base::SurplusRefinementFunctor srf(alphaWeight, adaptivityConfig.noPoints_,
                                         adaptivityConfig.threshold_);
      gridGen.refine(srf);

      if (!crossvalidationConfig.silent_) {
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

void LearnerSGDE::trainOnline(base::DataVector& plabels,
                              base::DataMatrix& ptestData, base::DataVector& ptestLabels) {
  labels = std::make_shared<base::DataVector>(plabels);
  size_t dim = trainData->getNcols();

  /*base::GridStorage& gridStorage = grid->getStorage();
  base::GridGenerator& gridGen = grid->getGenerator();
  base::DataVector rhs(grid->getSize());
  alpha->resize(grid->getSize());
  alpha->setAll(0.0);

  if (!crossvalidationConfig.silent_) {
    std::cout << "# LearnerSGDE: grid points " << grid->getSize() << std::endl;
  }*/
  
  size_t maxIterations = 1; //ToDo: pass as parameter
  size_t numIterations = 0;

  // refinement parameters - ToDo: pass to Learner as parameters (adaptivityConfig)
  size_t numRefSteps = adaptivityConfig.numRefinements_;
  size_t refPeriod = 10; //Ripley 25 - Banana 100 / 200
  size_t numPoints = adaptivityConfig.noPoints_;
  double threshold = adaptivityConfig.threshold_;

  // map number of already performed refinements to grids
  std::map<int, size_t> refCnts;
  
  bool doRefine = false;      // set true by monitor to trigger refinement
  size_t refSteps = 0;

  // for error plotting
  base::DataVector error;
  //base::DataVector gridSizes;
  // for plotting only
  //double acc = getAccuracy(*trainData, *labels, 0.0);
  //double acc = getAccuracy(ptestData, ptestLabels, 0.0);
  //error.append(1.0 - acc);

  while (numIterations < maxIterations) {			
    for (size_t i = 0; i < trainData->getNrows(); i++) {
      //std::cout << "iteration: " << i << std::endl;
      // Get next training sample x and its label y
      sgpp::base::DataVector x(dim);					
      trainData->getRow((size_t)i, x);
      double y = labels->get(i);
      int label = static_cast<int>(y);

      // check class label -> if it appears for the
      // first time, create new grid 
      if (grids.find(label) == grids.end()) {
        std::unique_ptr<base::Grid> uGrid;
        if (gridConfig.type_ == base::GridType::Linear) {
          uGrid = base::Grid::createLinearGrid(gridConfig.dim_);
        } 
        else if (gridConfig.type_ == base::GridType::ModLinear) {
          uGrid = base::Grid::createModLinearGrid(gridConfig.dim_);
        } 
        else {
          throw base::application_exception("LearnerSGDE::trainOnline : grid type is not supported");
        }

        uGrid->getGenerator().regular(gridConfig.level_);
        // move the grid to be shared
        std::shared_ptr<base::Grid> sGrid{std::move(uGrid)};
        // insert grid into grid collection
        grids.insert(std::pair<int, std::shared_ptr<base::Grid>>(label, sGrid));

        // create alpha vector for new label
        alphas.insert(std::pair<int, std::shared_ptr<base::DataVector>>(label, 
                      std::make_shared<base::DataVector>(sGrid->getSize())));
        alphas.at(label)->setAll(0.0);
        // create new queue in alpha storage to collect all 
        // upcoming alpha vectors corresponding to the new label 
        //alphaStorage.insert(std::make_pair(label, 
        //                    std::deque<std::shared_ptr<base::DataVector>>));
        
        appearances.insert(std::pair<int, size_t>(label, 0));
        refCnts.insert(std::pair<int, size_t>(label, 0));
      }
      appearances.at(label) += 1;

      grid = grids.at(label); // &*grids.at(label)
      alpha = alphas.at(label); // &*alphas.at(label)
      base::DataVector newAlpha(*alpha);
      //base::DataVector newAlpha(grid->getSize());
      //newAlpha.setAll(0.0);
      base::GridStorage& gridStorage = grid->getStorage();
      base::GridGenerator& gridGen = grid->getGenerator();
      base::DataVector rhs(grid->getSize());
      //base::GridStorage& gridStorage = grids.at(label)->getStorage();
      //base::GridGenerator& gridGen = grids.at(label)->getGenerator();
      //base::DataVector rhs(grids.at(label)->getSize());
      
      //std::shared_ptr<base::DataMatrix> dataSample = new base::DataMatrix(0, dim);
      //dataSample->append(x);
      base::DataMatrix dataSample(0, dim);
      dataSample.appendRow(x);
      std::unique_ptr<base::OperationMatrix> C = computeRegularizationMatrix(*grid);
      datadriven::DensitySystemMatrix SMatrix(*grid, dataSample, *C, lambdaReg);
      SMatrix.generateb(rhs);

      solver::ConjugateGradients myCG(solverConfig.maxIterations_, solverConfig.eps_);
      myCG.solve(SMatrix, newAlpha, rhs, false, false, solverConfig.threshold_);

      if (myCG.getResiduum() > solverConfig.threshold_) {
        throw base::operation_exception("LearnerSGDE - train: conjugate gradients is not converged");
      }

      // store new alpha vector
      //alphaStorage.at(label).push_back(alpha);
      // apply weighting -> determine new surplus vector 
      // corresponding to class label using new alpha vector and
      // the previously computed one -> use equal weighting depending
      // on the number of processed samples per class so far
      size_t alphaCnt = appearances.at(label);
      // reset
      alpha->mult(static_cast<double>(alphaCnt-1));
      // determine new surplus vector
      alpha->add(newAlpha);
      alpha->mult(1.0/static_cast<double>(alphaCnt));
      
      //alphas.at(label)->copyFrom(*alpha);

      //if ( (i > 0) && (i % refPeriod == 0) ) {
      if ( (appearances.at(label) > 0) && (appearances.at(label) % refPeriod == 0) ) {
        doRefine = true;
      }
      // Refinement
      if ( (refCnts.at(label) < numRefSteps) && doRefine ) {
        // Surplus refinement 
        // Weight surplus with function evaluation at grid points
        std::unique_ptr<base::OperationEval> opEval(op_factory::createOperationEval(*grid));
        base::DataVector p(dim);
        base::DataVector alphaWeight(alpha->getSize());
        for (size_t j = 0; j < grid->getSize(); j++) {
          gridStorage.getPoint(j).getStandardCoordinates(p);
          alphaWeight[j] = alpha->get(j) * opEval->eval(*alpha, p);
        }

        base::SurplusRefinementFunctor srf(alphaWeight, numPoints, threshold);
        //base::SurplusRefinementFunctor srf(*alpha, numPoints, threshold);
        //gridGen.refine(srf);
        
        // Impurity refinement
        // ToDo:

        // zero-crossing / data-based refinement
        std::vector<sgpp::base::Grid*> ref_grids;
        std::vector<sgpp::base::DataVector*> ref_alphas;
        //for (auto const& g : grids) {
        ref_grids.push_back(&*grids.at(-1));
        ref_grids.push_back(&*grids.at(1));
        ref_alphas.push_back(&*alphas.at(-1));
        ref_alphas.push_back(&*alphas.at(1));

        bool levelPenalize = false;  // Multiplies penalzing term for fine levels
        bool preCompute = true;      // Precomputes and caches evals for zrcr & grid
        sgpp::datadriven::MultiGridRefinementFunctor* func = nullptr;
        // Zero-crossing-based refinement
        sgpp::datadriven::ZeroCrossingRefinementFunctor funcZrcr(ref_grids, ref_alphas,
                                                                 numPoints,
                                                                 levelPenalize,
                                                                 preCompute);
        // Data-based refinement. Needs a problem dependent coeffA. The values
        // were determined by testing (aim at ~10 % of the training data is
        // to be marked relevant. Cross-validation or similar can/should be employed
        // to determine this value.
        std::vector<double> coeffA;
        coeffA.push_back(1.2);
        coeffA.push_back(1.2);
        base::DataMatrix* ref_trainData = trainData.get();
        base::DataVector* ref_trainLabels = labels.get();
        sgpp::datadriven::DataBasedRefinementFunctor funcData(ref_grids, ref_alphas,
                                                              ref_trainData,
                                                              ref_trainLabels,
                                                              numPoints,
                                                              levelPenalize,
                                                              coeffA);
        func = &funcZrcr;
        //func = &funcData;
        if (preCompute) {
          // precompute the evals. Needs to be done once per step, before
          // any refinement is done
          func->preComputeEvaluations();
        }
        if (label == -1) {
          func->setGridIndex(0);
        }
        else {
          func->setGridIndex(1);
        }
        grid->getGenerator().refine(*func);

        if (!crossvalidationConfig.silent_) {
          std::cout << "# LearnerSGDE (label " << label << "): ref " << refCnts.at(label)+1 << "/" << numRefSteps
                    << ": " << grids.at(label)->getSize() << std::endl;
        }

        //alpha->resize(grid->getSize());
        //alpha->setAll(0.0); 
        // keep computed alpha values and append zeros only for new points
        alpha->resizeZero(grid->getSize()); 

        //refSteps++;
        refCnts.at(label) += 1;
        doRefine = false;        
      }
      // for plotting only
      //if ((i > 0) && ((i+1) % 50 == 0)) {
      /*if ((i > 0) && ((i+1) % 10 == 0)) {
        //acc = getAccuracy(*trainData, *labels, 0.0);
        acc = getAccuracy(ptestData, ptestLabels, 0.0);
        error.append(1.0 - acc); 
      }*/
    }
    //refSteps = 0;
    numIterations++;

  }
  //write error evaluation to .csv
  /*std::ofstream output;
  output.open("ripley_err_rate_surplusRef_train_1.csv");
  //output.open("ripley_err_rate_impurity_train_1.csv");
  //output.open("banana_err_rate_surplusRef_train_10.csv");
  //output.open("banana_err_rate_impurity_train_1.csv");
  if (output.fail()) {
    std::cout << "failed to create .csv file!" << std::endl;  
  }
  else {
    for (size_t i = 0; i < error.getSize(); i++) {					
      output << error.get(i) << ";" << std::endl;
    }
    output.close();
  }*/
}

void LearnerSGDE::storeResults(base::DataMatrix& testDataset,
                               const base::DataVector& classesReference,
                               const double threshold) {

  base::DataVector classesComputed(testDataset.getNrows());
  predict(testDataset, classesComputed);

  std::ofstream output;
  //write computed classes to .csv
  //output.open("banana_impurity_predicted_train_1.csv");
  //output.open("banana_surplus_predicted_train_1.csv");
  //output.open("ripley_impurity_predicted_train_.csv");
  output.open("ripley_zero_predicted.csv");
  if (output.fail()) {
    std::cout << "failed to create .csv file!" << std::endl;  
  }
  else {
    for (size_t i = 0; i < classesComputed.getSize(); i++) {
      base::DataVector x(2);					
      testDataset.getRow((size_t)i, x);
      output << x[0] << ";" << x[1] << ";" << classesComputed[i] << std::endl;
    }
    output.close();
  }
  //write grid to .csv
  //output.open("banana_impurity_grid_1_train_1.csv");
  //output.open("banana_surplus_grid_1_train_1.csv");
  //output.open("ripley_impurity_grid_1_train_.csv");
  output.open("ripley_zero_grid_1.csv");
  if (output.fail()) {
    std::cout << "failed to create .csv file!" << std::endl;  
  }
  else {
    grid = grids.at(-1);
    base::GridStorage& storage = grid->getStorage();
    base::GridStorage::grid_map_iterator end_iter = storage.end();
    for (base::GridStorage::grid_map_iterator iter = storage.begin(); iter != end_iter; iter++) { 
      base::DataVector gpCoord(testDataset.getNcols());
      storage.getCoordinates(*(iter->first), gpCoord);
      for (size_t d = 0; d < gpCoord.getSize(); d++) {
        if (d < gpCoord.getSize()-1) {
          output << gpCoord[d] << ";";
        }
        else {
          output << gpCoord[d] << std::endl;
        }
      }
    }
    output.close();
  }
  //output.open("banana_impurity_grid_1_train_2.csv");
  //output.open("banana_surplus_grid_2_train_1.csv");
  //output.open("ripley_impurity_grid_2_train_.csv");
  output.open("ripley_zero_grid_2.csv");
  if (output.fail()) {
    std::cout << "failed to create .csv file!" << std::endl;  
  }
  else {
    grid = grids.at(1);
    base::GridStorage& storage = grid->getStorage();
    base::GridStorage::grid_map_iterator end_iter = storage.end();
    for (base::GridStorage::grid_map_iterator iter = storage.begin(); iter != end_iter; iter++) { 
      base::DataVector gpCoord(testDataset.getNcols());
      storage.getCoordinates(*(iter->first), gpCoord);
      for (size_t d = 0; d < gpCoord.getSize(); d++) {
        if (d < gpCoord.getSize()-1) {
          output << gpCoord[d] << ";";
        }
        else {
          output << gpCoord[d] << std::endl;
        }
      }
    }
    output.close();
  }

  //write density function evaluations to .csv
  /*double stepSize = 0.01;
  base::DataMatrix values(0,2);
  //std::cout << values.getNrows() << std::endl;
  base::DataVector range(101);
  for (size_t i = 0; i < 101; i++) {
    range.set(i, stepSize*(static_cast<double>(i)));
  }
  //std::cout << range.getSize() << std::endl;
  for (size_t i = 0; i < range.getSize(); i++) {
    for (size_t j = 0; j < range.getSize(); j++) {
      base::DataVector row(2);
      row.set(1, range.get(i));
      row.set(0, range.get(j));
      values.appendRow(row);
    }
  }
  //std::cout << values.getNrows() << std::endl;
  // evaluate each density function at all points from values
  // and write result to csv file
  for (auto const& g : grids) {
    output.open("density_fun_"+std::to_string(g.first)+"_evals.csv");
    std::unique_ptr<base::OperationEval> opEval(op_factory::createOperationEval(*g.second));
    for (size_t i = 0; i < values.getNrows(); i++) {
      // Get next test sample x 
      base::DataVector x(2);					
      values.getRow(i, x);
      //std::cout << x[0] << " , " << x[1] << std::endl;
      double res = opEval->eval(*alphas.at(g.first), x);
      output << res << ";" << std::endl;
    }
    output.close();
  }*/

}

double LearnerSGDE::getAccuracy(base::DataMatrix& testDataset,
                                const base::DataVector& classesReference,
                                const double threshold) {
  // evaluate test dataset
  base::DataVector classesComputed(testDataset.getNrows());
  predict(testDataset, classesComputed);

  return getAccuracy(classesComputed, classesReference, threshold);
}

double LearnerSGDE::getAccuracy(const base::DataVector& classesComputed,
                                const base::DataVector& classesReference,
                                const double threshold) {
  double result = -1.0;

  if (classesComputed.getSize() != classesReference.getSize()) {
    throw base::application_exception(
        "LearnerSGDE::getAccuracy: lengths of classes vectors do not match!");
  }

  size_t correct = 0;

  for (size_t i = 0; i < classesComputed.getSize(); i++) {
    //std::cout << "computed: " << classesComputed.get(i) << std::endl;
    //std::cout << "reference: " << classesReference.get(i) << std::endl;
    if ((classesComputed.get(i) >= threshold && classesReference.get(i) >= 0.0) ||
       (classesComputed.get(i) < threshold && classesReference.get(i) < 0.0)) {
      correct++;
    }
  }

  result = static_cast<double>(correct) / static_cast<double>(classesComputed.getSize());

  return result;
}

void LearnerSGDE::predict(base::DataMatrix& testData,
                          base::DataVector& computedLabels,
                          bool usePrior) {
  computedLabels.resize(testData.getNrows());
  size_t dim = testData.getNcols(); 

  int compLabel;
  double prior; 

  for (size_t i = 0; i < testData.getNrows(); i++) {
    // Get next test sample x 
    base::DataVector x(dim);					
    testData.getRow((size_t)i, x);
    // predict label using Bayes Theorem
    double max = 0.0;
    // compute each density functions for current test sample x
    for (auto const& g : grids) {
      std::unique_ptr<base::OperationEval> opEval(op_factory::createOperationEval(*g.second));
      double res = opEval->eval(*alphas.at(g.first), x);
      // determine prior
      if (usePrior) {
        prior = static_cast<double>(appearances.at(g.first)) / static_cast<double>(appearances.size());
      }
      else {
        prior = 1.0;
      }
      // multiply by prior
      res *= prior;
      // check if this class probability is larger than current maximum 
      if (res > max) {
        max = res;
        compLabel = g.first;
      }
    } 

    computedLabels.set(i, compLabel); 
  }
}

double LearnerSGDE::computeResidual(base::Grid& grid, base::DataVector& alpha,
                                    base::DataMatrix& test, double lambdaReg) {
  std::unique_ptr<base::OperationMatrix> C = computeRegularizationMatrix(grid);

  base::DataVector rhs(grid.getSize());
  base::DataVector res(grid.getSize());
  datadriven::DensitySystemMatrix SMatrix(grid, test, *C, lambdaReg);
  SMatrix.generateb(rhs);

  SMatrix.mult(alpha, res);

  for (size_t i = 0; i < res.getSize(); i++) {
    res[i] = res[i] - rhs[i];
  }
  return res.l2Norm();
}

std::unique_ptr<base::OperationMatrix> LearnerSGDE::computeRegularizationMatrix(base::Grid& grid) {
  std::unique_ptr<base::OperationMatrix> C;

  if (regularizationConfig.regType_ == datadriven::RegularizationType::Identity) {
    C = op_factory::createOperationIdentity(grid);
  } else if (regularizationConfig.regType_ == datadriven::RegularizationType::Laplace) {
    C = op_factory::createOperationLaplace(grid);
  } else {
    throw base::application_exception("LearnerSGDE::train : unknown regularization type");
  }

  return C;
}

void LearnerSGDE::splitset(std::vector<std::shared_ptr<base::DataMatrix> >& strain,
                           std::vector<std::shared_ptr<base::DataMatrix> >& stest) {
  std::shared_ptr<base::DataMatrix> mydata = std::make_shared<base::DataMatrix>(*trainData);
  base::DataVector p(trainData->getNcols());
  base::DataVector tmp(trainData->getNcols());

  size_t kfold = crossvalidationConfig.kfold_;

  std::vector<size_t> s(kfold);        // size of partition
  std::vector<size_t> ind(kfold + 1);  // index of partition
  size_t n = mydata->getNrows();       // size of data

  if (crossvalidationConfig.shuffle_) {
    if (crossvalidationConfig.seed_ == -1)
      srand(static_cast<unsigned int>(time(0)));
    else
      srand(crossvalidationConfig.seed_);

    for (size_t i = 0; i < mydata->getNrows(); i++) {
      size_t r = i + (static_cast<size_t>(rand()) % (mydata->getNrows() - i));
      mydata->getRow(i, p);
      mydata->getRow(r, tmp);
      mydata->setRow(r, p);
      mydata->setRow(i, tmp);
    }
  }

  // set size of partitions
  if (!crossvalidationConfig.silent_) std::cout << "# kfold: ";

  ind[0] = 0;

  for (size_t i = 0; i < kfold - 1; i++) {
    s[i] = n / kfold;
    ind[i + 1] = ind[i] + s[i];

    if (!crossvalidationConfig.silent_) std::cout << s[i] << " ";
  }

  ind[kfold] = n;
  s[kfold - 1] = n - (kfold - 1) * (n / kfold);

  if (!crossvalidationConfig.silent_) std::cout << s[kfold - 1] << std::endl;

  if (!crossvalidationConfig.silent_) {
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
