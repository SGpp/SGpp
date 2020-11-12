// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <sgpp/base/exception/application_exception.hpp>
#include <sgpp/base/exception/data_exception.hpp>
#include <sgpp/base/grid/Grid.hpp>
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
#include <sgpp/datadriven/application/SparseGridDensityEstimator.hpp>
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
#include <string>
#include <vector>

namespace sgpp {
namespace datadriven {

// --------------------------------------------------------------------------------------------
SparseGridDensityEstimatorConfiguration::SparseGridDensityEstimatorConfiguration() : json::JSON() {
  initConfig();
}

SparseGridDensityEstimatorConfiguration::SparseGridDensityEstimatorConfiguration(
    const std::string& fileName)
    : json::JSON(fileName) {
  // initialize structs with default values
  initConfig();
  // initialize structs from file
  // configure grid
  try {
    if (this->contains("grid_filename")) gridConfig.filename_ = (*this)["grid_filename"].get();
    if (this->contains("grid_dim")) gridConfig.dim_ = (*this)["grid_level"].getUInt();
    if (this->contains("grid_level"))
      gridConfig.level_ = static_cast<int>((*this)["grid_level"].getInt());
    if (this->contains("grid_type"))
      gridConfig.type_ = base::Grid::stringToGridType((*this)["grid_type"].get());
    if (this->contains("grid_maxDegree"))
      gridConfig.maxDegree_ = (*this)["grid_maxDegree"].getUInt();
    if (this->contains("grid_boundaryLevel"))
      gridConfig.boundaryLevel_ =
          static_cast<base::level_t>((*this)["grid_boundaryLevel"].getUInt());

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
    if (this->contains("solver_verbose"))
      solverConfig.verbose_ = (*this)["solver_verbose"].getBool();

    // configure regularization
    if (this->contains("regularization_type"))
      regularizationConfig.type_ = stringToRegularizationType((*this)["regularization_type"].get());

    // configure cross validation
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

    // configure learner
    if (this->contains("sgde_makePositive"))
      sgdeConfig.makePositive_ = (*this)["sgde_makePositive"].getBool();
    if (this->contains("sgde_makePositive_candidateSearchAlgorithm")) {
      std::string candidateSearchStr = (*this)["sgde_makePositive_candidateSearchAlgorithm"].get();
      if ((strcmp(candidateSearchStr.c_str(), "intersections") == 0)) {
        sgdeConfig.makePositive_candidateSearchAlgorithm_ =
            datadriven::MakePositiveCandidateSearchAlgorithm::Intersections;
      } else if ((strcmp(candidateSearchStr.c_str(), "joined") == 0)) {
        sgdeConfig.makePositive_candidateSearchAlgorithm_ =
            datadriven::MakePositiveCandidateSearchAlgorithm::IntersectionsJoin;
      } else if ((strcmp(candidateSearchStr.c_str(), "hybrid") == 0)) {
        sgdeConfig.makePositive_candidateSearchAlgorithm_ =
            datadriven::MakePositiveCandidateSearchAlgorithm::HybridFullIntersections;
      } else if ((strcmp(candidateSearchStr.c_str(), "fullGrid") == 0)) {
        sgdeConfig.makePositive_candidateSearchAlgorithm_ =
            datadriven::MakePositiveCandidateSearchAlgorithm::FullGrid;
      } else {
        throw sgpp::base::application_exception("candidate search algorithm is unknown");
      }
    }
    if (this->contains("sgde_makePositive_interpolationAlgorithm")) {
      std::string candidateSearchStr = (*this)["sgde_makePositive_interpolationAlgorithm"].get();
      if ((strcmp(candidateSearchStr.c_str(), "setToZero") == 0)) {
        sgdeConfig.makePositive_interpolationAlgorithm_ =
            datadriven::MakePositiveInterpolationAlgorithm::SetToZero;
      } else if ((strcmp(candidateSearchStr.c_str(), "interpolateBoundaries1d") == 0)) {
        sgdeConfig.makePositive_interpolationAlgorithm_ =
            datadriven::MakePositiveInterpolationAlgorithm::InterpolateBoundaries1d;
      } else if ((strcmp(candidateSearchStr.c_str(), "interpolateExp") == 0)) {
        sgdeConfig.makePositive_interpolationAlgorithm_ =
            datadriven::MakePositiveInterpolationAlgorithm::InterpolateExp;
      } else {
        throw sgpp::base::application_exception("interpolation algorithm is unknown");
      }
    }
    if (this->contains("sgde_makePositive_generateConsistentGrid"))
      sgdeConfig.makePositive_generateConsistentGrid_ =
          (*this)["sgde_makePositive_generateConsistentGrid"].getBool();
    if (this->contains("sgde_unitIntegrand"))
      sgdeConfig.unitIntegrand_ = (*this)["sgde_unitIntegrand"].getBool();
  } catch (json::json_exception& e) {
    std::cout << e.what() << std::endl;
  }
}

void SparseGridDensityEstimatorConfiguration::initConfig() {
  // set default config
  gridConfig.dim_ = 0;
  gridConfig.level_ = 6;
  gridConfig.type_ = base::GridType::Linear;
  gridConfig.maxDegree_ = 3;
  gridConfig.boundaryLevel_ = 0;

  // configure adaptive refinement
  adaptivityConfig.numRefinements_ = 0;
  adaptivityConfig.numRefinementPoints_ = 5;

  // configure solver
  solverConfig.type_ = solver::SLESolverType::CG;
  solverConfig.maxIterations_ = 1000;
  solverConfig.eps_ = 1e-10;
  solverConfig.threshold_ = 1e-14;
  solverConfig.verbose_ = false;

  // configure regularization
  regularizationConfig.type_ = datadriven::RegularizationType::Laplace;

  // configure cross validation
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

  // configure learner
  sgdeConfig.makePositive_ = false;
  sgdeConfig.makePositive_candidateSearchAlgorithm_ =
      datadriven::MakePositiveCandidateSearchAlgorithm::Intersections;
  sgdeConfig.makePositive_interpolationAlgorithm_ =
      datadriven::MakePositiveInterpolationAlgorithm::SetToZero;
  sgdeConfig.makePositive_generateConsistentGrid_ = true;

  sgdeConfig.unitIntegrand_ = false;
}

SparseGridDensityEstimatorConfiguration* SparseGridDensityEstimatorConfiguration::clone() {
  SparseGridDensityEstimatorConfiguration* clone =
      new SparseGridDensityEstimatorConfiguration(*this);
  return clone;
}

sgpp::datadriven::RegularizationType
SparseGridDensityEstimatorConfiguration::stringToRegularizationType(
    std::string& regularizationType) {
  if (regularizationType.compare("Identity") == 0) {
    return sgpp::datadriven::RegularizationType::Identity;
  } else if (regularizationType.compare("Laplace") == 0) {
    return sgpp::datadriven::RegularizationType::Laplace;
  } else {
    throw sgpp::base::application_exception("regularization type is unknown");
  }
}

sgpp::solver::SLESolverType SparseGridDensityEstimatorConfiguration::stringToSolverType(
    std::string& solverType) {
  if (solverType.compare("CG")) {
    return sgpp::solver::SLESolverType::CG;
  } else if (solverType.compare("BiCGSTAB")) {
    return sgpp::solver::SLESolverType::BiCGSTAB;
  } else {
    throw sgpp::base::application_exception("solver type is unknown");
  }
}

// --------------------------------------------------------------------------------------------
SparseGridDensityEstimator::SparseGridDensityEstimator(
    sgpp::base::RegularGridConfiguration& gridConfig,
    sgpp::base::AdaptivityConfiguration& adaptivityConfig,
    sgpp::solver::SLESolverConfiguration& solverConfig,
    sgpp::datadriven::RegularizationConfiguration& regularizationConfig,
    CrossvalidationForRegularizationConfiguration& crossvalidationConfig,
    SGDEConfiguration& sgdeConfig)
    : grid(nullptr),
      alpha(nullptr),
      samples(nullptr),
      gridConfig(gridConfig),
      adaptivityConfig(adaptivityConfig),
      solverConfig(solverConfig),
      regularizationConfig(regularizationConfig),
      crossvalidationConfig(crossvalidationConfig),
      sgdeConfig(sgdeConfig) {}

SparseGridDensityEstimator::SparseGridDensityEstimator(
    SparseGridDensityEstimatorConfiguration& learnerSGDEConfig)
    : SparseGridDensityEstimator(
          learnerSGDEConfig.gridConfig, learnerSGDEConfig.adaptivityConfig,
          learnerSGDEConfig.solverConfig, learnerSGDEConfig.regularizationConfig,
          learnerSGDEConfig.crossvalidationConfig, learnerSGDEConfig.sgdeConfig) {}

SparseGridDensityEstimator::SparseGridDensityEstimator(
    const SparseGridDensityEstimator& learnerSGDE) {
  grid = learnerSGDE.grid;
  alpha = learnerSGDE.alpha;
  samples = learnerSGDE.samples;
  gridConfig = learnerSGDE.gridConfig;
  adaptivityConfig = learnerSGDE.adaptivityConfig;
  solverConfig = learnerSGDE.solverConfig;
  regularizationConfig = learnerSGDE.regularizationConfig;
  crossvalidationConfig = learnerSGDE.crossvalidationConfig;
  sgdeConfig = learnerSGDE.sgdeConfig;
}

SparseGridDensityEstimator::SparseGridDensityEstimator(base::Grid& grid, base::DataVector& alpha,
                                                       base::DataMatrix& samples) {
  std::shared_ptr<base::Grid> sGrid(grid.clone());
  this->grid = sGrid;
  this->alpha = std::make_shared<base::DataVector>(alpha);
  this->samples = std::make_shared<base::DataMatrix>(samples);
}

SparseGridDensityEstimator::~SparseGridDensityEstimator() {}

// -----------------------------------------------------------------------------------------------

void SparseGridDensityEstimator::initialize(base::DataMatrix& psamples) {
  samples = std::make_shared<base::DataMatrix>(psamples);
  gridConfig.dim_ = psamples.getNcols();
  if (gridConfig.dim_ == 0) {
    throw base::application_exception(
        "LearnerSGDE::initialize - the data is corrupted -> no columns available");
  }
  grid = createRegularGrid();
  alpha = std::make_shared<base::DataVector>(grid->getSize());

  // optimize the regularization parameter
  double lambdaReg = 0.0;

  if (crossvalidationConfig.enable_) {
    lambdaReg = optimizeLambdaCV();
  } else {
    lambdaReg = crossvalidationConfig.lambda_;
  }

  // learn the data -> do the density estimation
  train(*grid, *alpha, *samples, lambdaReg);

  // make the density positive
  if (sgdeConfig.makePositive_) {
    op_factory::createOperationMakePositive(sgdeConfig.makePositive_candidateSearchAlgorithm_,
                                            sgdeConfig.makePositive_interpolationAlgorithm_,
                                            sgdeConfig.makePositive_generateConsistentGrid_)
        ->makePositive(*grid, *alpha, true);
  }

  // force the integral to be 1
  if (sgdeConfig.unitIntegrand_) {
    double vol = op_factory::createOperationQuadrature(*grid)->doQuadrature(*alpha);
    alpha->mult(1. / vol);
  }
}

// ---------------------------------------------------------------------------

double SparseGridDensityEstimator::pdf(base::DataVector& x) {
  std::unique_ptr<base::OperationEval> opEval(op_factory::createOperationEval(*grid));
  return opEval->eval(*alpha, x);
}

void SparseGridDensityEstimator::pdf(base::DataMatrix& points, base::DataVector& res) {
  std::unique_ptr<base::OperationMultipleEval> opEval(
      op_factory::createOperationMultipleEval(*grid, points));
  opEval->eval(*alpha, res);
}

double SparseGridDensityEstimator::mean(base::Grid& grid, base::DataVector& alpha) {
  std::unique_ptr<base::OperationFirstMoment> opFirstMoment(
      op_factory::createOperationFirstMoment(grid));
  return opFirstMoment->doQuadrature(alpha);
}

double SparseGridDensityEstimator::mean() { return mean(*grid, *alpha); }

double SparseGridDensityEstimator::variance(base::Grid& grid, base::DataVector& alpha) {
  std::unique_ptr<base::OperationFirstMoment> opFirstMoment(
      op_factory::createOperationFirstMoment(grid));
  double firstMoment = opFirstMoment->doQuadrature(alpha);
  std::unique_ptr<base::OperationSecondMoment> opSecondMoment(
      op_factory::createOperationSecondMoment(grid));
  double secondMoment = opSecondMoment->doQuadrature(alpha);
  // use Steiners translation theorem to compute the variance
  return secondMoment - firstMoment * firstMoment;
}

double SparseGridDensityEstimator::variance() { return variance(*grid, *alpha); }

void SparseGridDensityEstimator::cov(base::DataMatrix& cov, base::DataMatrix* bounds) {
  std::unique_ptr<datadriven::OperationCovariance> opCov(
      op_factory::createOperationCovariance(*grid));
  opCov->doQuadrature(*alpha, cov, bounds);
}

std::shared_ptr<base::DataVector> SparseGridDensityEstimator::getSamples(size_t dim) {
  std::shared_ptr<base::DataVector> isamples = std::make_shared<base::DataVector>(getNsamples());
  samples->getColumn(dim, *isamples);
  return isamples;
}

std::shared_ptr<base::DataMatrix> SparseGridDensityEstimator::getSamples() { return samples; }

size_t SparseGridDensityEstimator::getDim() { return gridConfig.dim_; }

size_t SparseGridDensityEstimator::getNsamples() { return samples->getNrows(); }

base::DataVector& SparseGridDensityEstimator::getSurpluses() { return *alpha; }

base::Grid& SparseGridDensityEstimator::getGrid() { return *grid; }

// ---------------------------------------------------------------------------

std::shared_ptr<base::Grid> SparseGridDensityEstimator::createRegularGrid() {
  // load grid
  std::shared_ptr<base::Grid> sGrid(base::Grid::createGrid(gridConfig));
  sGrid->getGenerator().regular(gridConfig.level_);
  return sGrid;
}

double SparseGridDensityEstimator::optimizeLambdaCV() {
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
                std::max(static_cast<double>(crossvalidationConfig.lambdaSteps_) / 10.0,
                         static_cast<double>(1.0))) ==
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

void SparseGridDensityEstimator::train(base::Grid& grid, base::DataVector& alpha,
                                       base::DataMatrix& train, double lambdaReg) {
  size_t dim = train.getNcols();

  base::GridStorage& gridStorage = grid.getStorage();
  base::GridGenerator& gridGen = grid.getGenerator();
  base::DataVector rhs(grid.getSize());
  alpha.resize(grid.getSize());
  alpha.setAll(0.0);

  if (!crossvalidationConfig.silent_) {
    std::cout << "# LearnerSGDE: grid points " << grid.getSize() << std::endl;
  }

  for (size_t ref = 0; ref <= adaptivityConfig.numRefinements_; ref++) {
    auto sMatrix = computeDensitySystemMatrix(grid, train, lambdaReg);
    sMatrix->generateb(rhs);

    if (!crossvalidationConfig.silent_) {
      std::cout << "# LearnerSGDE: Solving " << std::endl;
    }

    solver::ConjugateGradients myCG(solverConfig.maxIterations_, solverConfig.eps_);
    myCG.solve(*sMatrix, alpha, rhs, false, solverConfig.verbose_, solverConfig.threshold_);

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
        alphaWeight[i] = alpha[i] * opEval->eval(alpha, p);
      }

      base::SurplusRefinementFunctor srf(alphaWeight, adaptivityConfig.numRefinementPoints_,
                                         adaptivityConfig.refinementThreshold_);
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

double SparseGridDensityEstimator::computeResidual(base::Grid& grid, base::DataVector& alpha,
                                                   base::DataMatrix& test, double lambdaReg) {
  base::DataVector rhs(grid.getSize());
  base::DataVector res(grid.getSize());
  auto sMatrix = computeDensitySystemMatrix(grid, test, lambdaReg);
  sMatrix->generateb(rhs);

  sMatrix->mult(alpha, res);

  for (size_t i = 0; i < res.getSize(); i++) {
    res[i] = res[i] - rhs[i];
  }
  return res.l2Norm();
}

base::OperationMatrix* SparseGridDensityEstimator::computeLTwoDotProductMatrix(base::Grid& grid) {
  return op_factory::createOperationLTwoDotProduct(grid);
}

base::OperationMultipleEval* SparseGridDensityEstimator::computeMultipleEvalMatrix(
    base::Grid& grid, base::DataMatrix& train) {
  if (grid.getType() == base::GridType::Bspline ||
      grid.getType() == base::GridType::BsplineBoundary ||
      grid.getType() == base::GridType::BsplineClenshawCurtis ||
      grid.getType() == base::GridType::ModBsplineClenshawCurtis ||
      grid.getType() == base::GridType::ModBspline ||
      grid.getType() == base::GridType::PolyClenshawCurtis ||
      grid.getType() == base::GridType::PolyClenshawCurtisBoundary ||
      grid.getType() == base::GridType::ModPolyClenshawCurtis) {
    return op_factory::createOperationMultipleEvalNaive(grid, train);
  } else {
    return op_factory::createOperationMultipleEval(grid, train);
  }
}

base::OperationMatrix* SparseGridDensityEstimator::computeRegularizationMatrix(base::Grid& grid) {
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

std::shared_ptr<datadriven::DensitySystemMatrix>
SparseGridDensityEstimator::computeDensitySystemMatrix(base::Grid& grid, base::DataMatrix& train,
                                                       double lambdaReg) {
  auto A = computeLTwoDotProductMatrix(grid);
  auto B = computeMultipleEvalMatrix(grid, train);
  auto C = computeRegularizationMatrix(grid);

  return std::make_shared<datadriven::DensitySystemMatrix>(A, B, C, lambdaReg, train.getNrows());
}

void SparseGridDensityEstimator::splitset(std::vector<std::shared_ptr<base::DataMatrix> >& strain,
                                          std::vector<std::shared_ptr<base::DataMatrix> >& stest) {
  std::shared_ptr<base::DataMatrix> mydata = std::make_shared<base::DataMatrix>(*samples);
  base::DataVector p(samples->getNcols());
  base::DataVector tmp(samples->getNcols());

  size_t kfold = crossvalidationConfig.kfold_;

  std::vector<size_t> s(kfold);        // size of partition
  std::vector<size_t> ind(kfold + 1);  // index of partition
  size_t n = mydata->getNrows();       // size of data

  if (crossvalidationConfig.shuffle_) {
    if (crossvalidationConfig.seed_ == -1)
      srand(static_cast<unsigned int>(time(nullptr)));
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

SparseGridDensityEstimator* SparseGridDensityEstimator::margToDimX(size_t idim) {
  std::unique_ptr<datadriven::OperationDensityMargTo1D> opMarg(
      op_factory::createOperationDensityMargTo1D(*grid));

  base::Grid* margGrid = nullptr;
  base::DataVector* margAlpha = nullptr;
  opMarg->margToDimX(&(*alpha), margGrid, margAlpha, idim);

  // create new Learner
  base::DataMatrix margSamples(samples->getNrows(), 1);
  base::DataVector samples1d(samples->getNrows());
  samples->getColumn(idim, samples1d);
  margSamples.setColumn(0, samples1d);

  SparseGridDensityEstimator* ans =
      new SparseGridDensityEstimator(*margGrid, *margAlpha, margSamples);

  delete margGrid;
  delete margAlpha;

  return ans;
}

SparseGridDensityEstimator* SparseGridDensityEstimator::marginalize(size_t idim) {
  std::unique_ptr<datadriven::OperationDensityMarginalize> opMarg(
      op_factory::createOperationDensityMarginalize(*grid));

  base::Grid* margGrid = nullptr;
  base::DataVector margAlpha;
  opMarg->doMarginalize(*alpha, margGrid, margAlpha, static_cast<unsigned int>(idim));

  // create new Learner
  size_t numDims = samples->getNcols();
  base::DataMatrix margSamples(samples->getNrows(), numDims - 1);
  base::DataVector samples1d(samples->getNrows());
  size_t mdim = 0;
  for (size_t jdim = 0; jdim < numDims; jdim++) {
    if (jdim != idim) {
      samples->getColumn(jdim, samples1d);
      margSamples.setColumn(mdim, samples1d);
      mdim++;
    }
  }

  SparseGridDensityEstimator* ans =
      new SparseGridDensityEstimator(*margGrid, margAlpha, margSamples);

  delete margGrid;

  return ans;
}

}  // namespace datadriven
}  // namespace sgpp
