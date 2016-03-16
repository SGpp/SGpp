// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <sgpp/datadriven/application/LearnerSGDE.hpp>

#include <sgpp/base/exception/application_exception.hpp>
#include <sgpp/base/grid/generation/functors/SurplusRefinementFunctor.hpp>
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
LearnerSGDEConfiguration::LearnerSGDEConfiguration() : json::JSON() {
  // set default config
  gridConfig.dim_ = 0;
  gridConfig.level_ = 2;
  gridConfig.type_ = base::GridType::Linear;
  gridConfig.maxDegree_ = 1;
  gridConfig.boundaryLevel_ = 0;

  // configure adaptive refinement
  adaptivityConfig.numRefinements_ = 0;
  adaptivityConfig.noPoints_ = 0;

  // configure solver
  solverConfig.type_ = solver::SLESolverType::CG;
  solverConfig.maxIterations_ = 100;
  solverConfig.eps_ = 1e-10;
  solverConfig.threshold_ = 1e-10;

  // configure regularization
  regularizationConfig.regType_ = datadriven::RegularizationType::Laplace;

  // configure learner
  crossvalidationConfig.enable_ = false;
  crossvalidationConfig.kfold_ = 10;
  crossvalidationConfig.lambda_ = 1e-5;
  crossvalidationConfig.lambdaStart_ = 1e-1;
  crossvalidationConfig.lambdaEnd_ = 1e-10;
  crossvalidationConfig.lambdaSteps_ = 10;
  crossvalidationConfig.logScale_ = true;
  crossvalidationConfig.shuffle_ = true;
  crossvalidationConfig.seed_ = 1234567;
  crossvalidationConfig.silent_ = true;
}

LearnerSGDEConfiguration::LearnerSGDEConfiguration(const std::string& fileName)
    : json::JSON(fileName) {
  // initialize structs from file
  // configure grid
  try {
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
    if (this->contains("crossValidation_verbose"))
      crossvalidationConfig.silent_ = (*this)["crossValidation_verbose"].getBool();
  } catch (json::json_exception& e) {
    std::cout << e.what() << std::endl;
  }
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
      samples(nullptr),
      gridConfig(gridConfig),
      adaptivityConfig(adaptivityConfig),
      solverConfig(solverConfig),
      regularizationConfig(regularizationConfig),
      crossvalidationConfig(crossvalidationConfig) {}

LearnerSGDE::LearnerSGDE(LearnerSGDEConfiguration& learnerSGDEConfig)
    : LearnerSGDE(learnerSGDEConfig.gridConfig, learnerSGDEConfig.adaptivityConfig,
                  learnerSGDEConfig.solverConfig, learnerSGDEConfig.regularizationConfig,
                  learnerSGDEConfig.crossvalidationConfig) {}

LearnerSGDE::~LearnerSGDE() {}

void LearnerSGDE::initialize(base::DataMatrix& psamples) {
  base::DataMatrix mysamples(psamples);
  samples = std::make_shared<base::DataMatrix>(mysamples);
  size_t ndim = psamples.getNcols();
  grid = createRegularGrid(ndim);
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
}

// ---------------------------------------------------------------------------

double LearnerSGDE::pdf(base::DataVector& x) {
  return op_factory::createOperationEval(*grid)->eval(*alpha, x);
}

void LearnerSGDE::pdf(base::DataMatrix& points, base::DataVector& res) {
  op_factory::createOperationMultipleEval(*grid, points)->eval(*alpha, res);
}

double LearnerSGDE::mean() {
  return op_factory::createOperationFirstMoment(*grid)->doQuadrature(*alpha);
}

double LearnerSGDE::variance() {
  double secondMoment = op_factory::createOperationSecondMoment(*grid)->doQuadrature(*alpha);

  // use Steiners translation theorem to compute the variance
  double firstMoment = mean();
  double res = secondMoment - firstMoment * firstMoment;
  return res;
}

void LearnerSGDE::cov(base::DataMatrix& cov) { return; }

std::shared_ptr<base::DataVector> LearnerSGDE::getSamples(size_t dim) {
  std::shared_ptr<base::DataVector> isamples = std::make_shared<base::DataVector>(getNsamples());
  samples->getColumn(dim, *isamples);
  return isamples;
}

std::shared_ptr<base::DataMatrix> LearnerSGDE::getSamples() { return samples; }

size_t LearnerSGDE::getDim() { return samples->getNcols(); }

size_t LearnerSGDE::getNsamples() { return samples->getNrows(); }

std::shared_ptr<base::DataVector> LearnerSGDE::getSurpluses() { return alpha; }

std::shared_ptr<base::Grid> LearnerSGDE::getGrid() { return grid; }

// ---------------------------------------------------------------------------

std::shared_ptr<base::Grid> LearnerSGDE::createRegularGrid(size_t ndim) {
  // load grid
  std::unique_ptr<base::Grid> uGrid;
  if (gridConfig.type_ == base::GridType::Linear) {
    uGrid = base::Grid::createLinearGrid(ndim);
  } else if (gridConfig.type_ == base::GridType::LinearL0Boundary) {
    uGrid = base::Grid::createLinearBoundaryGrid(ndim, 0);
  } else if (gridConfig.type_ == base::GridType::LinearBoundary) {
    uGrid = base::Grid::createLinearBoundaryGrid(ndim, 1);
  } else {
    throw base::application_exception("LeanerSGDE::initialize : grid type is not supported");
  }

  uGrid->getGenerator().regular(gridConfig.level_);

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
      grid = createRegularGrid(getDim());
      alpha.setAll(0.0);

      // compute density
      train(*grid, alpha, *(kfold_train[j]), curLambda);
      // get L2 norm of residual for test set
      curMean = computeResidual(*grid, alpha, *(kfold_test[j]), 0.0);
      curMeanAcc += curMean;

      if (!crossvalidationConfig.silent_) {
        std::cout << "# " << curLambda << " " << i << " " << j << " " << curMeanAcc << " "
                  << curMean << std::endl;
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

void LearnerSGDE::train(base::Grid& grid, base::DataVector& alpha, base::DataMatrix& train,
                        double lambdaReg) {
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
    std::unique_ptr<base::OperationMatrix> C = computeRegularizationMatrix(grid);

    datadriven::DensitySystemMatrix SMatrix(grid, train, *C, lambdaReg);
    SMatrix.generateb(rhs);

    if (!crossvalidationConfig.silent_) {
      std::cout << "# LearnerSGDE: Solving " << std::endl;
    }

    solver::ConjugateGradients myCG(solverConfig.maxIterations_, solverConfig.eps_);
    myCG.solve(SMatrix, alpha, rhs, false, false, solverConfig.threshold_);

    if (ref < adaptivityConfig.numRefinements_) {
      if (!crossvalidationConfig.silent_) {
        std::cout << "# LearnerSGDE: Refine grid ... ";
      }

      // Weight surplus with function evaluation at grid points
      std::unique_ptr<base::OperationEval> opEval(op_factory::createOperationEval(grid));
      base::GridIndex* gp;
      base::DataVector p(dim);
      base::DataVector alphaWeight(alpha.getSize());

      for (size_t i = 0; i < grid.getSize(); i++) {
        gp = gridStorage.get(i);
        gp->getCoords(p);
        alphaWeight[i] = alpha[i] * opEval->eval(alpha, p);
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
  std::shared_ptr<base::DataMatrix> mydata = std::make_shared<base::DataMatrix>(*samples);
  base::DataVector p(samples->getNcols());
  base::DataVector tmp(samples->getNcols());

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
