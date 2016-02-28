// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <sgpp/datadriven/datamining/ModelFittingDensityEstimation.hpp>
#include <sgpp/base/operation/BaseOpFactory.hpp>

#include <sgpp/datadriven/algorithm/DensitySystemMatrix.hpp>
#include <sgpp/solver/sle/ConjugateGradients.hpp>
#include <sgpp/base/grid/generation/functors/SurplusRefinementFunctor.hpp>
#include <sgpp/base/exception/application_exception.hpp>
#include <sgpp/base/tools/json/json_exception.hpp>

using namespace SGPP::base;

namespace SGPP {
namespace datadriven {

ModelFittingDensityEstimation::ModelFittingDensityEstimation(
    datadriven::DataMiningConfiguration config)
    : datadriven::ModelFittingBase() {
  // load all the needed variables
  // ...
  // TODO: set default values
  // set default in get like .get("Linear")

  try {
    gridConfig.dim_ = 0;
    gridConfig.level_ = static_cast<int>(config["grid"]["level"].getInt());
    gridConfig.type_ = config.stringToGridType(config["grid"]["type"].get());

    // configure adaptive refinement
    adaptivityConfig.numRefinements_ =
        config["refinement"]["numSteps"].getUInt();
    adaptivityConfig.noPoints_ = config["refinement"]["numPoints"].getUInt();

    // configure solver
    solverConfig.type_ =
        config.stringToSolverType(config["solver"]["type"].get());
    solverConfig.maxIterations_ = config["solver"]["maxIterations"].getUInt();
    solverConfig.eps_ = config["solver"]["eps"].getDouble();
    solverConfig.threshold_ = config["solver"]["threshold"].getDouble();

    // configure regularization
    regularizationConfig.regType_ = config.stringToRegularizationType(
        config["regularization"]["type"].get());

    // configure learner
    sgdeConfig.doCrossValidation_ =
        config["crossValidation"]["doCrossValidation"].getBool();
    sgdeConfig.kfold_ = config["crossValidation"]["kfold"].getUInt();
    sgdeConfig.lambdaStart_ =
        config["crossValidation"]["lambdaStart"].getDouble();
    sgdeConfig.lambdaEnd_ = config["crossValidation"]["lambdaEnd"].getDouble();
    sgdeConfig.lambdaSteps_ =
        config["crossValidation"]["lambdaSteps"].getUInt();
    sgdeConfig.logScale_ = config["crossValidation"]["logScale"].getBool();
    sgdeConfig.shuffle_ = config["crossValidation"]["shuffle"].getBool();
    sgdeConfig.seed_ = config["crossValidation"]["seed"].getUInt();
    sgdeConfig.silent_ = config["crossValidation"]["verbose"].getBool();
  } catch (json::json_exception& e) {
    std::cout << e.what() << std::endl;
  }
}

ModelFittingDensityEstimation::~ModelFittingDensityEstimation() {}

void ModelFittingDensityEstimation::fit(datadriven::Dataset& dataset) {
  DataMatrix samples = dataset.getData();
  size_t numDims = samples.getNcols();

  initializeGrid(gridConfig);

  GridStorage* gridStorage = grid->getStorage();
  GridGenerator* gridGen = grid->createGridGenerator();
  DataVector rhs(grid->getStorage()->size());
  alpha->resize(grid->getStorage()->size());
  alpha->setAll(0.0);

  if (!sgdeConfig.silent_) {
    std::cout << "# LearnerSGDE: grid points " << grid->getSize() << std::endl;
  }

  for (size_t ref = 0; ref <= adaptivityConfig.numRefinements_; ref++) {
    std::shared_ptr<OperationMatrix> C =
        getRegularizationMatrix(regularizationConfig.regType_);
    datadriven::DensitySystemMatrix SMatrix(*grid, samples, *C, 1e-10);
    SMatrix.generateb(rhs);

    if (!sgdeConfig.silent_) {
      std::cout << "# LearnerSGDE: Solving " << std::endl;
    }

    solver::ConjugateGradients myCG(solverConfig.maxIterations_,
                                    solverConfig.eps_);
    myCG.solve(SMatrix, *alpha, rhs, false, false, solverConfig.threshold_);

    if (ref < adaptivityConfig.numRefinements_) {
      if (!sgdeConfig.silent_) {
        std::cout << "# LearnerSGDE: Refine grid ... ";
      }

      // Weight surplus with function evaluation at grid points
      OperationEval* opEval = op_factory::createOperationEval(*grid);
      GridIndex* gp;
      DataVector p(numDims);
      DataVector alphaWeight(alpha->getSize());

      for (size_t i = 0; i < gridStorage->size(); i++) {
        gp = gridStorage->get(i);
        gp->getCoords(p);
        alphaWeight[i] = alpha->get(i) * opEval->eval(*alpha, p);
      }

      delete opEval;
      opEval = NULL;

      base::SurplusRefinementFunctor srf(&alphaWeight,
                                         adaptivityConfig.noPoints_,
                                         adaptivityConfig.threshold_);
      gridGen->refine(&srf);

      if (!sgdeConfig.silent_) {
        std::cout << "# LearnerSGDE: ref " << ref << "/"
                  << adaptivityConfig.numRefinements_ - 1 << ": "
                  << grid->getStorage()->size() << std::endl;
      }

      alpha->resize(grid->getStorage()->size());
      rhs.resize(grid->getStorage()->size());
      alpha->setAll(0.0);
      rhs.setAll(0.0);
    }
  }
}

void ModelFittingDensityEstimation::refine() {
  // TODO implement
}

void ModelFittingDensityEstimation::update(datadriven::Dataset& dataset) {
  // TODO implement
}

} /* namespace datadriven */
} /* namespace SGPP */
//git please do not delete me
