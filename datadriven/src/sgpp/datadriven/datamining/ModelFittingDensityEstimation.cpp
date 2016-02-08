// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org


#include "ModelFittingDensityEstimation.hpp"

using namespace SGPP::base;


namespace SGPP {
  namespace datadriven {

    ModelFittingDensityEstimation::ModelFittingDensityEstimation(SGPP::datadriven::DataMiningConfiguration config) : SGPP::datadriven::ModelFittingBase(config) {
      config = *static_cast<SGPP::datadriven::DataMiningConfigurationDensityEstimation*>(config.clone());
    }

    ModelFittingDensityEstimation::~ModelFittingDensityEstimation() {
      // TODO Auto-generated destructor stub
    }

    void ModelFittingDensityEstimation::fit() {
      // TODO: set the dimensionality
      size_t dim = 1;

      GridStorage* gridStorage = grid->getStorage();
      GridGenerator* gridGen = grid->createGridGenerator();
      DataVector rhs(grid->getStorage()->size());
      alpha->resize(grid->getStorage()->size());
      alpha->setAll(0.0);

      if (!config.sgdeConfig.silent_) {
        std::cout << "# LearnerSGDE: grid points " << grid->getSize() << std::endl;
      }

      //  for (size_t ref = 0; ref <= config.adaptivityConfig.numRefinements_; ref++) {
      //    OperationMatrix* C = computeRegularizationMatrix(grid);

      //    SGPP::datadriven::DensitySystemMatrix SMatrix(grid, train, *C, lambdaReg);
      //    SMatrix.generateb(rhs);
      //
      //    if (!config.sgdeConfig.silent_) {
      //      cout << "# LearnerSGDE: Solving " << endl;
      //    }
      //
      //    SGPP::solver::ConjugateGradients myCG(config.solverConfig.maxIterations_,
      //                                          config.solverConfig.eps_);
      //    myCG.solve(SMatrix, alpha, rhs, false, false, config.solverConfig.threshold_);
      //
      //    if (ref < config.adaptivityConfig.numRefinements_) {
      //      if (!config.sgdeConfig.silent_) {
      //        cout << "# LearnerSGDE: Refine grid ... ";
      //      }
      //
      //      //Weight surplus with function evaluation at grid points
      //      OperationEval* opEval = SGPP::op_factory::createOperationEval(grid);
      //      GridIndex* gp;
      //      DataVector p(dim);
      //      DataVector alphaWeight(alpha->getSize());
      //
      //      for (size_t i = 0; i < gridStorage->size(); i++) {
      //        gp = gridStorage->get(i);
      //        gp->getCoords(p);
      //        alphaWeight[i] = alpha[i] * opEval->eval(alpha, p);
      //      }
      //
      //      delete opEval;
      //      opEval = NULL;
      //
      //      SurplusRefinementFunctor srf(&alphaWeight,
      //                                   config.adaptivityConfig.noPoints_,
      //                                   config.adaptivityConfig.threshold_);
      //      gridGen->refine(&srf);
      //
      //      if (!config.sgdeConfig.silent_) {
      //        cout << "# LearnerSGDE: ref " << ref << "/"
      //             << config.adaptivityConfig.numRefinements_ - 1 << ": "
      //             << grid->getStorage()->size() << endl;
      //      }
      //
      //      alpha->resize(grid->getStorage()->size());
      //      rhs.resize(grid->getStorage()->size());
      //      alpha->setAll(0.0);
      //      rhs.setAll(0.0);
      //    }
      //
      //    delete C;
      //  }

      return;
    }

  } /* namespace datadriven */
} /* namespace SGPP */
