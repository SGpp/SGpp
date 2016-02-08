// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org


#include "ModelFitterDensityEstimation.hpp"

using namespace SGPP::base;
using namespace std;

namespace SGPP {
namespace datadriven {

ModelFitterDensityEstimation::ModelFitterDensityEstimation() {
}

ModelFitterDensityEstimation::~ModelFitterDensityEstimation() {
  // TODO Auto-generated destructor stub
}

void ModelFitterDensityEstimation::fit() {
  size_t dim = train.getNcols();

  GridStorage* gridStorage = grid->getStorage();
  GridGenerator* gridGen = grid->createGridGenerator();
  DataVector rhs(grid->getStorage()->size());
  alpha->resize(grid->getStorage()->size());
  alpha->setAll(0.0);

  if (!learnerSGDEConfig.silent_) {
    cout << "# LearnerSGDE: grid points " << grid.getSize() << endl;
  }

  for (size_t ref = 0; ref <= adaptivityConfig.numRefinements_; ref++) {
    OperationMatrix* C = computeRegularizationMatrix(grid);

    SGPP::datadriven::DensitySystemMatrix SMatrix(grid, train, *C, lambdaReg);
    SMatrix.generateb(rhs);

    if (!learnerSGDEConfig.silent_) {
      cout << "# LearnerSGDE: Solving " << endl;
    }

    SGPP::solver::ConjugateGradients myCG(solverConfig.maxIterations_,
                                          solverConfig.eps_);
    myCG.solve(SMatrix, alpha, rhs, false, false, solverConfig.threshold_);

    if (ref < adaptivityConfig.numRefinements_) {
      if (!learnerSGDEConfig.silent_) {
        cout << "# LearnerSGDE: Refine grid ... ";
      }

      //Weight surplus with function evaluation at grid points
      OperationEval* opEval = SGPP::op_factory::createOperationEval(grid);
      GridIndex* gp;
      DataVector p(dim);
      DataVector alphaWeight(alpha->getSize());

      for (size_t i = 0; i < gridStorage->size(); i++) {
        gp = gridStorage->get(i);
        gp->getCoords(p);
        alphaWeight[i] = alpha[i] * opEval->eval(alpha, p);
      }

      delete opEval;
      opEval = NULL;

      SurplusRefinementFunctor srf(&alphaWeight,
                                   adaptivityConfig.noPoints_, adaptivityConfig.threshold_);
      gridGen->refine(&srf);

      if (!learnerSGDEConfig.silent_) {
        cout << "# LearnerSGDE: ref " << ref << "/"
             << adaptivityConfig.numRefinements_ - 1 << ": "
             << grid->getStorage()->size() << endl;
      }

      alpha->resize(grid->getStorage()->size());
      rhs.resize(grid->getStorage()->size());
      alpha->setAll(0.0);
      rhs.setAll(0.0);
    }

    delete C;
  }

  return;
}

} /* namespace datadriven */
} /* namespace SGPP */
