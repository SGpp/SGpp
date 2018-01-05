// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <sgpp/base/operation/BaseOpFactory.hpp>
#include <sgpp/datadriven/DatadrivenOpFactory.hpp>
#include <sgpp/datadriven/algorithm/SystemMatrixLeastSquaresIdentity.hpp>
#include <sgpp/datadriven/datamining/modules/fitting/ModelFittingLeastSquaresFista.hpp>
#include <sgpp/solver/SLESolver.hpp>

#include <sgpp/base/exception/application_exception.hpp>
#include <sgpp/base/grid/generation/functors/SurplusRefinementFunctor.hpp>
#include <sgpp/base/operation/hash/OperationMultipleEval.hpp>

#include <sgpp/solver/sle/fista/Fista.hpp>
#include <sgpp/solver/sle/fista/RidgeFunction.hpp>
#include <sgpp/solver/sle/fista/LassoFunction.hpp>


#include <sgpp/datadriven/datamining/modules/hpo/OperationMultipleEvalMatrix.hpp>

// TODO(lettrich): allow different refinement types
// TODO(lettrich): allow different refinement criteria

using sgpp::base::Grid;
using sgpp::base::DataMatrix;
using sgpp::base::DataVector;
using sgpp::base::SurplusRefinementFunctor;

using sgpp::base::application_exception;

using sgpp::solver::SLESolver;

namespace sgpp {
namespace datadriven {

ModelFittingLeastSquaresFista::ModelFittingLeastSquaresFista(const FitterConfigurationLeastSquares& config)
    : ModelFittingBase{}, refinementsPerformed{0} {
  this->config = std::unique_ptr<FitterConfiguration>(
      std::make_unique<FitterConfigurationLeastSquares>(config));
  solver = std::unique_ptr<SLESolver>{buildSolver(this->config->getSolverFinalConfig())};
}

// TODO(lettrich): exceptions have to be thrown if not valid.
double ModelFittingLeastSquaresFista::evaluate(const DataVector& sample) const {
  auto opEval = std::unique_ptr<base::OperationEval>{op_factory::createOperationEval(*grid)};
  return opEval->eval(alpha, sample);
}

// TODO(lettrich): exceptions have to be thrown if not valid.
void ModelFittingLeastSquaresFista::evaluate(DataMatrix& samples, DataVector& results) {
  //auto opMultEval = std::unique_ptr<base::OperationMultipleEval>{
  //    op_factory::createOperationMultipleEval(*grid, samples, config->getMultipleEvalConfig())};
  DataMatrix paritymatrix{samples};
  
  
  int ncols = paritymatrix.getNcols();
  paritymatrix.resizeRowsCols(paritymatrix.getNrows(), (ncols*ncols+5)*ncols/6 +1); //EDIT: bias
  //alpha = DataVector{paritymatrix.getNcols()};
  //int parities[3][300];
    for(int p = 0; p < samples.getNrows(); p++){
        for(int k = 0; k < samples.getNcols(); k++){
          paritymatrix.set(p,k,samples.get(p,k));
        }
    }
  
  for(int p = 0; p < paritymatrix.getNrows(); p++){
    int cnt = ncols;
    int cnt2 =(ncols +1)*ncols/2;
    for(int i = 0; i < ncols-1; i++){
      for(int k = i+1; k < ncols; k++){
        double first = paritymatrix.get(p, i) * paritymatrix.get(p, k);
        paritymatrix.set(p, cnt, first);
      //  parities[0][cnt]=i;
    //    parities[1][cnt]=k;
        cnt++;
        for(int m = k+1; m < ncols; m++){
          paritymatrix.set(p, cnt2, first*paritymatrix.get(p, m));
      //    parities[0][cnt2]=i;
        //  parities[1][cnt2]=k;
       //   parities[2][cnt2]=m;
          cnt2++;
        }
      }
    }
    paritymatrix.set(p, paritymatrix.getNcols()-1, 1); //EDIT: bias
  }
  std::cout<< "size: "<<paritymatrix.getNrows()<<" x "<<paritymatrix.getNcols()<<", "<<paritymatrix.get(100, 100)<<", "<<paritymatrix.get(100, 10)<<", "<<samples.get(100, 10)<<std::endl;
  OperationMultipleEvalMatrix opMultEval{*grid, paritymatrix};
  //EDIT: no evaluation
  opMultEval.eval(alpha, results);
}

void ModelFittingLeastSquaresFista::fit(Dataset& newDataset) {
  // clear model
  resetState();
  grid.reset();
  dataset = &newDataset;

  // build grid
  auto& gridConfig = config->getGridConfig();
  gridConfig.dim_ = dataset->getDimension();
  grid = std::unique_ptr<Grid>{buildGrid(config->getGridConfig())};
  // build surplus vector
  alpha = DataVector{dataset->getDimension()};
  //alpha = DataVector{grid->getSize()};

  assembleSystemAndSolve(config->getSolverFinalConfig(), alpha);
}

bool ModelFittingLeastSquaresFista::refine() {
  if (grid != nullptr) {
    if (refinementsPerformed < config->getRefinementConfig().numRefinements_) {
      // create refinement functor
      SurplusRefinementFunctor refinementFunctor(alpha, config->getRefinementConfig().noPoints_,
                                                 config->getRefinementConfig().threshold_);
      // refine grid
      auto noPoints = grid->getSize();
      grid->getGenerator().refine(refinementFunctor);
      if (grid->getSize() > noPoints) {
        // Tell the SLE manager that the grid changed (for interal data structures)
        alpha.resizeZero(grid->getSize());

        assembleSystemAndSolve(config->getSolverRefineConfig(), alpha);
        refinementsPerformed++;
        return true;
      } else {
        return false;
      }
    } else {
      return false;
    }

  } else {
    throw application_exception(
        "ModelFittingLeastSquaresFista: Can't refine before initial grid is created");
    return false;
  }
}

void ModelFittingLeastSquaresFista::update(Dataset& newDataset) {
  if (grid != nullptr) {
    resetState();
    // reassign dataset
    dataset = &newDataset;
    // create sytem matrix
    assembleSystemAndSolve(config->getSolverFinalConfig(), alpha);
  } else {
    fit(newDataset);
  }
}

DMSystemMatrixBase* ModelFittingLeastSquaresFista::buildSystemMatrix(
    Grid& grid, DataMatrix& trainDataset, double lambda,
    OperationMultipleEvalConfiguration& mutipleEvalconfig) const {
  auto systemMatrix = new SystemMatrixLeastSquaresIdentity(grid, trainDataset, lambda);
  systemMatrix->setImplementation(mutipleEvalconfig);

  return systemMatrix;
}

void ModelFittingLeastSquaresFista::resetState() { refinementsPerformed = 0; }

void ModelFittingLeastSquaresFista::assembleSystemAndSolve(const SLESolverConfiguration& solverConfig,
                                                      DataVector& alpha) const {
  auto systemMatrix = std::unique_ptr<DMSystemMatrixBase>(
      buildSystemMatrix(*grid, dataset->getData(), config->getRegularizationConfig().lambda_,
                        config->getMultipleEvalConfig()));

  DataVector b{grid->getSize()};
  systemMatrix->generateb(dataset->getTargets(), b);

  reconfigureSolver(*solver, solverConfig);
  //solver->solve(*systemMatrix, alpha, b, true, verboseSolver, DEFAULT_RES_THRESHOLD);
  
  //sgpp::solver::RidgeFunction g{config->getRegularizationConfig().lambda_};
  solver::LassoFunction g{config->getRegularizationConfig().lambda_};
  solver::Fista<solver::LassoFunction> fista{g};
  
  DataMatrix paritymatrix{dataset->getData()};
  int ncols = paritymatrix.getNcols();
  paritymatrix.resizeRowsCols(paritymatrix.getNrows(), (ncols*ncols+5)*ncols/6 +1); //EDIT: bias
  
  for(int p = 0; p < dataset->getData().getNrows(); p++){
    for(int k = 0; k < dataset->getData().getNcols(); k++){
      paritymatrix.set(p,k,dataset->getData().get(p,k));
    }
  }
  
  alpha = DataVector{paritymatrix.getNcols()};
  int parities[3][300];
  
  for(int p = 0; p < paritymatrix.getNrows(); p++){
    int cnt = ncols;
    int cnt2 =(ncols +1)*ncols/2; //edit
    for(int i = 0; i < ncols-1; i++){
      for(int k = i+1; k < ncols; k++){
        double first = paritymatrix.get(p, i) * paritymatrix.get(p, k);
        paritymatrix.set(p, cnt, first);
        parities[0][cnt]=i;
        parities[1][cnt]=k;
        cnt++;
        for(int m = k+1; m < ncols; m++){
          paritymatrix.set(p, cnt2, first*paritymatrix.get(p, m));
          parities[0][cnt2]=i;
          parities[1][cnt2]=k;
          parities[2][cnt2]=m;
          cnt2++;
        }
      }
    }
    paritymatrix.set(p, paritymatrix.getNcols()-1, 1); //EDIT: bias
  }
  //auto opMultEval = std::unique_ptr<base::OperationMultipleEval>{
    //  op_factory::createOperationMultipleEval(*grid, dataset->getData(), config->getMultipleEvalConfig())};
  OperationMultipleEvalMatrix opMultEval{*grid, paritymatrix};
  fista.solve(opMultEval, alpha, dataset->getTargets(), 100, DEFAULT_RES_THRESHOLD);
  
  for(int i=0;i<paritymatrix.getNcols();i++){
      std::cout<< "Alpha "<<i<<": "<<alpha.get(i) << ": "<<parities[0][i]<<";"<<parities[1][i]<<";"<<parities[2][i]<<std::endl;
  }
}

}  // namespace datadriven
}  // namespace sgpp
