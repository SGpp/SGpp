/* ****************************************************************************
* Copyright (C) 2009 Technische Universitaet Muenchen                         *
* This file is part of the SG++ project. For conditions of distribution and   *
* use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
**************************************************************************** */
// @author Alexander Heinecke (Alexander.Heinecke@mytum.de)

#include "finance/algorithm/BlackScholesParabolicPDESolverSystemEuroAmer.hpp"
#include "base/exception/algorithm_exception.hpp"
#include "base/grid/generation/functors/SurplusCoarseningFunctor.hpp"
#include "base/grid/generation/functors/SurplusRefinementFunctor.hpp"
#include "base/operation/BaseOpFactory.hpp"
#include "pde/operation/PdeOpFactory.hpp"
#include "finance/operation/FinanceOpFactory.hpp"
#include <cmath>

namespace sg {
  namespace finance {

    BlackScholesParabolicPDESolverSystemEuroAmer::BlackScholesParabolicPDESolverSystemEuroAmer(sg::base::Grid& SparseGrid, sg::base::DataVector& alpha, sg::base::DataVector& mu,
        sg::base::DataVector& sigma, sg::base::DataMatrix& rho, double r, double TimestepSize, std::string OperationMode,
        double dStrike, std::string option_type,
        bool bLogTransform, bool useCoarsen, double coarsenThreshold, std::string adaptSolveMode,
        int numCoarsenPoints, double refineThreshold, std::string refineMode, sg::base::GridIndex::level_type xLevel) {
      this->BoundGrid = &SparseGrid;
      this->alpha_complete = &alpha;

      this->alpha_complete_old = new sg::base::DataVector(*this->alpha_complete);
      this->alpha_complete_tmp = new sg::base::DataVector(*this->alpha_complete);
      this->oldGridStorage = new sg::base::GridStorage(*(this->BoundGrid)->getStorage());
      this->secondGridStorage = new sg::base::GridStorage(*(this->BoundGrid)->getStorage());

      this->InnerGrid = NULL;
      this->alpha_inner = NULL;
      this->tOperationMode = OperationMode;
      this->TimestepSize = TimestepSize;
      this->TimestepSize_old = TimestepSize;
      this->BoundaryUpdate = new sg::base::DirichletUpdateVector(SparseGrid.getStorage());
      this->GridConverter = new sg::base::DirichletGridConverter();
      this->r = r;
      this->mus = &mu;
      this->sigmas = &sigma;
      this->rhos = &rho;
      this->BSalgoDims = this->BoundGrid->getAlgorithmicDimensions();
      this->nExecTimesteps = 0;

      // throw exception if grid dimensions not equal algorithmic dimensions
      if (this->BSalgoDims.size() != this->BoundGrid->getStorage()->dim()) {
        throw sg::base::algorithm_exception("BlackScholesParabolicPDESolverSystemEuropean::BlackScholesParabolicPDESolverSystemEuropean : Number of algorithmic dimensions is not equal to the number of grid's dimensions.");
      }

      // test if number of dimensions in the coefficients match the numbers of grid dimensions (mu and sigma)
      if (this->BoundGrid->getStorage()->dim() != this->mus->getSize() || this->BoundGrid->getStorage()->dim() != this->sigmas->getSize()) {
        throw sg::base::algorithm_exception("BlackScholesParabolicPDESolverSystemEuropean::BlackScholesParabolicPDESolverSystemEuropean : Dimension of mu and sigma parameters don't match the grid's dimensions!");
      }

      // test if number of dimensions in the coefficients match the numbers of grid dimensions (rho)
      if (this->BoundGrid->getStorage()->dim() != this->rhos->getNrows() || this->BoundGrid->getStorage()->dim() != this->rhos->getNcols()) {
        throw sg::base::algorithm_exception("BlackScholesParabolicPDESolverSystemEuropean::BlackScholesParabolicPDESolverSystemEuropean : Row or col of rho parameter don't match the grid's dimensions!");
      }

      // test if all algorithmic dimensions are inside the grid's dimensions
      for (size_t i = 0; i < this->BSalgoDims.size(); i++) {
        if (this->BSalgoDims[i] >= this->BoundGrid->getStorage()->dim()) {
          throw sg::base::algorithm_exception("BlackScholesParabolicPDESolverSystemEuropean::BlackScholesParabolicPDESolverSystemEuropean : Minimum one algorithmic dimension is not inside the grid's dimensions!");
        }
      }

      // test if there are double algorithmic dimensions
      std::vector<size_t> tempAlgoDims(this->BSalgoDims);

      for (size_t i = 0; i < this->BSalgoDims.size(); i++) {
        size_t dimCount = 0;

        for (size_t j = 0; j < tempAlgoDims.size(); j++) {
          if (this->BSalgoDims[i] == tempAlgoDims[j]) {
            dimCount++;
          }
        }

        if (dimCount > 1) {
          throw sg::base::algorithm_exception("BlackScholesParabolicPDESolverSystemEuropean::BlackScholesParabolicPDESolverSystemEuropean : There is minimum one doubled algorithmic dimension!");
        }
      }

      // build the coefficient vectors for the operations
      this->gammaCoef = new sg::base::DataMatrix(this->BSalgoDims.size(), this->BSalgoDims.size());
      this->deltaCoef = new sg::base::DataVector(this->BSalgoDims.size());

      // create the inner grid
      this->GridConverter->buildInnerGridWithCoefs(*this->BoundGrid, *this->alpha_complete, &this->InnerGrid, &this->alpha_inner);

      if (bLogTransform == false) {
        buildDeltaCoefficients();
        buildGammaCoefficients();

        //Create needed operations, on inner grid
        this->OpDeltaInner = sg::op_factory::createOperationDelta(*this->InnerGrid, *this->deltaCoef);
        this->OpGammaInner = sg::op_factory::createOperationGamma(*this->InnerGrid, *this->gammaCoef);
        // Create needed operations, on boundary grid
        this->OpDeltaBound = sg::op_factory::createOperationDelta(*this->BoundGrid, *this->deltaCoef);
        this->OpGammaBound = sg::op_factory::createOperationGamma(*this->BoundGrid, *this->gammaCoef);
      }
      // create needed operations that are different in case of a log-transformed Black-Scholoes equation
      else {
        buildDeltaCoefficientsLogTransform();
        buildGammaCoefficientsLogTransform();

        // operations on boundary grid
        this->OpDeltaBound = sg::op_factory::createOperationDeltaLog(*this->BoundGrid, *this->deltaCoef);
        this->OpGammaBound = sg::op_factory::createOperationGammaLog(*this->BoundGrid, *this->gammaCoef);
        //operations on inner grid
        this->OpDeltaInner = sg::op_factory::createOperationDeltaLog(*this->InnerGrid, *this->deltaCoef);
        this->OpGammaInner = sg::op_factory::createOperationGammaLog(*this->InnerGrid, *this->gammaCoef);
      }

      // Create operations, independent bLogTransform
      this->OpLTwoInner = sg::op_factory::createOperationLTwoDotProduct(*this->InnerGrid);
      this->OpLTwoBound = sg::op_factory::createOperationLTwoDotProduct(*this->BoundGrid);

      // right hand side if System
      this->rhs = NULL;

      // set coarsen settings
      this->useCoarsen = useCoarsen;
      this->coarsenThreshold = coarsenThreshold;
      this->refineThreshold = refineThreshold;
      this->adaptSolveMode = adaptSolveMode;
      this->numCoarsenPoints = numCoarsenPoints;
      this->refineMode = refineMode;
      this->refineMaxLevel = refineMaxLevel;

      // init Number of AverageGridPoins
      this->numSumGridpointsInner = 0;
      this->numSumGridpointsComplete = 0;

      // init option type and strike
      this->dStrike = dStrike;
      this->option_type = option_type;

      // save coordinate transformations
      this->b_log_transform = bLogTransform;

#ifdef HEDGE
      sg::base::BoundingBox* grid_bb = this->BoundGrid->getBoundingBox();
      sg::base::DimensionBoundary* myBoundaries = new sg::base::DimensionBoundary[grid_bb->getDimensions()];

      for (size_t d = 0; d < grid_bb->getDimensions(); d++) {
        if (bLogTransform == true) {
          double interval_width = exp(grid_bb->getBoundary(d).rightBoundary) - exp(grid_bb->getBoundary(d).leftBoundary);
          double hedge_offset = (interval_width - (interval_width * HEDGE_WIDTH_PERCENT)) / 2.0;
          myBoundaries[d].leftBoundary = exp(grid_bb->getBoundary(d).leftBoundary) + hedge_offset;
          myBoundaries[d].rightBoundary = exp(grid_bb->getBoundary(d).rightBoundary) - hedge_offset;
        } else {
          double hedge_offset = (grid_bb->getIntervalWidth(d) - (grid_bb->getIntervalWidth(d) * HEDGE_WIDTH_PERCENT)) / 2.0;
          myBoundaries[d].leftBoundary = grid_bb->getBoundary(d).leftBoundary + hedge_offset;
          myBoundaries[d].rightBoundary = grid_bb->getBoundary(d).rightBoundary - hedge_offset;
        }

        myBoundaries[d].bDirichletLeft = true;
        myBoundaries[d].bDirichletRight = true;
      }

      sg::base::BoundingBox* myHedgeBB = new sg::base::BoundingBox(grid_bb->getDimensions(), myBoundaries);
      // hedging
      myHedge = new sg::finance::Hedging(*myHedgeBB, HEDGE_POINTS_PER_DIM, HEDGE_EPS, bLogTransform);

      delete myHedgeBB;
      delete[] myBoundaries;
#endif
    }

    BlackScholesParabolicPDESolverSystemEuroAmer::~BlackScholesParabolicPDESolverSystemEuroAmer() {
      delete this->OpDeltaBound;
      delete this->OpGammaBound;
      delete this->OpLTwoBound;
      delete this->OpDeltaInner;
      delete this->OpGammaInner;
      delete this->OpLTwoInner;
      delete this->gammaCoef;
      delete this->deltaCoef;
      delete this->BoundaryUpdate;
      delete this->GridConverter;

      if (this->InnerGrid != NULL) {
        delete this->InnerGrid;
      }

      if (this->alpha_inner != NULL) {
        delete this->alpha_inner;
      }

      if (this->rhs != NULL) {
        delete this->rhs;
      }

      delete this->alpha_complete_old;
      delete this->alpha_complete_tmp;

#ifdef HEDGE
      delete myHedge;
#endif
    }

    void BlackScholesParabolicPDESolverSystemEuroAmer::applyLOperatorComplete(sg::base::DataVector& alpha, sg::base::DataVector& result) {
      sg::base::DataVector temp(alpha.getSize());

      result.setAll(0.0);

      // Apply the riskfree rate
      if (this->r != 0.0) {
        this->OpLTwoBound->mult(alpha, temp);
        result.axpy((-1.0)*this->r, temp);
      }

      // Apply the delta method
      this->OpDeltaBound->mult(alpha, temp);
      result.add(temp);

      // Apply the gamma method
      this->OpGammaBound->mult(alpha, temp);
      result.sub(temp);
    }

    void BlackScholesParabolicPDESolverSystemEuroAmer::applyLOperatorInner(sg::base::DataVector& alpha, sg::base::DataVector& result) {
      sg::base::DataVector temp(alpha.getSize());

      result.setAll(0.0);

      // Apply the riskfree rate
      if (this->r != 0.0) {
        this->OpLTwoInner->mult(alpha, temp);
        result.axpy((-1.0)*this->r, temp);
      }

      // Apply the delta method
      this->OpDeltaInner->mult(alpha, temp);
      result.add(temp);

      // Apply the gamma method
      this->OpGammaInner->mult(alpha, temp);
      result.sub(temp);
    }

    void BlackScholesParabolicPDESolverSystemEuroAmer::applyMassMatrixComplete(sg::base::DataVector& alpha, sg::base::DataVector& result) {
      sg::base::DataVector temp(alpha.getSize());

      result.setAll(0.0);

      // Apply the mass matrix
      this->OpLTwoBound->mult(alpha, temp);

      result.add(temp);
    }

    void BlackScholesParabolicPDESolverSystemEuroAmer::applyMassMatrixInner(sg::base::DataVector& alpha, sg::base::DataVector& result) {
      sg::base::DataVector temp(alpha.getSize());

      result.setAll(0.0);

      // Apply the mass matrix
      this->OpLTwoInner->mult(alpha, temp);

      result.add(temp);
    }

    void BlackScholesParabolicPDESolverSystemEuroAmer::finishTimestep() {
      this->nExecTimesteps++;

      // Replace the inner coefficients on the boundary grid
      this->GridConverter->updateBoundaryCoefs(*this->alpha_complete, *this->alpha_inner);

#ifndef NOBOUNDARYDISCOUNT

      // Adjust the boundaries with the riskfree rate
      if (this->r != 0.0) {
        if (this->tOperationMode == "ExEul" || this->tOperationMode == "AdBas") {
          this->BoundaryUpdate->multiplyBoundary(*this->alpha_complete, exp(((-1.0) * (this->r * this->TimestepSize))));
        }
      }

#endif

      // check if we are doing an American put -> handle early exercise
      if (this->option_type == "std_amer_put") {
        sg::base::OperationHierarchisation* myHierarchisation = sg::op_factory::createOperationHierarchisation(*this->BoundGrid);
        myHierarchisation->doDehierarchisation(*this->alpha_complete);
        size_t dim = this->BoundGrid->getStorage()->dim();
        sg::base::BoundingBox* myBB = new sg::base::BoundingBox(*(this->BoundGrid->getBoundingBox()));

        double* dblFuncValues = new double[dim];

        for (size_t i = 0; i < this->BoundGrid->getStorage()->size(); i++) {
          std::string coords = this->BoundGrid->getStorage()->get(i)->getCoordsStringBB(*myBB);
          std::stringstream coordsStream(coords);

          double tmp;

          // read coordinates
          for (size_t j = 0; j < dim; j++) {
            coordsStream >> tmp;

            dblFuncValues[j] = tmp;
          }

          tmp = 0.0;

          if (this->b_log_transform == true) {
            for (size_t j = 0; j < dim; j++) {
              tmp += exp(dblFuncValues[j]);
            }
          } else {
            for (size_t j = 0; j < dim; j++) {
              tmp += dblFuncValues[j];
            }
          }

          //      (*this->alpha_complete)[i] = std::max<double>((*this->alpha_complete)[i], (std::max<double>(this->dStrike-((tmp/static_cast<double>(dim))), 0.0))*exp(((-1.0)*(this->r*static_cast<double>(this->nExecTimesteps)*this->TimestepSize))));
          (*this->alpha_complete)[i] = std::max<double>((*this->alpha_complete)[i], (std::max<double>(this->dStrike - ((tmp / static_cast<double>(dim))), 0.0)));
        }

        delete[] dblFuncValues;

        myHierarchisation->doHierarchisation(*this->alpha_complete);
        delete myHierarchisation;
        delete myBB;
      }
    }

    void BlackScholesParabolicPDESolverSystemEuroAmer::coarsenAndRefine(bool isLastTimestep) {
      // add number of Gridpoints
      this->numSumGridpointsInner += this->InnerGrid->getSize();
      this->numSumGridpointsComplete += this->BoundGrid->getSize();

      if (this->useCoarsen == true && isLastTimestep == false) {
        ///////////////////////////////////////////////////
        // Start integrated refinement & coarsening
        ///////////////////////////////////////////////////

        size_t originalGridSize = this->BoundGrid->getStorage()->size();

        // Coarsen the grid
        sg::base::GridGenerator* myGenerator = this->BoundGrid->createGridGenerator();

        //std::cout << "Coarsen Threshold: " << this->coarsenThreshold << std::endl;
        //std::cout << "Grid Size: " << originalGridSize << std::endl;

        if (this->adaptSolveMode == "refine" || this->adaptSolveMode == "coarsenNrefine") {
          size_t numRefines = myGenerator->getNumberOfRefinablePoints();
          sg::base::SurplusRefinementFunctor* myRefineFunc = new sg::base::SurplusRefinementFunctor(this->alpha_complete, numRefines, this->refineThreshold);

          if (this->refineMode == "maxLevel") {
            myGenerator->refineMaxLevel(myRefineFunc, this->refineMaxLevel);
            this->alpha_complete->resizeZero(this->BoundGrid->getStorage()->size());
          }

          if (this->refineMode == "classic") {
            myGenerator->refine(myRefineFunc);
            this->alpha_complete->resizeZero(this->BoundGrid->getStorage()->size());
          }

          delete myRefineFunc;
        }

        if (this->adaptSolveMode == "coarsen" || this->adaptSolveMode == "coarsenNrefine") {
          size_t numCoarsen = myGenerator->getNumberOfRemovablePoints();
          sg::base::SurplusCoarseningFunctor* myCoarsenFunctor = new sg::base::SurplusCoarseningFunctor(this->alpha_complete, numCoarsen, this->coarsenThreshold);
          myGenerator->coarsenNFirstOnly(myCoarsenFunctor, this->alpha_complete, originalGridSize);
          delete myCoarsenFunctor;
        }

        delete myGenerator;

        ///////////////////////////////////////////////////
        // End integrated refinement & coarsening
        ///////////////////////////////////////////////////

        // rebuild the inner grid + coefficients
        this->GridConverter->rebuildInnerGridWithCoefs(*this->BoundGrid, *this->alpha_complete, &this->InnerGrid, &this->alpha_inner);
      }

#ifdef HEDGE
      std::stringstream filename_ext;
      filename_ext << ((this->nExecTimesteps)*this->TimestepSize);

      myHedge->calc_hedging(*this->BoundGrid, *this->alpha_complete, filename_ext.str());
#endif
    }

    void BlackScholesParabolicPDESolverSystemEuroAmer::startTimestep() {
#ifndef NOBOUNDARYDISCOUNT

      // Adjust the boundaries with the riskfree rate
      if (this->r != 0.0) {
        if (this->tOperationMode == "CrNic" || this->tOperationMode == "ImEul") {
          this->BoundaryUpdate->multiplyBoundary(*this->alpha_complete, exp(((-1.0) * (this->r * this->TimestepSize))));
        }
      }

#endif
    }

    void BlackScholesParabolicPDESolverSystemEuroAmer::buildGammaCoefficients() {
      size_t dim = this->BSalgoDims.size();

      for (size_t i = 0; i < dim; i++) {
        for (size_t j = 0; j < dim; j++) {
          // handle diagonal
          if (i == j) {
            this->gammaCoef->set(i, j, 0.5 * ((this->sigmas->get(this->BSalgoDims[i])*this->sigmas->get(this->BSalgoDims[j]))*this->rhos->get(this->BSalgoDims[i], this->BSalgoDims[j])));
          } else {
            this->gammaCoef->set(i, j, ((this->sigmas->get(this->BSalgoDims[i])*this->sigmas->get(this->BSalgoDims[j]))*this->rhos->get(this->BSalgoDims[i], this->BSalgoDims[j])));
          }
        }
      }
    }

    void BlackScholesParabolicPDESolverSystemEuroAmer::buildDeltaCoefficients() {
      size_t dim = this->BSalgoDims.size();
      double covar_sum = 0.0;

      for (size_t i = 0; i < dim; i++) {
        covar_sum = 0.0;

        for (size_t j = 0; j < dim; j++) {
          // handle diagonal
          if (i == j) {
            covar_sum += ((this->sigmas->get(this->BSalgoDims[i]) * this->sigmas->get(this->BSalgoDims[j])) * this->rhos->get(this->BSalgoDims[i], this->BSalgoDims[j]));
          } else {
            covar_sum += (0.5 * ((this->sigmas->get(this->BSalgoDims[i]) * this->sigmas->get(this->BSalgoDims[j])) * this->rhos->get(this->BSalgoDims[i], this->BSalgoDims[j])));
          }
        }

        this->deltaCoef->set(i, this->mus->get(this->BSalgoDims[i]) - covar_sum);
      }
    }

    void BlackScholesParabolicPDESolverSystemEuroAmer::buildGammaCoefficientsLogTransform() {
      size_t dim = this->BSalgoDims.size();

      for (size_t i = 0; i < dim; i++) {
        for (size_t j = 0; j < dim; j++) {
          // handle diagonal
          if (i == j) {
            this->gammaCoef->set(i, j, 0.5 * ((this->sigmas->get(this->BSalgoDims[i])*this->sigmas->get(this->BSalgoDims[j]))*this->rhos->get(this->BSalgoDims[i], this->BSalgoDims[j])));
          } else {
            this->gammaCoef->set(i, j, ((this->sigmas->get(this->BSalgoDims[i])*this->sigmas->get(this->BSalgoDims[j]))*this->rhos->get(this->BSalgoDims[i], this->BSalgoDims[j])));
          }
        }
      }
    }

    void BlackScholesParabolicPDESolverSystemEuroAmer::buildDeltaCoefficientsLogTransform() {
      size_t dim = this->BSalgoDims.size();

      for (size_t i = 0; i < dim; i++) {
        this->deltaCoef->set(i, this->mus->get(this->BSalgoDims[i]) - (0.5 * (this->sigmas->get(this->BSalgoDims[i])*this->sigmas->get(this->BSalgoDims[i]))));
      }
    }

  }
}
