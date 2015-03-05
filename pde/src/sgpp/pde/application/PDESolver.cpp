// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <sgpp/pde/application/PDESolver.hpp>
#include <sgpp/base/grid/Grid.hpp>
#include <sgpp/base/exception/application_exception.hpp>
#include <sgpp/base/tools/StdNormalDistribution.hpp>
#include <sgpp/base/grid/generation/functors/SurplusRefinementFunctor.hpp>
#include <sgpp/base/grid/generation/functors/SurplusCoarseningFunctor.hpp>
#include <sgpp/base/operation/BaseOpFactory.hpp>

#include <sstream>
#include <iostream>
#include <string>
#include <vector>
#include <fstream>
#include <cmath>

#include <sgpp/globaldef.hpp>


namespace SGPP {
  namespace pde {

    PDESolver::PDESolver(): levels(0), dim(0), myBoundingBox(nullptr), myGridStorage(nullptr), myGrid(nullptr) {
      //initializers may be wrong - David
      bGridConstructed = false;
    }

    PDESolver::~PDESolver() {
      if (bGridConstructed) {
        delete myGrid;
      }
    }

    void PDESolver::getGridNormalDistribution(SGPP::base::DataVector& alpha, std::vector<float_t>& norm_mu,
        std::vector<float_t>& norm_sigma) {
      if (bGridConstructed) {
        float_t tmp;
        float_t value;
        SGPP::base::StdNormalDistribution myNormDistr;

        for (size_t i = 0; i < this->myGrid->getStorage()->size(); i++) {
          std::string coords = this->myGridStorage->get(i)->getCoordsStringBB(*(this->myBoundingBox));
          std::stringstream coordsStream(coords);

          value = 1.0;

          for (size_t j = 0; j < this->dim; j++) {
            coordsStream >> tmp;

            value *= myNormDistr.getDensity(tmp, norm_mu[j], norm_sigma[j]);
          }

          alpha[i] = value;
        }
      } else {
        throw new SGPP::base::application_exception(
          "PDESolver::getGridNormalDistribution : The grid wasn't initialized before!");
      }
    }

    void PDESolver::deleteGrid() {
      if (bGridConstructed) {
        delete myGrid;
        bGridConstructed = false;
        myBoundingBox = NULL;
        myGridStorage = NULL;
      } else {
        throw new SGPP::base::application_exception("PDESolver::deleteGrid : The grid wasn't initialized before!");
      }
    }

    void PDESolver::setGrid(const std::string& serializedGrid) {
      if (bGridConstructed) {
        delete myGrid;
        bGridConstructed = false;
        myBoundingBox = NULL;
        myGridStorage = NULL;
      }

      myGrid = SGPP::base::Grid::unserialize(serializedGrid);

      myBoundingBox = myGrid->getBoundingBox();
      myGridStorage = myGrid->getStorage();

      dim = myGrid->getStorage()->dim();
      levels = 0;

      bGridConstructed = true;
    }

    std::string PDESolver::getGrid() const {
      std::string gridSer = "";

      if (bGridConstructed) {
        // Serialize the grid
        myGrid->serialize(gridSer);
      } else {
        throw new SGPP::base::application_exception("PDESolver::getGrid : The grid wasn't initialized before!");
      }

      return gridSer;
    }

    void PDESolver::refineInitialGridSurplus(SGPP::base::DataVector& alpha, int numRefinePoints, float_t dThreshold) {
      size_t nRefinements;

      if (numRefinePoints < 0) {
        nRefinements = myGrid->createGridGenerator()->getNumberOfRefinablePoints();
      } else {
        nRefinements = numRefinePoints;
      }

      if (bGridConstructed) {
        SGPP::base::SurplusRefinementFunctor* myRefineFunc = new SGPP::base::SurplusRefinementFunctor(&alpha, nRefinements,
            dThreshold);

        myGrid->createGridGenerator()->refine(myRefineFunc);

        delete myRefineFunc;

        alpha.resize(myGridStorage->size());
      } else {
        throw new SGPP::base::application_exception(
          "PDESolver::refineIntialGridSurplus : The grid wasn't initialized before!");
      }
    }

    void PDESolver::refineInitialGridSurplusSubDomain(SGPP::base::DataVector& alpha, int numRefinePoints, float_t dThreshold,
        std::vector<float_t>& norm_mu, std::vector<float_t>& norm_sigma) {
      size_t nRefinements;

      if (numRefinePoints < 0) {
        nRefinements = myGrid->createGridGenerator()->getNumberOfRefinablePoints();
      } else {
        nRefinements = numRefinePoints;
      }

      if (bGridConstructed) {
        SGPP::base::DataVector stdNormDist(alpha.getSize());

        // calculate multidimensional normal distribution and apply to alpha on it
        this->getGridNormalDistribution(stdNormDist, norm_mu, norm_sigma);
        //printSparseGrid(stdNormDist, "normalDistribution.grid.gnuplot", true);
        stdNormDist.componentwise_mult(alpha);
        //printSparseGrid(stdNormDist, "normalDistribution_refine.grid.gnuplot", true);

        SGPP::base::SurplusRefinementFunctor* myRefineFunc = new SGPP::base::SurplusRefinementFunctor(&stdNormDist,
            nRefinements, dThreshold);

        myGrid->createGridGenerator()->refine(myRefineFunc);

        delete myRefineFunc;

        alpha.resize(myGridStorage->size());
      } else {
        throw new SGPP::base::application_exception(
          "PDESolver::refineIntialGridSurplusSubDomain : The grid wasn't initialized before!");
      }
    }

    void PDESolver::refineInitialGridSurplusToMaxLevel(SGPP::base::DataVector& alpha, float_t dThreshold,
        SGPP::base::GridStorage::index_type::level_type maxLevel) {
      if (bGridConstructed) {
        size_t nRefinements = myGrid->createGridGenerator()->getNumberOfRefinablePointsToMaxLevel(maxLevel);

        SGPP::base::SurplusRefinementFunctor* myRefineFunc = new SGPP::base::SurplusRefinementFunctor(&alpha, nRefinements,
            dThreshold);

        myGrid->createGridGenerator()->refineMaxLevel(myRefineFunc, maxLevel);

        delete myRefineFunc;

        alpha.resize(myGridStorage->size());
      } else {
        throw new SGPP::base::application_exception(
          "PDESolver::refineInitialGridSurplusToMaxLevel : The grid wasn't initialized before!");
      }
    }

    void PDESolver::refineInitialGridSurplusToMaxLevelSubDomain(SGPP::base::DataVector& alpha, float_t dThreshold,
        SGPP::base::GridStorage::index_type::level_type maxLevel, std::vector<float_t>& norm_mu,
        std::vector<float_t>& norm_sigma) {
      if (bGridConstructed) {
        size_t nRefinements = myGrid->createGridGenerator()->getNumberOfRefinablePointsToMaxLevel(maxLevel);

        SGPP::base::DataVector stdNormDist(alpha.getSize());

        // calculate multidimensional normal distribution and apply to alpha on it
        this->getGridNormalDistribution(stdNormDist, norm_mu, norm_sigma);
        //printSparseGrid(stdNormDist, "normalDistribution.grid.gnuplot", true);
        stdNormDist.componentwise_mult(alpha);
        //printSparseGrid(stdNormDist, "normalDistribution_refine.grid.gnuplot", true);

        SGPP::base::SurplusRefinementFunctor* myRefineFunc = new SGPP::base::SurplusRefinementFunctor(&stdNormDist,
            nRefinements, dThreshold);

        myGrid->createGridGenerator()->refineMaxLevel(myRefineFunc, maxLevel);

        delete myRefineFunc;

        alpha.resize(myGridStorage->size());
      } else {
        throw new SGPP::base::application_exception(
          "PDESolver::refineInitialGridSurplusToMaxLevelSubDomain : The grid wasn't initialized before!");
      }
    }

    void PDESolver::coarsenInitialGridSurplus(SGPP::base::DataVector& alpha, float_t dThreshold) {
      if (bGridConstructed) {
        SGPP::base::GridGenerator* myGenerator = myGrid->createGridGenerator();
        size_t numCoarsen = myGenerator->getNumberOfRemovablePoints();
        size_t originalGridSize = myGrid->getStorage()->size();
        SGPP::base::SurplusCoarseningFunctor* myCoarsenFunctor = new SGPP::base::SurplusCoarseningFunctor(&alpha,
            numCoarsen, dThreshold);

        myGenerator->coarsenNFirstOnly(myCoarsenFunctor, &alpha, originalGridSize);

        delete myCoarsenFunctor;
        delete myGenerator;
      } else {
        throw new SGPP::base::application_exception(
          "PDESolver::coarsenInitialGridSurplus : The grid wasn't initialized before!");
      }
    }

    void PDESolver::printLevelIndexGrid(std::string tfilename) const {
      SGPP::base::GridPrinter myPrinter(*this->myGrid);
      myPrinter.printLevelIndexGrid(tfilename);
    }

    void PDESolver::printGrid(SGPP::base::DataVector& alpha, float_t PointesPerDimension, std::string tfilename) const {
      SGPP::base::GridPrinter myPrinter(*this->myGrid);
      myPrinter.printGrid(alpha, tfilename, static_cast<size_t>(PointesPerDimension));
    }

    void PDESolver::printGridDomain(SGPP::base::DataVector& alpha, float_t PointesPerDimension,
                                    SGPP::base::BoundingBox& GridArea, std::string tfilename) const {
      SGPP::base::GridPrinter myPrinter(*this->myGrid);
      myPrinter.printGridDomain(alpha, tfilename, GridArea, static_cast<size_t>(PointesPerDimension));
    }

    void PDESolver::printSparseGrid(SGPP::base::DataVector& alpha, std::string tfilename, bool bSurplus) const {
      SGPP::base::GridPrinter myPrinter(*this->myGrid);
      myPrinter.printSparseGrid(alpha, tfilename, bSurplus);
    }

    void PDESolver::printSparseGridExpTransform(SGPP::base::DataVector& alpha, std::string tfilename, bool bSurplus) const {
      SGPP::base::GridPrinter myPrinter(*this->myGrid);
      myPrinter.printSparseGridExpTransform(alpha, tfilename, bSurplus);
    }

    float_t PDESolver::evaluatePoint(std::vector<float_t>& evalPoint, SGPP::base::DataVector& alpha) {
      float_t result = 0.0;

      if (bGridConstructed) {
        SGPP::base::OperationEval* myEval = SGPP::op_factory::createOperationEval(*myGrid);
        result = myEval->eval(alpha, evalPoint);
        delete myEval;
      } else {
        throw new SGPP::base::application_exception("PDESolver::evaluatePoint : A grid wasn't constructed before!");
      }

      return result;
    }

    void PDESolver::evaluateCuboid(SGPP::base::DataVector& alpha, SGPP::base::DataVector& OptionPrices,
                                   SGPP::base::DataMatrix& EvaluationPoints) {
      if (bGridConstructed) {
        if (OptionPrices.getSize() != EvaluationPoints.getNrows()) {
          throw new SGPP::base::application_exception(
            "PDESolver::evaluateCuboid : The size of the price vector doesn't match the size of the evaluation points' vector!");
        }

        SGPP::base::OperationMultipleEval* myOpMultEval = SGPP::op_factory::createOperationMultipleEval(*myGrid,
            EvaluationPoints);
        myOpMultEval->mult(alpha, OptionPrices);
        delete myOpMultEval;
      } else {
        throw new SGPP::base::application_exception("PDESolver::evaluateCuboid : A grid wasn't constructed before!");
      }
    }

    size_t PDESolver::getNumberGridPoints() const {
      if (bGridConstructed) {
        return myGridStorage->size();
      } else {
        throw new SGPP::base::application_exception("PDESolver::getNumberGridPoints : A grid wasn't constructed before!");
      }
    }

    size_t PDESolver::getNumberInnerGridPoints() const {
      if (bGridConstructed) {
        return myGridStorage->getNumInnerPoints();
      } else {
        throw new SGPP::base::application_exception("PDESolver::getNumberGridPoints : A grid wasn't constructed before!");
      }
    }

    size_t PDESolver::getNumberDimensions() const {
      if (bGridConstructed) {
        return myGridStorage->dim();
      } else {
        throw new SGPP::base::application_exception("PDESolver::getNumberDimensions : A grid wasn't constructed before!");
      }
    }

  }
}