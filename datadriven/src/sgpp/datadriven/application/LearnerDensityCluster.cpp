// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <sgpp/base/grid/Grid.hpp>
#include "LearnerDensityCluster.hpp"
#include <sgpp/solver/sle/ConjugateGradients.hpp>
#include <sgpp/pde/operation/PdeOpFactory.hpp>
#include <sgpp/base/operation/hash/OperationMatrix.hpp>
#include <sgpp/base/operation/hash/OperationIdentity.hpp>
#include <sgpp/base/exception/data_exception.hpp>
#include <sgpp/base/datatypes/DataVector.hpp>
#include <sgpp/base/operation/hash/OperationMultipleEval.hpp>
#include <sgpp/base/operation/BaseOpFactory.hpp>
#include <sgpp/base/exception/factory_exception.hpp>

#include <iostream>
#include <algorithm>


#include "alglib/stdafx.h"
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "alglib/alglibmisc.h"

using namespace alglib;

#include <iostream>
#include <algorithm>
#include <utility>

#include "Graph.hpp"
#include <iostream>
#include <string>
#include <stdio.h>
#include <time.h>
#include <limits>

#include <sgpp/globaldef.hpp>


namespace SGPP {

  namespace datadriven {

#if USE_DOUBLE_PRECISION==1

    LearnerDensityCluster::LearnerDensityCluster() : LearnerBase(false, false) {

    }

    LearnerDensityCluster::LearnerDensityCluster(bool isVerbose) : LearnerBase(false, isVerbose) {

    }

    LearnerDensityCluster::~LearnerDensityCluster() {
      delete minimumPoint_;
      delete clusterAssignments_;
      clusterAssignments_ = NULL;
      delete gridVals_;
      gridVals_ = NULL;
      deleteNeighbors();
    }

    void LearnerDensityCluster::deleteNeighbors() {
      if (neighbors_ != NULL) {
        for (int i = 0; i < neighborsRows_; i++) {
          delete [] neighbors_[i];
          neighbors_[i] = NULL;
        }
      }

      delete [] neighbors_;
      neighbors_ = NULL;
    }

    const std::string currentDateTime() {
      time_t     now = time(0);
      struct tm  tstruct;
      char       buf[80];
      tstruct = *localtime(&now);
      // Visit http://en.cppreference.com/w/cpp/chrono/c/strftime
      // for more information about date/time format
      strftime(buf, sizeof(buf), "%Y-%m-%d.%X", &tstruct);

      return buf;
    }


    LearnerTiming LearnerDensityCluster::train(SGPP::base::DataMatrix& testDataset,  SGPP::base::DataVector& classes,
        const SGPP::base::RegularGridConfiguration& GridConfig, const SGPP::solver::SLESolverConfiguration& SolverConfig,
        const float_t lambda) {

      LearnerTiming result;
      result.timeComplete_ = 0.0;
      result.timeMultComplete_ = 0.0;
      result.timeMultCompute_ = 0.0;
      result.timeMultTransComplete_ = 0.0;
      result.timeMultTransCompute_ = 0.0;
      result.timeRegularization_ = 0.0;
      result.GFlop_ = 0.0;
      result.GByte_ = 0.0;


      calculateGridValues(testDataset, GridConfig, SolverConfig, lambda);
      cluster(testDataset, GridConfig);

      return result;
    }

    void LearnerDensityCluster::calculateGridValues(SGPP::base::DataMatrix& testDataset, const SGPP::base::RegularGridConfiguration& GridConfig,
        const SGPP::solver::SLESolverConfiguration& SolverConfig, const float_t lambda) {
      if (GridConfig.type_ == SGPP::base::Periodic) {
        grid_ = SGPP::base::Grid::createPeriodicGrid(GridConfig.dim_);
        grid_->createGridGenerator()->regular(GridConfig.level_);
      } else if (GridConfig.type_ == SGPP::base::Linear) {
        grid_ = SGPP::base::Grid::createLinearGrid(GridConfig.dim_);
        grid_->createGridGenerator()->regular(GridConfig.level_);
      } else if (GridConfig.type_ == SGPP::base::LinearBoundary) {
        grid_ = SGPP::base::Grid::createLinearBoundaryGrid(GridConfig.dim_);
        grid_->createGridGenerator()->regular(GridConfig.level_);
      } else {
        throw base::data_exception ("Not supported grid.");
      }

      //Equation: (A + lI)a = f
      //A = (phi_i,phi_j)_L2  - LTwoDot
      //l = lambda
      //I = identity matrix
      //a = alpha
      //f_i = 1/M * sum_j=1^M (phi_i(x_j))
      //M = number of data points (testDataset)
      if (isVerbose_)
        std::cout << currentDateTime() <<  " Start clustering"  << std::endl;

      SGPP::solver::SLESolver* myCG;
      myCG = new SGPP::solver::ConjugateGradients(SolverConfig.maxIterations_, SolverConfig.eps_);


      SGPP::base::OperationIdentity C_;
      SGPP::datadriven::DensitySystemMatrix DMatrix(*grid_, testDataset, C_, lambda);

      SGPP::base::DataVector rhs(grid_->getStorage()->size());
      DMatrix.generateb(rhs);

      if (isVerbose_)
        std::cout  << currentDateTime() << " System of linear equations created" << std::endl;

      delete alpha_;
      alpha_ = new SGPP::base::DataVector(grid_->getStorage()->size());
      alpha_->setAll(0);
      myCG->solve(DMatrix, *alpha_, rhs, false, isVerbose_, -1.0);

      if (isVerbose_)
        std::cout << currentDateTime() << " Equation solved" << std::endl;

      delete gridVals_;
      gridVals_ = new SGPP::base::DataVector(testDataset.getNrows());

      SGPP::base::OperationMultipleEval* MultEval = SGPP::op_factory::createOperationMultipleEval(*grid_, testDataset);
      MultEval->mult(*alpha_, *gridVals_);

      if (isVerbose_)
        std::cout << currentDateTime() << " Grid evaluated" << std::endl;

      delete myCG;
      myCG = NULL;
    }

    bool LearnerDensityCluster::pairCompare(const std::pair<int, float_t>& firstElem, const std::pair<int, float_t>& secondElem) {
      return firstElem.second > secondElem.second;
    }

    void LearnerDensityCluster::cluster(SGPP::base::DataMatrix& testDataset, const SGPP::base::RegularGridConfiguration& GridConfig) {
      ae_int_t nx = GridConfig.dim_;
      ae_int_t ny = 2;
      ae_int_t normtype = 3;
      kdtree kdt;
      real_1d_array x = real_1d_array();
      real_2d_array r = "[[]]";
      real_1d_array dist = "[]";

      if (neighbors_ == NULL) {
        float_t* points = new float_t[(GridConfig.dim_ + 2) * testDataset.getNrows()];

        for (size_t i = 0; i < testDataset.getNrows(); i++) {
          SGPP::base::DataVector point(GridConfig.dim_);
          testDataset.getRow(i, point);
          points[(GridConfig.dim_ + 2) * (i + 1) - 2] = float_t(i);
          points[(GridConfig.dim_ + 2) * (i + 1) - 1] = (*gridVals_)[i];
          std::copy(point.getPointer(), point.getPointer() + GridConfig.dim_, points + (i) * (GridConfig.dim_ + 2));
        }

        real_2d_array a = real_2d_array();
        a.setcontent(testDataset.getNrows(), GridConfig.dim_ + 2, points);
        delete [] points;
        points = NULL;

        kdtreebuild(a, nx, ny, normtype, kdt);
      }

      Graph G(int(testDataset.getNrows()));

      bool (LearnerDensityCluster::*thresholdFunction)(SGPP::base::DataMatrix & testDataset, int, int);
      thresholdFunction = &LearnerDensityCluster::constantThreshold;

      switch (thresholdType) {
        case SGPP::datadriven::Constant:
          thresholdFunction = &LearnerDensityCluster::constantThreshold;
          break;

        case SGPP::datadriven::Relative:
          thresholdFunction = &LearnerDensityCluster::relativeThreshold;
          break;

        case SGPP::datadriven::Difference:
          thresholdFunction = &LearnerDensityCluster::differenceThreshold;
          break;

        case SGPP::datadriven::Minima:
          thresholdFunction = &LearnerDensityCluster::relativeThreshold;
          break;

        default:
          break;
      }

      delete minimumPoint_;
      minimumPoint_ = new std::vector<bool>(testDataset.getNrows());

      if (thresholdType == SGPP::datadriven::Minima) {
        for (size_t i = 0; i < testDataset.getNrows(); i++) {
          if (neighbors_ == NULL) {
            SGPP::base::DataVector point(GridConfig.dim_);
            testDataset.getRow(i, point);

            x.setcontent(GridConfig.dim_, point.getPointer());
            kdtreequeryknn(kdt, x, numberOfNeighbors + 1);
            kdtreequeryresultsxy(kdt, r);
            kdtreequeryresultsdistances(kdt, dist);
          } else {
            float_t* neighbors = new float_t [(numberOfNeighbors + 1) * (GridConfig.dim_ + 1)];

            for (int k = 0; k < numberOfNeighbors + 1; k++) {
              neighbors[(GridConfig.dim_ + 1) * (k + 1) - 1] = neighbors_[i][k];
            }

            r.setcontent(numberOfNeighbors + 1, GridConfig.dim_ + 1, neighbors);

            delete [] neighbors;
            neighbors = NULL;
          }

          std::vector<std::pair<int, float_t> > densityNeighbours;
          // check the threshold

          for (int j = 1; j < numberOfNeighbors + 1; j++) {
            densityNeighbours.push_back(std::make_pair(j, gridVals_->get(int(r[j][GridConfig.dim_]))));
          }

          std::sort(densityNeighbours.begin(), densityNeighbours.end(), (LearnerDensityCluster::pairCompare));
          float_t min = densityNeighbours.at(numberOfNeighbors - 1).second;
          float_t max = densityNeighbours.at(0).second;

          if (gridVals_->get(i) < min)
            min = gridVals_->get(i);

          if (gridVals_->get(i) > max)
            max = gridVals_->get(i);

          minimumPoint_->at(i) = minimumPoint_->at(i) =  (gridVals_->get(i) - min) / (max - min) < threshold;
        }
      }

      for (size_t i = 0; i < testDataset.getNrows(); i++) {
        if (neighbors_ == NULL) {
          SGPP::base::DataVector point(GridConfig.dim_);
          testDataset.getRow(i, point);

          x.setcontent(GridConfig.dim_, point.getPointer());
          kdtreequeryknn(kdt, x, numberOfNeighbors + 1);
          kdtreequeryresultsxy(kdt, r);
          kdtreequeryresultsdistances(kdt, dist);
        } else {

          float_t* neighbors = new float_t [(numberOfNeighbors + 1) * (GridConfig.dim_ + 1)];

          for (int k = 0; k < numberOfNeighbors + 1; k++) {
            neighbors[(GridConfig.dim_ + 1) * (k + 1) - 1] = neighbors_[i][k];
          }

          r.setcontent(numberOfNeighbors + 1, GridConfig.dim_ + 1, neighbors);

          delete [] neighbors;
          neighbors = NULL;
        }

        if (thresholdType == SGPP::datadriven::Minima && minimumPoint_->at(i)) {
          for (int j = 1; j < numberOfNeighbors + 1; j++) {
            if (minimumPoint_->at(int(r[j][GridConfig.dim_])))
              continue;

            G.addEdge(static_cast<int>(i), int(r[j][GridConfig.dim_]));
            break;
          }
        } else {

          int numberOfAddedEdges = 0;

          // check the threshold
          for (int j = 1; j < numberOfNeighbors + 1; j++) {
            if (thresholdType == SGPP::datadriven::Minima && minimumPoint_->at(int(r[j][GridConfig.dim_])))
              break;

            if ((*this.*thresholdFunction)(testDataset, int(i), int(r[j][GridConfig.dim_]))) {
              G.addEdge(int(r[0][GridConfig.dim_]), int(r[j][GridConfig.dim_]));
              numberOfAddedEdges ++;
            }
          }

          if (numberOfAddedEdges == 0)
            G.addEdge(int(r[0][GridConfig.dim_]), int(r[1][GridConfig.dim_]));
        }

        G.addEdge(int(i), int(i));
      }


      if (isVerbose_)
        std::cout << currentDateTime() << " Graph created" << std::endl;

      std::vector<int> component = G.getComponents();

      if (isVerbose_)
        std::cout << currentDateTime() << " Connected components calculated"  << std::endl;


      delete clusterAssignments_;
      clusterAssignments_ = new SGPP::base::DataVector(component);
    }


    bool LearnerDensityCluster::constantThreshold(SGPP::base::DataMatrix& testDataset, int i, int j) {
      return (gridVals_->get(i) >= eps || gridVals_->get(j) >= eps);
    }
    bool LearnerDensityCluster::relativeThreshold(SGPP::base::DataMatrix& testDataset, int i, int j) {
      return (std::abs(std::min(gridVals_->get(i), gridVals_->get(j))) / std::abs(std::max(gridVals_->get(i), gridVals_->get(j))) >= eps);
    }
    bool LearnerDensityCluster::differenceThreshold(SGPP::base::DataMatrix& testDataset, int i, int j) {
      return (std::abs(gridVals_->get(i) - gridVals_->get(j)) <= eps);
    }

    SGPP::datadriven::DMSystemMatrixBase* LearnerDensityCluster::createDMSystem(SGPP::base::DataMatrix& trainDataset, float_t lambda) {
      throw base::factory_exception("createDMSystem is not implemented yet.");
    }

    SGPP::base::DataVector* LearnerDensityCluster::getClusterAssignments() {
      return clusterAssignments_;
    }

    void LearnerDensityCluster::setClusterConfiguration(const SGPP::datadriven::DensityBasedClusteringConfiguration& ClusterConfig) {
      eps = ClusterConfig.eps;
      numberOfNeighbors = ClusterConfig.numberOfNeighbors;
      thresholdType = ClusterConfig.thresholdType;
      threshold = ClusterConfig.threshold;
    }

    void LearnerDensityCluster::precalculateGridValues(const char* filename, SGPP::base::DataMatrix& testDataset, const SGPP::base::RegularGridConfiguration& GridConfig,
        const SGPP::solver::SLESolverConfiguration& SolverConfig, const float_t lambda) {
      delete gridVals_;
      gridVals_ = NULL;
      calculateGridValues(testDataset, GridConfig, SolverConfig, lambda);
      saveArray(filename, int(gridVals_->getSize()), gridVals_->getPointer());

      if (isVerbose_)
        std::cout << currentDateTime() << " Grid values saved" << std::endl;
    }

    void LearnerDensityCluster::loadPrecalculatedValues(const char* filename) {

      int len  = 0;
      float_t* data = loadArray(filename, &len);
      delete gridVals_;
      gridVals_ = NULL;
      gridVals_ = new SGPP::base::DataVector(data, len);

      if (isVerbose_)
        std::cout << currentDateTime() << " Grid values loaded" << std::endl;

      delete [] data;
      data =  NULL;
    }

    void LearnerDensityCluster::precalculateNeighbors(const char* filename, SGPP::base::DataMatrix& testDataset, int n) {
      float_t* points = new float_t[(testDataset.getNcols() + 2) * testDataset.getNrows()];

      for (size_t i = 0; i < testDataset.getNrows(); i++) {
        SGPP::base::DataVector point(testDataset.getNcols());
        testDataset.getRow(i, point);
        points[(testDataset.getNcols() + 2) * (i + 1) - 2] = float_t(i);
        //points[(testDataset.getNcols()+2)*(i+1)-1] = (*gridVals_)[i];
        std::copy(point.getPointer(), point.getPointer() + testDataset.getNcols(), points + (i) * (testDataset.getNcols() + 2));
      }

      ae_int_t nx = testDataset.getNcols();
      ae_int_t ny = 2;
      ae_int_t normtype = 3;
      kdtree kdt;
      real_1d_array x = real_1d_array();
      real_2d_array r = "[[]]";
      real_1d_array dist = "[]";

      real_2d_array a = real_2d_array();
      a.setcontent(testDataset.getNrows(), testDataset.getNcols() + 2, points);
      delete points;
      points = NULL;

      kdtreebuild(a, nx, ny, normtype, kdt);
      deleteNeighbors();//delete neighbors_; //TODO
      neighborsRows_ = static_cast<int>(testDataset.getNrows());
      neighborsCols_ = n;
      delete neighbors_;
      neighbors_ = new int* [testDataset.getNrows()];

      for (size_t i = 0; i < testDataset.getNrows(); i++) {
        SGPP::base::DataVector point(testDataset.getNcols());
        testDataset.getRow(i, point);

        x.setcontent(testDataset.getNcols(), point.getPointer());
        kdtreequeryknn(kdt, x, n);
        kdtreequeryresultsxy(kdt, r);
        kdtreequeryresultsdistances(kdt, dist);
        neighbors_[i] = new int[n];

        for (int j = 0; j < n; j++) {
          neighbors_[i][j] = int(r[j][testDataset.getNcols()]);
        }
      }

      saveArray(filename, static_cast<int>(testDataset.getNrows()), n, neighbors_);

      if (isVerbose_)
        std::cout << currentDateTime() << " Neighbors saved" << std::endl;
    }

    void LearnerDensityCluster::loadPrecalculatedNeighbors(const char* filename) {
      int row = 0;
      int col = 0;
      deleteNeighbors();//delete neighbors_; //TODO
      neighbors_ = loadArray(filename, &row, &col);
      neighborsRows_ = row;
      neighborsCols_ = col;

      if (isVerbose_)
        std::cout << currentDateTime() << " Neighbors loaded" << std::endl;
    }


    void LearnerDensityCluster::saveArray(const char* tFilename, int len, float_t data[]) {
      fstream f;
      f.open(tFilename, std::ios::out);
      f.write((char*) &len, sizeof(int));
      f.write((char*) data, sizeof(float_t) * (len));

      f.close();

    }

    float_t* LearnerDensityCluster::loadArray(const char* tFilename, int* len) {
      fstream f;
      f.open(tFilename, std::ios::in);
      f.read((char*) len, sizeof(int));
      float_t* data = new float_t[*len];
      f.read((char*) data, sizeof(float_t) * (*len));
      f.close();
      return data;
    }

    void LearnerDensityCluster::saveArray(const char* tFilename, int row, int col, int** data) {
      fstream f;
      f.open(tFilename, std::ios::out);
      f.write((char*) &row, sizeof(int));
      f.write((char*) &col, sizeof(int));

      for (int i = 0; i < row; i++) {
        f.write((char*) data[i], sizeof(int) * (col));
      }

      f.close();
    }

    int** LearnerDensityCluster::loadArray(const char* tFilename, int* row, int* col) {
      fstream f;
      f.open(tFilename, std::ios::in);
      f.read((char*) row, sizeof(int));
      f.read((char*) col, sizeof(int));
      int** data = new int* [*row];

      for (int i = 0; i < *row; i++)
        data[i] = new int[*col];

      for (int i = 0; i < *row; i++) {
        f.read((char*) data[i], sizeof(int) * (*col));
      }

      f.close();
      return data;
    }

    void LearnerDensityCluster::postprocessing(SGPP::base::DataMatrix& trainDataset, SGPP::base::DataVector& components, int n, SGPP::base::DataVector& newComponents) {
      //SGPP::base::DataVector components = getClusterAssignments();
      int numberOfClusters = int(components.max());

      if (n == 0 ||  n > numberOfClusters + 1)
        newComponents = components;

      int* dict = new int[numberOfClusters + 1];

      for (int i = 0; i < numberOfClusters + 1; i++) {
        dict[i] = 0;
      }

      for (size_t i = 0; i < components.getSize(); i++) {
        dict[static_cast<int>(components[i])]++;
      }

      int* sortedDict = new int[numberOfClusters + 1];
      std::memcpy(sortedDict, dict, sizeof(int) * (numberOfClusters + 1));
      std::sort(sortedDict, sortedDict + numberOfClusters + 1, std::greater<int>());

      int threshold = sortedDict[n] + 1;
      delete [] sortedDict;


      bool* aboveThreshold = new bool[trainDataset.getNrows()];
      std::vector<int>* indexUnderThreshold = new std::vector<int>();

      for (size_t i = 0; i < components.getSize(); i++) {
        if (dict[static_cast<int>(components[i])] < threshold) {
          indexUnderThreshold->push_back(static_cast<int>(i));
        }

        aboveThreshold[i] = dict[int(components[i])] >= threshold;
      }

      newComponents = components;

      for (size_t i = 0; i < indexUnderThreshold->size(); i++) {
        SGPP::base::DataVector point(trainDataset.getNcols());
        trainDataset.getRow(indexUnderThreshold->at(i), point);


        int index = indexUnderThreshold->at(i);
        int indexNeighbor = 0; //r[1][trainDataset.getNcols()];//neighbors_[i][1];//

        for (int j = 1; j < neighborsCols_; j++) {
          if (aboveThreshold[neighbors_[index][j]]) {
            indexNeighbor = neighbors_[index][j];
            break;
          }
        }

        newComponents.set(index, components.get(indexNeighbor));
      }


      delete [] aboveThreshold;
      aboveThreshold = NULL;
      delete indexUnderThreshold;
      indexUnderThreshold = NULL;
    }
#endif
  }
}
