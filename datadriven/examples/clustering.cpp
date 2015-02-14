// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at 
// sgpp.sparsegrids.org

#include <iostream>
// All SG++ headers
//#include "sgpp_base.hpp"

// Or, better!, include only those that are required
#include "base/datatypes/DataVector.hpp"
#include "base/datatypes/DataMatrix.hpp"
#include "base/grid/Grid.hpp"
#include "base/grid/GridStorage.hpp"
#include "base/grid/generation/GridGenerator.hpp"
#include "base/operation/OperationEval.hpp"
#include "base/operation/BaseOpFactory.hpp"
#include "datadriven/application/LearnerDensityCluster.hpp"

using namespace std;
using namespace SGPP;
using namespace SGPP::base;



int main() {
#if USE_DOUBLE_PRECISION==1

  SGPP::datadriven::LearnerDensityCluster* clust = new SGPP::datadriven::LearnerDensityCluster(false);

  float_t raw_data[] = { 0.01, 0, 0.3, 0.7, 0.2, 0.78,  0.6, 0.82,  0.71, 0.18, 0.5};
  DataMatrix data(raw_data, 10, 1);

  DataVector classes(14);
  classes.setAll(0);
  SGPP::base::RegularGridConfiguration gridConf;
  gridConf.dim_ = 1;
  gridConf.level_ = 4;
  gridConf.type_ = SGPP::base::Periodic;

  SGPP::solver::SLESolverConfiguration solvConf;
  solvConf.eps_ = 0.0001;
  solvConf.maxIterations_ = 100;
  solvConf.threshold_ = -1.0;
  solvConf.type_ = SGPP::solver::CG;


  SGPP::datadriven::DensityBasedClusteringConfiguration clustConf;
  clustConf.eps = 10;
  clustConf.numberOfNeighbors = 4;
  clustConf.thresholdType = SGPP::datadriven::Constant;

  clust->setClusterConfiguration(clustConf);
  clust->train(data, classes, gridConf, solvConf, 0.01);

  DataVector res = clust->getClusterAssignments();
  delete clust;
  std::cout << res.toString() << std::endl;
#else
  std::cout << __FILE__ << ": This example only builds if SG++ uses double precision floating point numbers." << std::endl;
#endif
}
