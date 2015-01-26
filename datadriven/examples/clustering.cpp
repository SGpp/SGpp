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
using namespace sg;
using namespace sg::base;



int main() {

	sg::datadriven::LearnerDensityCluster* clust = new sg::datadriven::LearnerDensityCluster(false);

	double raw_data[] = { 0.01,0,0.3, 0.7,0.2, 0.78,  0.6,0.82,  0.71,0.18,0.5};
	DataMatrix data(raw_data, 10,1);

	DataVector classes(14);
	classes.setAll(0);
	sg::base::RegularGridConfiguration gridConf;
	gridConf.dim_ = 1;
	gridConf.level_ = 4;
	gridConf.type_ = sg::base::Periodic;

	sg::solver::SLESolverConfiguration solvConf;
	solvConf.eps_ = 0.0001;
	solvConf.maxIterations_ = 100;
	solvConf.threshold_ = -1.0;
	solvConf.type_ = sg::solver::CG;


	sg::datadriven::DensityBasedClusteringConfiguration clustConf;
	clustConf.eps = 10;
	clustConf.numberOfNeighbors = 4;
	clustConf.thresholdType = sg::datadriven::Constant;

	clust->setClusterConfiguration(clustConf);
	clust->train(data,classes,gridConf,solvConf,0.01);

	DataVector res = clust->getClusterAssignments();
 	delete clust;
 	std::cout << res.toString() << std::endl;
}
