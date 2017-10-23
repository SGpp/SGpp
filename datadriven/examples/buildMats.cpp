// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <sgpp/globaldef.hpp>

#include <sgpp/datadriven/algorithm/DBMatDensityConfiguration.hpp>
#include <sgpp/datadriven/application/LearnerSGDEOnOff.hpp>
#include <sgpp/datadriven/algorithm/DBMatOffline.hpp>
#include <sgpp/datadriven/algorithm/DBMatOfflineFactory.hpp>
#include <sgpp/datadriven/tools/ARFFTools.hpp>

#include <string>

using sgpp::base::DataMatrix;
using sgpp::base::DataVector;

std::vector<std::vector<size_t>> getAllInteractions(size_t res) {
	size_t geodim = res;
	std::vector<std::vector<size_t>> vec = std::vector<std::vector<size_t>>();

	for(size_t i = 0; i < geodim; i++){
		for(size_t j = 0; j < geodim-1; j++){
			std::vector<size_t> xdir = std::vector<size_t>();

			xdir.push_back(i*geodim+j);

			xdir.push_back(i*geodim+j+1);

			vec.push_back(xdir);
		}
	}for(size_t i = 0; i < geodim-1; i++){
		for(size_t j = 0; j < geodim; j++){
			std::vector<size_t> ydir = std::vector<size_t>();

			ydir.push_back(i*geodim+j);

			ydir.push_back((i+1)*geodim+j);

			vec.push_back(ydir);
		}
	}
	//1d vector for all dimensions
    for(size_t i = 0; i < geodim*geodim; i++){
		std::vector<size_t> tmp = std::vector<size_t>();
		tmp.push_back(i);
		vec.push_back(tmp);
	}
	//add empty vector
	std::vector<size_t> empty = std::vector<size_t>();
	vec.push_back(empty);


	return vec;

}
std::vector<std::vector<size_t>> getConvs(size_t res) {
	size_t geodim = res;
	std::vector<std::vector<size_t>> vec = std::vector<std::vector<size_t>>();

	for(size_t i = 0; i < geodim-1; i+=2){
		for(size_t j = 0; j < geodim-1; j+=2){
			std::vector<size_t> xdir1 = std::vector<size_t>();
			std::vector<size_t> xdir2 = std::vector<size_t>();
			std::vector<size_t> ydir1 = std::vector<size_t>();
			std::vector<size_t> ydir2 = std::vector<size_t>();

			xdir1.push_back(i*geodim+j);
			ydir1.push_back(i*geodim+j);

			xdir1.push_back(i*geodim+j+1);
			ydir2.push_back(i*geodim+j+1);

			xdir2.push_back((i+1)*geodim+j);
			ydir1.push_back((i+1)*geodim+j);

			xdir2.push_back((i+1)*geodim+j+1);
			ydir2.push_back((i+1)*geodim+j+1);


			vec.push_back(xdir1);
			vec.push_back(xdir2);
			vec.push_back(ydir1);
			vec.push_back(ydir2);
		}
	}

    for(size_t i = 0; i < geodim*geodim; i++){
		std::vector<size_t> tmp = std::vector<size_t>();
		tmp.push_back(i);
		vec.push_back(tmp);
	}
	std::vector<size_t> empty = std::vector<size_t>();
	vec.push_back(empty);


	return vec;

}

/**
 * This example shows how to perform offline/online-classification using sparse
 * grid density estimation and matrix decomposition methods. It creates an
 * instance of LearnerSGDEOnOff and runs the function train() where the
 * main functionality is implemented.
 *
 * Currently, only binary classification with class labels -1 and 1 is possible.
 *
 * The example provides the option to execute several runs over differently
 * ordered data and perform a 5-fold cross-validation within each run.
 * Therefore,
 * already randomly ordered and partitioned data is required.
 * Average results from several runs might be more reliable in an
 * online-learning
 * scenario, because the ordering of the data points seen by the learner
 * can affect the result.
 */

int main() {
  int lvl = 3;
  for(size_t res = 28; res<=28; res+=2){
    std::string filename = "mats/" + std::to_string(res) + "x" + std::to_string(res) + "_ModLin_NN_Inter_lvl"+std::to_string(lvl)+"_Chol.out";
    std::cout << "Setting up " << filename << std::endl;
	/**
	* The grid configuration.
	*/
	std::cout << "# create grid config" << std::endl;
	sgpp::base::RegularGridConfiguration gridConfig;
	gridConfig.dim_ = res*res;
	gridConfig.level_ = lvl;
	//gridConfig.type_ = sgpp::base::GridType::Linear;
	gridConfig.type_ = sgpp::base::GridType::ModLinear;

	/**
	* Configure regularization.
	*/
	std::cout << "# create regularization config" << std::endl;
	sgpp::datadriven::RegularizationConfiguration regularizationConfig;
	regularizationConfig.regType_ = sgpp::datadriven::RegularizationType::Identity;

	/**
	* Select the desired decomposition type for the offline step.
	* Note: Refinement/Coarsening only possible for Cholesky decomposition.
	*/
	sgpp::datadriven::DBMatDecompostionType dt;
	std::string decompType;
	// choose "LU decomposition"
	// dt = DBMatDecompostionType::DBMatDecompLU;
	// decompType = "LU decomposition";
	// choose"Eigen decomposition"
	// dt = DBMatDecompostionType::DBMatDecompEigen;
	// decompType = "Eigen decomposition";
	// choose "Cholesky decomposition"
	dt = sgpp::datadriven::DBMatDecompostionType::Chol;
	decompType = "Cholesky decomposition";
	//      dt = sgpp::datadriven::DBMatDecompostionType::IChol;
	//      decompType = "Incomplete Cholesky decomposition";
	//		dt = sgpp::datadriven::DBMatDecompostionType::DenseIchol;
	//		decompType = "Incomplete Cholesky decomposition on Dense Matrix";
	std::cout << "Decomposition type: " << decompType << std::endl;

	/**
	* Configure adaptive refinement (if Cholesky is chosen). As refinement
	* monitor the periodic monitor or the convergence monitor
	* can be chosen. Possible refinement indicators are
	* surplus refinement, data-based refinement, zero-crossings-based
	* refinement.
	*/
	std::cout << "# create adaptive refinement configuration" << std::endl;
	std::string refMonitor;
	// select periodic monitor - perform refinements in fixed intervals
	refMonitor = "periodic";
	std::cout << "Refinement monitor: " << refMonitor << std::endl;
	std::string refType;
	// select surplus refinement
	// refType = "surplus";
	// select data-based refinement
	// refType = "data";
	// select zero-crossings-based refinement
	refType = "zero";
	std::cout << "Refinement type: " << refType << std::endl;
	sgpp::base::AdpativityConfiguration adaptConfig;
	/**
	* Specify number of refinement steps and the max number
	* of grid points to refine each step.
	*/
	adaptConfig.numRefinements_ = 0;
	adaptConfig.noPoints_ = 7;
	adaptConfig.threshold_ = 0.0;  // only required for surplus refinement

	// initial regularization parameter lambda
	double lambda = 0.01;


	// configuration
	sgpp::datadriven::DBMatDensityConfiguration dconf(gridConfig, adaptConfig,
			                                        regularizationConfig.regType_, lambda, dt);

	sgpp::datadriven::DBMatOffline *offline = sgpp::datadriven::DBMatOfflineFactory::buildOfflineObject(dconf);
    offline->setInter(getAllInteractions(res));
	std::cout << "Building Matrix..." << std::endl;
	offline->buildMatrix();
	std::cout << "Matrix build.\nBegin decomposition..." << std::endl;
	offline->decomposeMatrix();
	//offline->printMatrix();
    offline->store(filename);
  }
}


