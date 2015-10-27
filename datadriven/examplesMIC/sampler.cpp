#include <iostream>
#include <string.h>

#include "sgpp/datadriven/application/MetaLearner.hpp"
#include "sgpp/datadriven/operation/hash/simple/DatadrivenOperationCommon.hpp"

int main(int argc, char** argv) {

	int baseLevel = 7;

	//std::string fileName = "debugging.arff";
	std::string fileName = "DR5_train.arff";
	//  std::string fileName = "friedman2_90000.arff";
	//  std::string fileName = "bigger.arff";

	sg::base::RegularGridConfiguration gridConfig;
	sg::solver::SLESolverConfiguration SLESolverConfigRefine;
	sg::solver::SLESolverConfiguration SLESolverConfigFinal;
	sg::base::AdpativityConfiguration adaptConfig;

	// setup grid
	gridConfig.dim_ = 0; //dim is inferred from the data
	gridConfig.level_ = baseLevel;
	gridConfig.type_ = sg::base::Linear;

	// setup adaptivity
	adaptConfig.maxLevelType_ = false;
	adaptConfig.noPoints_ = 80;
	adaptConfig.numRefinements_ = 0;
	adaptConfig.percent_ = 200.0;
	adaptConfig.threshold_ = 0.0;

	// setup solver during refinement
	SLESolverConfigRefine.eps_ = 0;
	SLESolverConfigRefine.maxIterations_ = 5;
	SLESolverConfigRefine.threshold_ = -1.0;
	SLESolverConfigRefine.type_ = sg::solver::CG;

	// setup solver for final step
	SLESolverConfigFinal.eps_ = 0;
	SLESolverConfigFinal.maxIterations_ = 5;
	SLESolverConfigFinal.threshold_ = -1.0;
	SLESolverConfigFinal.type_ = sg::solver::CG;

	double lambda = 0.000001;

	bool verbose = true;
	SGPP::datadriven::MetaLearner learner(gridConfig, SLESolverConfigRefine,
			SLESolverConfigFinal, adaptConfig, lambda, verbose);

	//Configuration for intrisics-based streaming implementation
	SGPP::datadriven::OperationMultipleEvalConfiguration configuration(
	SGPP::datadriven::OperationMultipleEvalType::STREAMING,
	SGPP::datadriven::OperationMultipleEvalSubType::DEFAULT);

//	//Configuration for subspace-based evaluation (with subspace-skipping)
//	SGPP::datadriven::OperationMultipleEvalConfiguration configuration(
//	SGPP::datadriven::OperationMultipleEvalType::SUBSPACELINEAR,
//	SGPP::datadriven::OperationMultipleEvalSubType::COMBINED);

	// only execute learning (no comparisons or tests, for performance measurements)
    learner.learn(configuration, fileName);

    // execute learning with the specified configuration and use the implementation from base as comparison
    // result grids are compared by sampling the domain (again with a grid) and comparing the evaluated values
//    learner.learnAndCompare(configuration, fileName, 8);

	return EXIT_SUCCESS;
}
