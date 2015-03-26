#include <iostream>

#include "sgpp/datadriven/application/MetaLearner.hpp"
#include "sgpp/datadriven/operation/hash/DatadrivenOperationCommon.hpp"

int main(int argc, char **argv) {

//	int maxLevel = 9;
	int maxLevel = 9;

//	std::string fileName = "debugging.arff";
	//std::string fileName = "friedman_4d_2000.arff";
	//std::string fileName = "friedman2_90000.arff";
	std::string fileName = "bigger.arff";

	//sg::base::RegularGridConfiguration gridConfig;
	sg::solver::SLESolverConfiguration SLESolverConfigRefine;
	sg::solver::SLESolverConfiguration SLESolverConfigFinal;
	sg::base::AdpativityConfiguration adaptConfig;

	// Set Adaptivity
	adaptConfig.maxLevelType_ = false;
	adaptConfig.noPoints_ = 80;
	adaptConfig.numRefinements_ = 0;
	adaptConfig.percent_ = 200.0;
	adaptConfig.threshold_ = 0.0;

	// Set solver during refinement
	SLESolverConfigRefine.eps_ = 0;
	SLESolverConfigRefine.maxIterations_ = 0;
	SLESolverConfigRefine.threshold_ = -1.0;
	SLESolverConfigRefine.type_ = sg::solver::CG;

	// Set solver for final step
	SLESolverConfigFinal.eps_ = 0;
	SLESolverConfigFinal.maxIterations_ = 10;
	SLESolverConfigFinal.threshold_ = -1.0;
	SLESolverConfigFinal.type_ = sg::solver::CG;

	std::string metaInformation = "refine: "
			+ std::to_string((unsigned long long) adaptConfig.numRefinements_)
			+ " points: "
			+ std::to_string((unsigned long long) adaptConfig.noPoints_)
			+ " iterations: "
			+ std::to_string(
					(unsigned long long) SLESolverConfigRefine.maxIterations_);

	double lambda = 0.000001;

	bool verbose = true;
	SGPP::datadriven::MetaLearner learner(SLESolverConfigRefine,
			SLESolverConfigFinal, adaptConfig, maxLevel, lambda, verbose);

	//learner.learn(kernelType, fileName);
	//learner.learnReference(fileName);

	//buggy are:
	//subspace simple - 0
	//subspacelinear combined - 60
	//streaming default - 1600 (13 without avx)
	//streaming ocl - 13

	SGPP::datadriven::OperationMultipleEvalConfiguration configuration;
	configuration.type = SGPP::datadriven::OperationMultipleEvalType::STREAMING;
	configuration.subType = SGPP::datadriven::OperationMultipleEvalSubType::OCL;
	learner.learn(configuration, fileName);
	//learner.learnReference(fileName);

	//learner.learnAndTest(fileName, testFileName, isBinaryClassificationProblem);
	//learner.learnAndCompare(configuration, fileName, 4, 1.0);
	//learner.writeStatisticsFile("statistics.csv", "test");

	return EXIT_SUCCESS;
}
