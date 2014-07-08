#include <iostream>
#include <fstream>
#include "test/FriedmanTestData.hpp"
#include "sgpp_base.hpp"
#include "sgpp_datadriven.hpp"

/*
 * Usage:
 * main condition numPoints type
 */

int main(int argc, char **args) {

	using namespace sg::base;
	using namespace std;

	if (argc != 4) {
		std::cout << "Not enough arguments." << std::endl;
		return 1;
	}

	/* Parameters */
	// 0 = every 10 SGD steps
	// 1 = if smoothed error decline < 10%
	int refinementCondition = atoi(args[1]);

	// Number of points to refine
	int refinementNumPoints = atoi(args[2]);

	// Type of refinement indicator
	// 0-3
	int refinementType = atoi(args[3]);

	/* Constants */
	size_t dim = 3;
	int level = 3;
	int trainSize = 1000;
	int seed = 1208108;

	size_t numIterations = 10;
	size_t batchSize = 10;
	double lambda = 0.00001;
	double gamma = 0.001;

	int numRuns = 1;

	/* Output files */
	fstream ferr0, ferr1, ferr2, fgrid, fcoor;
	ferr0.open ("ferr0", ios::out | ios::trunc);
	ferr1.open ("ferr1", ios::out | ios::trunc);
	ferr2.open ("ferr2", ios::out | ios::trunc);
	fgrid.open ("fgrid", ios::out | ios::trunc);
	fcoor.open ("fcoor", ios::out | ios::trunc);

	ostream *outputStreams[5] = {&ferr0, &ferr1, &ferr2, &fgrid, &fcoor};

	sg::datadriven::LearnerRegularizationType r = sg::datadriven::Identity;
	sg::datadriven::LearnerOnlineSGD *learner =
			new sg::datadriven::LearnerOnlineSGD(r, true, false);

	/*
	 * TODO
	 */
	sg::test::FriedmanTestData testData(dim, level, trainSize);
	testData.generate(seed);
	sg::base::DataMatrix trainData = testData.getTrainData();
	sg::base::DataVector classes = testData.getClasses();

	sg::base::RegularGridConfiguration gridConfig;
	gridConfig.dim_ = dim;
	gridConfig.level_ = level;
	gridConfig.type_ = sg::base::ModLinear;

	// Initialize HashRefinement
	HashRefinement *hashRefinement = new HashRefinement();

	learner->train(trainData, classes, gridConfig, numIterations, batchSize, lambda, gamma, *hashRefinement,
			refinementCondition, refinementType, refinementNumPoints, numRuns, outputStreams);

	ferr0.close();
	ferr1.close();
	ferr2.close();
	fgrid.close();
	fcoor.close();
}
