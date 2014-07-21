#include <iostream>
#include <fstream>
#include "test/FriedmanTestData.hpp"
#include "sgpp_base.hpp"
#include "sgpp_datadriven.hpp"

/*
 * Usage:
 * main condition numPoints type
 *
 * condition: the condition for refinement
 * 0 = every 10 SGD steps
 * 1 = if smoothed error decline < 10%
 *
 * numPoints: maximal number of points to be refined
 *
 * type: type of refinement operator
 * 0 = Surplus refinement
 * 1 = Weighted error on last minibatch
 * 2 = Weighted error on complete training dataset
 * 3 = Discouting error indicator
 */

int main(int argc, char **args) {

	using namespace sg::base;
	using namespace std;

	if (argc != 4) {
		std::cout << "Not enough arguments." << std::endl;
		return 1;
	}

	/* Parameters */
	int refinementCondition = atoi(args[1]);
	int refinementNumPoints = atoi(args[2]);
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
	ferr0.open("ferr0", ios::out | ios::trunc);
	ferr1.open("ferr1", ios::out | ios::trunc);
	ferr2.open("ferr2", ios::out | ios::trunc);
	fgrid.open("fgrid", ios::out | ios::trunc);
	fcoor.open("fcoor", ios::out | ios::trunc);

	ostream *outputStreams[5] = { &ferr0, &ferr1, &ferr2, &fgrid, &fcoor };

	/* Learner instances (set verbosity here) */
	sg::datadriven::LearnerRegularizationType r = sg::datadriven::Identity;
	sg::datadriven::LearnerOnlineSGD *learner =
			new sg::datadriven::LearnerOnlineSGD(r, false, false);

	/*
	 * Friedman
	 */
	/*
	 sg::test::FriedmanTestData testData(dim, level, trainSize);
	 testData.generate(seed);

	 sg::base::DataMatrix trainData = testData.getTrainData();
	 sg::base::DataVector classes = testData.getClasses();
	 */

	sg::datadriven::ARFFTools arff = sg::datadriven::ARFFTools();

	std::string mainARFFFilename = "data/DR5_MGS_ugrizeC_train.60000.arff";
	std::string testARFFFilename = "data/DR5_MGS_ugrizeC_test.60000.arff";

	size_t mainNumRows = arff.getNumberInstances(mainARFFFilename);
	size_t mainNumCols = arff.getDimension(mainARFFFilename);

	size_t testNumRows = arff.getNumberInstances(testARFFFilename);
	size_t testNumCols = arff.getDimension(testARFFFilename);

	sg::base::DataMatrix mainTrainData = sg::base::DataMatrix(mainNumRows, mainNumCols);
	sg::base::DataVector mainClasses = sg::base::DataVector(mainNumRows);

	sg::base::DataMatrix testTrainData = sg::base::DataMatrix(testNumRows, testNumCols);
	sg::base::DataVector testClasses = sg::base::DataVector(testNumRows);

	arff.readTrainingData(mainARFFFilename, mainTrainData);
	arff.readClasses(mainARFFFilename, mainClasses);

	arff.readTrainingData(testARFFFilename, testTrainData);
	arff.readClasses(testARFFFilename, testClasses);

	sg::base::RegularGridConfiguration gridConfig;
	gridConfig.dim_ = dim;
	gridConfig.level_ = level;
	gridConfig.type_ = sg::base::ModLinear;

	// Initialize HashRefinement
	HashRefinement *hashRefinement = new HashRefinement();

	learner->train(mainTrainData, mainClasses, testTrainData, testClasses,
			gridConfig, numIterations, batchSize, lambda, gamma,
			*hashRefinement, refinementCondition, refinementType,
			refinementNumPoints, numRuns, outputStreams);

	ferr0.close();
	ferr1.close();
	ferr2.close();
	fgrid.close();
	fcoor.close();
}
