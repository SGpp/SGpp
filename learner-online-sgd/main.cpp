#include <iostream>
#include <fstream>
#include "sgpp_base.hpp"
#include "sgpp_datadriven.hpp"

const double VALIDATION_FACTOR = 0.2;

double execLearnerPartial(sg::datadriven::Learner* learner,
		sg::base::DataMatrix& mainTrainData, sg::base::DataVector& mainClasses,
		size_t part, const sg::base::RegularGridConfiguration& GridConfig,
		const sg::solver::SLESolverConfiguration& SolverConfig,
		const double lambda);

int main(int argc, char **args) {

	using namespace sg::base;
	using namespace std;

	map<string, string> argsMap;

	for (int i = 1; i < argc; i += 2) {
		argsMap[args[i]] = args[i + 1];
	}

	if (argsMap["mode"] == "LEARNER_ONLINE_SGD") {
		/*
		 * Output files
		 */
		fstream ferr0, ferr1, ferr2, fgrid, fcoor;
		ferr0.open((argsMap["experimentDir"] + std::string("/ferr0")).c_str(),
				ios::out | ios::trunc);
		ferr1.open((argsMap["experimentDir"] + std::string("/ferr1")).c_str(),
				ios::out | ios::trunc);
		ferr2.open((argsMap["experimentDir"] + std::string("/ferr2")).c_str(),
				ios::out | ios::trunc);
		fgrid.open((argsMap["experimentDir"] + std::string("/fgrid")).c_str(),
				ios::out | ios::trunc);
		fcoor.open((argsMap["experimentDir"] + std::string("/fcoor")).c_str(),
				ios::out | ios::trunc);

		ostream *outputStreams[5] = { &ferr0, &ferr1, &ferr2, &fgrid, &fcoor };

		/*
		 * Get data
		 */

		sg::datadriven::ARFFTools arff = sg::datadriven::ARFFTools();

		size_t mainNumRows = arff.getNumberInstances(
				argsMap["mainARFFFilename"]);
		size_t mainNumCols = arff.getDimension(argsMap["mainARFFFilename"]);

		size_t testNumRows = arff.getNumberInstances(
				argsMap["testARFFFilename"]);
		size_t testNumCols = arff.getDimension(argsMap["testARFFFilename"]);

		sg::base::DataMatrix mainTrainData = sg::base::DataMatrix(mainNumRows,
				mainNumCols);
		sg::base::DataVector mainClasses = sg::base::DataVector(mainNumRows);

		sg::base::DataMatrix testTrainData = sg::base::DataMatrix(testNumRows,
				testNumCols);
		sg::base::DataVector testClasses = sg::base::DataVector(testNumRows);

		arff.readTrainingData(argsMap["mainARFFFilename"], mainTrainData);
		arff.readClasses(argsMap["mainARFFFilename"], mainClasses);

		arff.readTrainingData(argsMap["testARFFFilename"], testTrainData);
		arff.readClasses(argsMap["testARFFFilename"], testClasses);

		/*
		 * Grid Configuration
		 */

		sg::base::RegularGridConfiguration gridConfig;
		gridConfig.dim_ = (size_t) atoi(argsMap["dim"].c_str());
		gridConfig.level_ = atoi(argsMap["level"].c_str());
		gridConfig.type_ = sg::base::ModLinear;

		/*
		 * Adaptivity Configuration
		 */

		sg::datadriven::LearnerOnlineSGDRefinementConfiguration refineConfig;
		refineConfig.refinementCondition = argsMap["refinementCondition"];
		refineConfig.refinementType = argsMap["refinementType"];
		refineConfig.numIterations = (size_t) atoi(
				argsMap["numIterations"].c_str());
		refineConfig.numMinibatchError = (size_t) atoi(
				argsMap["numMinibatchError"].c_str());
		refineConfig.refinementNumPoints = atoi(
				argsMap["refinementNumPoints"].c_str());

		/*
		 * Other
		 */

		double lambda_ = atof(argsMap["regularizationLambda"].c_str());
		double gamma_ = atof(argsMap["CGStepSizeGamma"].c_str());

		int batchSize = atoi(argsMap["batchSize"].c_str());
		int numRuns = atoi(argsMap["numRuns"].c_str());

		string errorType = argsMap["errorType"];

		/*
		 * Peform training
		 */
		sg::datadriven::LearnerRegularizationType r = sg::datadriven::Identity;
		sg::datadriven::LearnerOnlineSGD *learner =
				new sg::datadriven::LearnerOnlineSGD(r, false, false);

		HashRefinement* hashRef = new HashRefinement();
		learner->train(mainTrainData, mainClasses, testTrainData, testClasses,
				gridConfig, refineConfig, *hashRef, batchSize, lambda_, gamma_,
				numRuns, errorType, outputStreams);
		/*
		 * Close output files
		 */
		ferr0.close();
		ferr1.close();
		ferr2.close();
		fgrid.close();
		fcoor.close();

	} else if (argsMap["mode"] == "FIND_LAMBDA_WITH_LEARNER_BASE") {

		const int CG_IMAX = 10000;
		const double CG_EPS = 0.0001;

		/*
		 * Get data
		 */

		sg::datadriven::ARFFTools arff = sg::datadriven::ARFFTools();

		size_t mainNumRows = arff.getNumberInstances(
				argsMap["mainARFFFilename"]);
		size_t mainNumCols = arff.getDimension(argsMap["mainARFFFilename"]);

		sg::base::DataMatrix mainTrainData = sg::base::DataMatrix(mainNumRows,
				mainNumCols);
		sg::base::DataVector mainClasses = sg::base::DataVector(mainNumRows);

		arff.readTrainingData(argsMap["mainARFFFilename"], mainTrainData);
		arff.readClasses(argsMap["mainARFFFilename"], mainClasses);

		/*
		 * Configuration
		 */

		sg::base::RegularGridConfiguration GridConfig;
		GridConfig.dim_ = (size_t) atoi(argsMap["dim"].c_str());
		GridConfig.level_ = atoi(argsMap["level"].c_str());
		GridConfig.type_ = sg::base::ModLinear;

		sg::solver::SLESolverConfiguration SolverConfig;
		SolverConfig.type_ = sg::solver::SLESolverType::CG;
		SolverConfig.eps_ = CG_EPS;
		SolverConfig.maxIterations_ = CG_IMAX;
		SolverConfig.threshold_ = 10;

		double lambda_ = atof(argsMap["regularizationLambda"].c_str());

		std::cout << "Lambda: " << lambda_ << std::endl;

		/*
		 * Learn
		 */

		sg::datadriven::LearnerRegularizationType r = sg::datadriven::Identity;
		sg::datadriven::Learner* learner = new sg::datadriven::Learner(r, false,
				false);

		double accuracy = 0;

		if (argsMap["crossValidation"] == "FALSE") {
			accuracy = execLearnerPartial(learner, mainTrainData, mainClasses,
					0, GridConfig, SolverConfig, lambda_);
		} else {
			double sum = 0;
			double current_accuracy = 0;
			int numParts = floor(1 / VALIDATION_FACTOR);

			for (int i = 0; i < numParts; i++) {
				current_accuracy = execLearnerPartial(learner, mainTrainData, mainClasses,
						i, GridConfig, SolverConfig, lambda_);
				sum += current_accuracy;
				std::cout << "Accuracy (part " << (i+1) << "): " << current_accuracy << std::endl;
			}
			accuracy = sum / numParts;
		}

		std::cout << "Accuracy: " << accuracy << std::endl;
	}
}

/*
 * Cross Validation
 */

double execLearnerPartial(sg::datadriven::Learner* learner,
		sg::base::DataMatrix& trainData, sg::base::DataVector& classes,
		size_t part, const sg::base::RegularGridConfiguration& GridConfig,
		const sg::solver::SLESolverConfiguration& SolverConfig,
		const double lambda) {

	int numRows = trainData.getNrows();
	int numCols = trainData.getNcols();

	/*
	 * Split data
	 */

	size_t sizePartial = floor(numRows * (1-VALIDATION_FACTOR));
	size_t sizeValidation = numRows - sizePartial;

	sg::base::DataMatrix partialTrainData = sg::base::DataMatrix(sizePartial,
			numCols);
	sg::base::DataVector partialClasses = sg::base::DataVector(sizePartial);

	sg::base::DataMatrix validationTrainData = sg::base::DataMatrix(
			sizeValidation, numCols);

	sg::base::DataVector validationClasses = sg::base::DataVector(
			sizeValidation);

	sg::base::DataVector tmpRow = sg::base::DataVector(numCols);
	double tmpVal = 0;

	for (size_t i = 0; i < numRows; i++) {
		trainData.getRow(i, tmpRow);
		tmpVal = classes.get(i);

		if (i < part * sizeValidation) {
			partialTrainData.setRow(i, tmpRow);
			partialClasses.set(i, tmpVal);

		} else if (i < (part + 1) * sizeValidation) {
			validationTrainData.setRow(i - part * sizeValidation, tmpRow);
			validationClasses.set(i - part * sizeValidation, tmpVal);
		} else {
			partialTrainData.setRow(i - sizeValidation, tmpRow);
			partialClasses.set(i - sizeValidation, tmpVal);
		}
	}

	learner->train(partialTrainData, partialClasses, GridConfig,
			SolverConfig, lambda);

	return learner->getAccuracy(validationTrainData, validationClasses,
			0.0);
}
