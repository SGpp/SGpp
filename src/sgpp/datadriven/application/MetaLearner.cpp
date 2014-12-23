#include "MetaLearner.hpp"

#include <random>

#include <fstream>
#include <iomanip>
#include <cstdio>

#include "LearnerLeastSquaresIdentity.hpp"
#include "LearnerVectorizedSubspaces.hpp"

using namespace sg::base;

namespace sg {
namespace datadriven {

MetaLearner::MetaLearner(sg::solver::SLESolverConfiguration solverConfig,
		sg::solver::SLESolverConfiguration solverFinalStep, sg::base::AdpativityConfiguration adaptivityConfiguration,
		size_t baseLevel, double lambda, bool verbose) {
	this->csvSep = "& ";
	this->solverConfig = solverConfig;
	this->solverFinalStep = solverFinalStep;
	this->adaptivityConfiguration = adaptivityConfiguration;
	this->baseLevel = baseLevel;
	this->lambda = lambda;
	this->verbose = verbose;
	this->instances = 0;
	this->dim = 0;
}

/* Bug: leads to linker errors - don't know why
 MetaLearner::~MetaLearner() {
 if (this->myLearner != nullptr) {
 delete this->myLearner;
 }

 if (this->referenceLearner != nullptr) {
 delete this->referenceLearner;
 }
 }
 */

void MetaLearner::learn(sg::base::OperationMultipleEval *kernel, std::string datasetFileName) {

	this->dim = arffTools.getDimension(datasetFileName);
	this->instances = arffTools.getNumberInstances(datasetFileName);
	if (verbose) {
		std::cout << "instances: " << this->instances << std::endl;
	}

	// Set grid config
	this->gridConfig.dim_ = this->dim;
	this->gridConfig.type_ = sg::base::Linear;
	this->gridConfig.level_ = static_cast<int>(this->baseLevel);

	DataVector classesVector(instances);
	arffTools.readClasses(datasetFileName, classesVector);

	DataMatrix trainingData(instances, dim);
	arffTools.readTrainingData(datasetFileName, trainingData);

	bool isRegression = true; // treat everything as if it were a regression, as classification is not fully supported by Learner
	//TODO kernel must be configurable here
	LearnerVectorizedSubspaces *myLearner = new LearnerVectorizedSubspaces(kernel, isRegression, this->verbose);

	LearnerTiming timings = myLearner->train(trainingData, classesVector, this->gridConfig, this->solverConfig,
			this->solverFinalStep, this->adaptivityConfiguration, false, this->lambda);

	this->myTiming = timings;
	this->ExecTimesOnStep = myLearner->getRefinementExecTimes();

	//myLearner->dumpFunction("mygridFunction", 50);
	if (this->myLearner != nullptr) {
		delete this->myLearner;
	}
	this->myLearner = myLearner;
}

void MetaLearner::learnReference(std::string fileName) {

	this->dim = arffTools.getDimension(fileName);
	this->instances = arffTools.getNumberInstances(fileName);
	if (verbose) {
		std::cout << "instances: " << this->instances << std::endl;
	}

	// Set grid config
	this->gridConfig.dim_ = this->dim;
	this->gridConfig.type_ = sg::base::Linear;
	this->gridConfig.level_ = static_cast<int>(this->baseLevel);

	DataVector classesVector(instances);
	arffTools.readClasses(fileName, classesVector);

	DataMatrix trainingData(instances, dim);
	arffTools.readTrainingData(fileName, trainingData);

	bool isRegression = true; // treat everything as if it were a regression, as classification is not fully supported by Learner
	//TODO kernel must be configurable here
	LearnerLeastSquaresIdentity *referenceLearner = new LearnerLeastSquaresIdentity(nullptr, isRegression, this->verbose);

	LearnerTiming timings = referenceLearner->train(trainingData, classesVector, gridConfig, solverConfig,
			solverFinalStep, adaptivityConfiguration, false, lambda);
	this->referenceTiming = timings;
	this->ExecTimesOnStepReference = referenceLearner->getRefinementExecTimes();

	//referenceLearner->dumpFunction("referenceGridFunction", 50);
	if (this->referenceLearner != nullptr) {
		delete this->referenceLearner;
	}
	this->referenceLearner = referenceLearner;
}

//learn and test against test dataset and measure hits/mse
void MetaLearner::learnAndTest(sg::base::OperationMultipleEval *kernel, std::string datasetFileName,
		std::string testFileName, bool isBinaryClassification) {

	//always to this first
	this->learn(kernel, datasetFileName);

	size_t testDim = arffTools.getDimension(testFileName);
	size_t testInstances = arffTools.getNumberInstances(testFileName);

	DataVector testClassesVector(testInstances);
	DataMatrix testTrainingData(testInstances, testDim);

	arffTools.readClasses(testFileName, testClassesVector);
	arffTools.readTrainingData(testFileName, testTrainingData);

	if (verbose && testDim != dim) {
		std::cout << "dim of test dataset and training dataset doesn't match" << std::endl;
	}

	if (verbose) {
		std::cout << "computing classes of test dataset" << std::endl;
	}
	DataVector computedClasses = this->myLearner->predict(testTrainingData);
	if (verbose) {
		std::cout << "classes computed" << std::endl;
	}

	if (!isBinaryClassification) {
		double mse = 0.0;
		for (size_t i = 0; i < computedClasses.getSize(); i++) {
			double diff = testClassesVector.get(i) - computedClasses.get(i);
			mse += diff * diff;
		}
		mse = mse / (double) testInstances;

		if (verbose) {
			std::cout << "mse: " << mse << std::endl;
		}
	} else {
		int hits = 0;
		for (size_t i = 0; i < computedClasses.getSize(); i++) {
			int classToken = 0;

			if (computedClasses.get(i) > 0) {
				classToken = 1;
			} else {
				classToken = -1;
			}
			if (classToken == testClassesVector.get(i)) {
				hits += 1;
			}
		}
		if (verbose) {
			std::cout << "hits (%): " << ((double) hits / (double) testInstances) << std::endl;
		}
	}
}

//learn and test against the streaming implemenation
double MetaLearner::learnAndCompare(sg::base::OperationMultipleEval *kernel, std::string datasetFileName,
		size_t gridGranularity, double tolerance) {
	//always do this first
	this->learn(kernel, datasetFileName);
	this->learnReference(datasetFileName);

	DataMatrix testTrainingData(0, this->dim);

	double increment = 1.0 / static_cast<double>(gridGranularity);

	//cout << std::setprecision(30);
	size_t index = 0;
	DataVector testPoint(dim);
	for (size_t i = 0; i < dim; i++) {
		testPoint[i] = 0.0 + increment;
	}

	/*cout << "tp: ";
	 for (size_t i = 0; i < dim; i++) {
	 if (i > 0) {
	 cout << ", ";
	 }
	 cout << testPoint[i];
	 }
	 cout << endl;*/

	testTrainingData.appendRow(testPoint);

	unsigned int testInstanceCounter = 0;
	while (index < dim) {
		// make 1.0 impossible
		double prior = testPoint[index] + increment;
		if (prior < 1.0) {

			testPoint[index] += increment;
			for (size_t i = 0; i < index; i++) {
				testPoint[i] = 0.0 + increment; //skip border
			}

			/*cout << "tp: ";
			 for (size_t i = 0; i < dim; i++) {
			 if (i > 0) {
			 cout << ", ";
			 }
			 cout << testPoint[i];
			 }
			 cout << endl;
			 */
			testTrainingData.appendRow(testPoint);
			testInstanceCounter += 1;
			if (verbose && testInstanceCounter % 1000000 == 0) {
				std::cout << "testInstanceCounter (still generating): " << testInstanceCounter << std::endl;
			}
			index = 0;
		} else {
			index += 1;
		}
	}

	if (verbose) {
		std::cout << "testInstanceCounter: " << testInstanceCounter << std::endl;
		std::cout << "predicting..." << std::endl;
	}

	DataVector computedClasses = this->myLearner->predict(testTrainingData);
	if (verbose) {
		std::cout << "predicting... (reference)" << std::endl;
	}
	DataVector referenceClasses = this->referenceLearner->predict(testTrainingData);

	/*bool ok = true;
	 for (size_t i = 0; i < computedClasses.getSize(); i++) {
	 if (abs(computedClasses.get(i) - referenceClasses.get(i)) > tolerance) {
	 cout << "computed: " << computedClasses.get(i) << " but reference is: " << referenceClasses.get(i) << endl;
	 ok = false;
	 //break;
	 }
	 }

	 if (ok) {
	 cout << "test passed" << endl;
	 } else {
	 cout << "test failed" << endl;
	 }*/

	double squareSum = 0.0;

	for (size_t i = 0; i < computedClasses.getSize(); i++) {
		//TODO use fabs instead of abs?
		double temp = abs(static_cast<int>(computedClasses.get(i) - referenceClasses.get(i)));
		temp *= temp;
		squareSum += temp;

		//TODO use fabs instead of abs?
		if (verbose && abs(static_cast<int>(computedClasses.get(i) - referenceClasses.get(i))) > tolerance) {
			std::cout << "computed: " << computedClasses.get(i) << " but reference is: " << referenceClasses.get(i)
					<< std::endl;
		}
	}
	squareSum = sqrt(squareSum);

	if (verbose) {
		std::cout << "sqrt: " << squareSum << std::endl;
	}

	return squareSum;
}

void MetaLearner::writeRefinementResults(std::string fileName, std::string fileHeader,
		std::vector<std::pair<std::string, std::vector<std::pair<size_t, double> > > > datasetDetails,
		std::vector<std::pair<std::string, std::vector<std::pair<size_t, double> > > > datasetDetailsReference,
		bool referenceComparison) {
	std::ofstream myFile;
	myFile.open(fileName, std::ios::out);
	myFile << "# " << fileHeader << std::endl;
	myFile << "refinement";
	for (auto datasetDetail : datasetDetails) {
		myFile << csvSep << datasetDetail.first << "/mine";
		if (referenceComparison) {
			myFile << csvSep << datasetDetail.first << "/reference";
			myFile << csvSep << datasetDetail.first << "/speedup";
		}
	}

	myFile << std::endl;
	size_t refinements = datasetDetails[0].second.size();
	for (size_t i = 0; i < refinements; i++) {
		bool first = true;
		//for (auto datasetDetail: datasetDetails) {
		for (size_t datasetIndex = 0; datasetIndex < datasetDetails.size(); datasetIndex++) {
			//pair<size_t, double> stepTuple = datasetDetail.second[i];
			std::pair<size_t, double> stepTuple = datasetDetails[datasetIndex].second[i];
			double executionTime = stepTuple.second;
			if (first) {
				// refinement steps
				myFile << stepTuple.first;
				first = false;
			}
			// execution time at refinement step
			myFile << csvSep << executionTime;
			if (referenceComparison) {
				std::pair<std::string, std::vector<std::pair<size_t, double> > > datasetTuple =
						datasetDetailsReference[datasetIndex];
				std::pair<size_t, double> stepTupleReference = datasetTuple.second[i];
				double executionTimeReference = stepTupleReference.second;
				myFile << csvSep << executionTimeReference;
				// speedup
				myFile << csvSep << executionTimeReference / executionTime;
			}
		}
		myFile << std::endl;
	}

	myFile.close();
}

void MetaLearner::refinementAndOverallPerformance(std::vector<sg::base::OperationMultipleEval *> kernels,
		std::vector<std::string> datasets, std::vector<std::string> experimentHeaders, std::string metaInformation,
		std::string experimentName, bool referenceComparison) {
	for (sg::base::OperationMultipleEval *kernel : kernels) {
		std::vector<std::pair<std::string, std::vector<std::pair<size_t, double> > > > refinementDetails;
		std::vector<std::pair<std::string, std::vector<std::pair<size_t, double> > > > refinementDetailsReference;

		std::string kernelMetaInformation = metaInformation + " kernel: " + kernel->getImplementationName();

		std::ofstream myFile;
		std::string experimentFile = "results/data/overall_" + experimentName + "_" + kernel->getImplementationName()
				+ ".csv";
		remove(experimentFile.c_str());
		std::string experimentDetailsFile = "results/data/refinement_" + experimentName + "_"
				+ kernel->getImplementationName() + ".csv";
		myFile.open(experimentFile, std::ios::out | std::ios::app);
		myFile << "# " << kernelMetaInformation << std::endl;
		myFile << "experiment" << csvSep << "duration" << csvSep << "durationLastIt";
		if (referenceComparison) {
			myFile << csvSep << "durationRef";
			myFile << csvSep << "durationLastItRef";
			myFile << csvSep << "speedup";
			myFile << csvSep << "speedupLastIt";
		}
		myFile << std::endl;

		for (size_t i = 0; i < datasets.size(); i++) {
			this->learn(kernel, datasets[i]);
			refinementDetails.push_back(make_pair(datasets[i], this->ExecTimesOnStep));

			if (referenceComparison) {
				this->learnReference(datasets[i]);
				refinementDetailsReference.push_back(make_pair(datasets[i], this->ExecTimesOnStepReference));
			}

			myFile << experimentHeaders[i] << csvSep << this->myTiming.timeComplete_;

			size_t lastStep = this->ExecTimesOnStep.size() - 1;
			myFile << csvSep << this->ExecTimesOnStep[lastStep].second;
			if (referenceComparison) {
				myFile << csvSep << this->referenceTiming.timeComplete_;
				myFile << csvSep << this->ExecTimesOnStepReference[lastStep].second;
				myFile << csvSep << this->referenceTiming.timeComplete_ / this->myTiming.timeComplete_;
				myFile << csvSep
						<< this->ExecTimesOnStepReference[lastStep].second / this->ExecTimesOnStep[lastStep].second;
			}
			myFile << std::endl;
		}
		myFile.close();

		this->writeRefinementResults(experimentDetailsFile, kernelMetaInformation, refinementDetails,
				refinementDetailsReference, referenceComparison);
	}
}

void MetaLearner::regularGridSpeedup(sg::base::OperationMultipleEval *kernel, std::vector<size_t> dimList,
		std::vector<size_t> levelList, size_t instances, std::string metaInformation, std::string experimentName) {

	std::string kernelMetaInformation = metaInformation + " kernel: " + kernel->getImplementationName();

	std::ofstream myFile;
	std::string experimentFile = "results/data/" + experimentName + "_" + kernel->getImplementationName() + ".csv";
	remove(experimentFile.c_str());

	myFile.open(experimentFile, std::ios::out | std::ios::app);
	myFile << "# " << kernelMetaInformation << std::endl;

	myFile << "dim" << csvSep;
	myFile << "level" << csvSep;
	myFile << "duration" << csvSep;
	myFile << "durationRef" << csvSep;
	myFile << "speedup" << std::endl;

	for (size_t dim : dimList) {
		for (size_t level : levelList) {
			double duration;
			double durationReference;
			this->testRegular(kernel, dim, level, instances, duration, durationReference);

			myFile << dim << csvSep << level << csvSep << duration << csvSep << durationReference << csvSep;
			myFile << (durationReference / duration) << std::endl;
		}
		myFile << std::endl;
	}
	myFile.close();

}

void MetaLearner::appendToPerformanceRun(std::string fileName, std::string changingRowName, std::string currentValues,
		std::vector<sg::base::OperationMultipleEval *> kernels, std::vector<std::string> datasets,
		std::vector<std::string> datasetNames, std::string metaInformation, bool removeOld) {

	if (removeOld) {
		std::ifstream fileForTest(fileName.c_str());
		if (!fileForTest.good()) {
			fileForTest.close();
			std::cout << "info: no old file to delete" << std::endl;
		} else {
			fileForTest.close();
			int returnCode = remove(fileName.c_str());
			if (returnCode == 0) {
				std::cout << "old result file deleted" << std::endl;
			} else {
				std::cout << "error: should remove old file, but couldn't, return code: " << returnCode << std::endl;
				throw;
			}
		}

		std::ofstream myFile;
		myFile.open(fileName, std::ios::out | std::ios::app);
		myFile << "# " << metaInformation << std::endl;
		myFile << changingRowName;
		for (sg::base::OperationMultipleEval *kernel : kernels) {
			for (std::string datasetName : datasetNames) {
				myFile << csvSep << kernel->getImplementationName() << "/" << datasetName;
			}
		}
		myFile << std::endl;
		myFile.close();
	}

	std::ofstream myFile;
	myFile.open(fileName, std::ios::out | std::ios::app);
	std::cout << "filename: " << fileName << std::endl;
	myFile << currentValues;
	for (sg::base::OperationMultipleEval *kernel : kernels) {
		for (std::string dataset : datasets) {
			this->learn(kernel, dataset);
			myFile << csvSep << this->myTiming.timeComplete_;
		}
	}
	myFile << std::endl;
	myFile.close();
}

double fRand(double fMin, double fMax) {
	double f = (double) rand() / RAND_MAX;
	return fMin + f * (fMax - fMin);
}

void MetaLearner::testRegular(sg::base::OperationMultipleEval *kernel, size_t dim, size_t level, size_t instances,
		double &duration, double &durationReference) {

	srand(static_cast<unsigned int>(time(NULL)));
	this->dim = dim;

	// Set grid config
	this->gridConfig.dim_ = dim;
	this->gridConfig.type_ = sg::base::Linear;
	this->gridConfig.level_ = static_cast<int>(level);

	DataMatrix testTrainingData(instances, dim);

	// std::uniform_real_distribution<double> unif(0.0, 1.0);
	// std::default_random_engine re;

	// for (size_t i = 0; i < instances; i++) {
	//   for (size_t d = 0; d < dim; d++) {
	//     testTrainingData.set(i, d, unif(re));
	//   }
	// }

	//std::uniform_real_distribution<double> unif(0.0, 1.0);
	//std::default_random_engine re;

	for (size_t i = 0; i < instances; i++) {
		for (size_t d = 0; d < dim; d++) {
			//testTrainingData.set(i, d, unif(re));
			testTrainingData.set(i, d, fRand(0.0, 1.0));
		}
	}

	bool isRegression = true;
	LearnerVectorizedSubspaces *learner = new LearnerVectorizedSubspaces(kernel, isRegression, this->verbose);

	//TODO plug default kernel here instead of "kernel"
	LearnerLeastSquaresIdentity *learnerReference = new LearnerLeastSquaresIdentity(kernel, isRegression,
			this->verbose);

	duration = learner->testRegular(this->gridConfig, testTrainingData);
	durationReference = learnerReference->testRegular(this->gridConfig, testTrainingData);
}

}
}

