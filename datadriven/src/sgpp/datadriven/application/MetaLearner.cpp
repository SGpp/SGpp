// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <sgpp/datadriven/application/MetaLearner.hpp>

#include <random>

#include <fstream>
#include <iomanip>
#include <cstdio>

#include <sgpp/datadriven/application/LearnerLeastSquaresIdentity.hpp>
#include <sgpp/datadriven/operation/hash/simple/DatadrivenOperationCommon.hpp>
#include <sgpp/datadriven/tools/ARFFTools.hpp>
#include <sgpp/datadriven/tools/Dataset.hpp>

using namespace SGPP::base;

#include <sgpp/globaldef.hpp>

namespace SGPP {
namespace datadriven {

MetaLearner::MetaLearner(SGPP::base::RegularGridConfiguration gridConfig,
SGPP::solver::SLESolverConfiguration solverConfig,
SGPP::solver::SLESolverConfiguration solverFinalStep, SGPP::base::AdpativityConfiguration adaptivityConfiguration,
        float_t lambda, bool verbose) {
    this->csvSep = "& ";
    this->gridConfig = gridConfig;
    this->solverConfig = solverConfig;
    this->solverFinalStep = solverFinalStep;
    this->adaptivityConfiguration = adaptivityConfiguration;
    this->lambda = lambda;
    this->verbose = verbose;
    this->instances = 0;
}

void MetaLearner::learn(SGPP::datadriven::OperationMultipleEvalConfiguration& operationConfiguration,
        std::string &datasetFileName, bool isRegression) {
    std::ifstream t(datasetFileName);
    if (!t.is_open()) {
        throw;
    }
    std::stringstream buffer;
    buffer << t.rdbuf();
    std::string bufferString = buffer.str();
    this->learnString(operationConfiguration, bufferString);
}

void MetaLearner::learnString(SGPP::datadriven::OperationMultipleEvalConfiguration& operationConfiguration,
        std::string &content, bool isRegression) {

    Dataset dataset = ARFFTools::readARFFFromString(content);

    this->gridConfig.dim_ = dataset.getDimension();
    this->instances = dataset.getNumberInstances();

    if (verbose) {
        std::cout << "instances: " << this->instances << std::endl;
    }

    DataVector &classesVector = dataset.getClasses();
    DataMatrix &trainingData = dataset.getTrainingData();

//    bool isRegression = true; // treat everything as if it were a regression, as classification is not fully supported by Learner
    LearnerLeastSquaresIdentity* myLearner = new LearnerLeastSquaresIdentity(isRegression, this->verbose);
    myLearner->setImplementation(operationConfiguration);

    LearnerTiming timings = myLearner->train(trainingData, classesVector, this->gridConfig, this->solverConfig,
            this->solverFinalStep, this->adaptivityConfiguration, false, this->lambda);

    this->myTiming = timings;
    this->ExecTimesOnStep = myLearner->getRefinementExecTimes();

    if (this->myLearner != nullptr) {
        delete this->myLearner;
    }
    this->myLearner = myLearner;
}

Grid &MetaLearner::getLearnedGrid() {
    if (this->myLearner == nullptr) {
        throw;
    }
    return this->myLearner->getGrid();
}

void MetaLearner::learnReference(std::string &datasetFileName, bool isRegression) {
    std::ifstream t(datasetFileName);
    if (!t.is_open()) {
        throw;
    }
    std::stringstream buffer;
    buffer << t.rdbuf();
    std::string bufferString = buffer.str();
    this->learnReferenceString(bufferString);
}

void MetaLearner::learnReferenceString(std::string &content, bool isRegression) {

    Dataset dataset = ARFFTools::readARFFFromString(content);
    this->gridConfig.dim_ = dataset.getDimension();
    this->instances = dataset.getNumberInstances();

    if (verbose) {
        std::cout << "instances: " << this->instances << std::endl;
    }

    DataVector &classesVector = dataset.getClasses();
    DataMatrix &trainingData = dataset.getTrainingData();

//    // treat everything as if it were a regression, as classification is not fully supported by Learner
//    bool isRegression = true;
    LearnerLeastSquaresIdentity* referenceLearner = new LearnerLeastSquaresIdentity(isRegression, this->verbose);
    SGPP::datadriven::OperationMultipleEvalConfiguration operationConfiguration(OperationMultipleEvalType::DEFAULT,
            OperationMultipleEvalSubType::DEFAULT, "STREAMING");
    referenceLearner->setImplementation(operationConfiguration);

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

void MetaLearner::learnAndTest(SGPP::datadriven::OperationMultipleEvalConfiguration& operationConfiguration,
        std::string &datasetFileName, std::string &testFileName, bool isRegression) {
    std::ifstream dataFile(datasetFileName);
    std::stringstream bufferData;
    bufferData << dataFile.rdbuf();
    std::ifstream testFile(datasetFileName);
    std::stringstream bufferTest;
    bufferTest << testFile.rdbuf();
    std::string bufferDataString = bufferData.str();
    std::string bufferTestString = bufferTest.str();
    this->learnAndTestString(operationConfiguration, bufferDataString, bufferTestString, isRegression);
}

//learn and test against test dataset and measure hits/mse
void MetaLearner::learnAndTestString(SGPP::datadriven::OperationMultipleEvalConfiguration& operationConfiguration,
        std::string &dataContent, std::string &testContent, bool isRegression) {

    //always to this first
    this->learnString(operationConfiguration, dataContent);

    Dataset testDataset = ARFFTools::readARFFFromString(testContent);
    size_t testDim = testDataset.getDimension();
    size_t testInstances = testDataset.getNumberInstances();

    DataVector &testClassesVector = testDataset.getClasses();
    DataMatrix &testTrainingData = testDataset.getTrainingData();

    if (verbose && testDim != this->gridConfig.dim_) {
        std::cout << "dim of test dataset and training dataset doesn't match" << std::endl;
    }

    if (verbose) {
        std::cout << "computing classes of test dataset" << std::endl;
    }

    DataVector computedClasses(testTrainingData.getNrows());
    this->myLearner->predict(testTrainingData, computedClasses);

    if (verbose) {
        std::cout << "classes computed" << std::endl;
    }

    if (isRegression) {
        float_t mse = 0.0;

        for (size_t i = 0; i < computedClasses.getSize(); i++) {
            float_t diff = testClassesVector.get(i) - computedClasses.get(i);
            mse += diff * diff;
        }

        mse = mse / (float_t) testInstances;

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
            std::cout << "hits (%): " << ((float_t) hits / (float_t) testInstances) << std::endl;
        }
    }
}

//learn and test against the streaming implemenation
float_t MetaLearner::learnAndCompare(SGPP::datadriven::OperationMultipleEvalConfiguration& operationConfiguration,
        std::string &datasetFileName, size_t gridGranularity) {
    std::ifstream t(datasetFileName);
    if (!t.is_open()) {
        throw;
    }
    std::stringstream buffer;
    buffer << t.rdbuf();
    std::string bufferString = buffer.str();
    return this->learnAndCompareString(operationConfiguration, bufferString, gridGranularity);
}

//learn and test against the streaming implemenation
float_t MetaLearner::learnAndCompareString(SGPP::datadriven::OperationMultipleEvalConfiguration& operationConfiguration,
        std::string &content, size_t gridGranularity) {
    //always do this first
    this->learnString(operationConfiguration, content);
    this->learnReferenceString(content);

    DataMatrix testTrainingData(0, this->gridConfig.dim_);

    float_t increment = 1.0 / static_cast<float_t>(gridGranularity);

    size_t index = 0;
    DataVector testPoint(this->gridConfig.dim_);

    for (size_t i = 0; i < this->gridConfig.dim_; i++) {
        testPoint[i] = 0.0 + increment;
    }

    testTrainingData.appendRow(testPoint);

    unsigned int testInstanceCounter = 0;

    while (index < this->gridConfig.dim_) {
        // make 1.0 impossible
        float_t prior = testPoint[index] + increment;

        if (prior < 1.0) {

            testPoint[index] += increment;

            for (size_t i = 0; i < index; i++) {
                testPoint[i] = 0.0 + increment; //skip border
            }

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
        std::cout << "testInstanceCounter: " << (testInstanceCounter + 1) << std::endl;
        std::cout << "predicting..." << std::endl;
    }

    DataVector computedClasses(testTrainingData.getNrows());
    this->myLearner->predict(testTrainingData, computedClasses);

    if (verbose) {
        std::cout << "predicting... (reference)" << std::endl;
    }

    DataVector referenceClasses(testTrainingData.getNrows());
    this->referenceLearner->predict(testTrainingData, referenceClasses);

    float_t squareSum = 0.0;

    // SGPP::base::DataVector *myAlpha = this->myLearner->alpha_;
    // for (size_t i = 0; i < myAlpha->getSize();i++) {
    //    std::cout << "alpha[ " << i << "]=" << (*myAlpha)[i] << ", ";
    // }
    // std::cout << std::endl;

    double maxDiff = -1.0;
    double firstValue = -1.0;
    double secondValue = -1.0;
    for (size_t i = 0; i < computedClasses.getSize(); i++) {
        float_t diff = fabs(computedClasses.get(i) - referenceClasses.get(i));

        if (diff > maxDiff) {
            maxDiff = diff;
            firstValue = computedClasses.get(i);
            secondValue = referenceClasses.get(i);
        }

        diff *= diff;
        squareSum += diff;
    }

    squareSum = sqrt(squareSum);

    if (verbose) {
        std::cout << "sqrt: " << squareSum << std::endl;
        std::cout << "maxDiff: " << maxDiff << " value: " << firstValue << " second: " << secondValue << std::endl;
    }

    return squareSum;
}

void MetaLearner::writeRefinementResults(std::string fileName, std::string fileHeader,
        std::vector<std::pair<std::string, std::vector<std::pair<size_t, float_t> > > > datasetDetails,
        std::vector<std::pair<std::string, std::vector<std::pair<size_t, float_t> > > > datasetDetailsReference,
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

        for (size_t datasetIndex = 0; datasetIndex < datasetDetails.size(); datasetIndex++) {
            std::pair<size_t, float_t> stepTuple = datasetDetails[datasetIndex].second[i];
            float_t executionTime = stepTuple.second;

            if (first) {
                // refinement steps
                myFile << stepTuple.first;
                first = false;
            }

            // execution time at refinement step
            myFile << csvSep << executionTime;

            if (referenceComparison) {
                std::pair<std::string, std::vector<std::pair<size_t, float_t> > > datasetTuple =
                        datasetDetailsReference[datasetIndex];
                std::pair<size_t, float_t> stepTupleReference = datasetTuple.second[i];
                float_t executionTimeReference = stepTupleReference.second;
                myFile << csvSep << executionTimeReference;
                // speedup
                myFile << csvSep << executionTimeReference / executionTime;
            }
        }

        myFile << std::endl;
    }

    myFile.close();
}

void MetaLearner::refinementAndOverallPerformance(
        std::vector<SGPP::datadriven::OperationMultipleEvalConfiguration*> operationConfigurations,
        std::vector<std::string> datasets, std::vector<std::string> experimentHeaders, std::string metaInformation,
        std::string experimentName, bool referenceComparison) {
    for (SGPP::datadriven::OperationMultipleEvalConfiguration * operationConfiguration : operationConfigurations) {
        std::vector<std::pair<std::string, std::vector<std::pair<size_t, float_t> > > > refinementDetails;
        std::vector<std::pair<std::string, std::vector<std::pair<size_t, float_t> > > > refinementDetailsReference;

        std::string kernelMetaInformation = metaInformation + " kernel: " + operationConfiguration->getName();

        std::ofstream myFile;
        std::string experimentFile = "results/data/overall_" + experimentName + "_" + operationConfiguration->getName()
                + ".csv";
        remove(experimentFile.c_str());
        std::string experimentDetailsFile = "results/data/refinement_" + experimentName + "_"
                + operationConfiguration->getName() + ".csv";
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
            this->learn(*operationConfiguration, datasets[i]);
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

void MetaLearner::regularGridSpeedup(SGPP::datadriven::OperationMultipleEvalConfiguration& operationConfiguration,
        std::vector<size_t> dimList, std::vector<size_t> levelList, size_t instances, std::string metaInformation,
        std::string experimentName) {

    std::string kernelMetaInformation = metaInformation + " kernel: " + operationConfiguration.getName();

    std::ofstream myFile;
    std::string experimentFile = "results/data/" + experimentName + "_" + operationConfiguration.getName() + ".csv";
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
            float_t duration;
            float_t durationReference;
            this->testRegular(operationConfiguration, dim, level, instances, duration, durationReference);

            myFile << dim << csvSep << level << csvSep << duration << csvSep << durationReference << csvSep;
            myFile << (durationReference / duration) << std::endl;
        }

        myFile << std::endl;
    }

    myFile.close();

}

void MetaLearner::appendToPerformanceRun(std::string fileName, std::string changingRowName, std::string currentValues,
        std::vector<SGPP::datadriven::OperationMultipleEvalConfiguration*> operationConfigurations,
        std::vector<std::string> datasets, std::vector<std::string> datasetNames, std::string metaInformation,
        bool removeOld) {

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

        for (SGPP::datadriven::OperationMultipleEvalConfiguration * operationConfiguration : operationConfigurations) {
            for (std::string datasetName : datasetNames) {
                myFile << csvSep << operationConfiguration->getName() << "/" << datasetName;
            }
        }

        myFile << std::endl;
        myFile.close();
    }

    std::ofstream myFile;
    myFile.open(fileName, std::ios::out | std::ios::app);
    std::cout << "filename: " << fileName << std::endl;
    myFile << currentValues;

    for (SGPP::datadriven::OperationMultipleEvalConfiguration * operationConfiguration : operationConfigurations) {
        for (std::string dataset : datasets) {
            this->learn(*operationConfiguration, dataset);
            myFile << csvSep << this->myTiming.timeComplete_;
        }
    }

    myFile << std::endl;
    myFile.close();
}

float_t fRand(float_t fMin, float_t fMax) {
    float_t f = (float_t) rand() / RAND_MAX;
    return fMin + f * (fMax - fMin);
}

void MetaLearner::testRegular(SGPP::datadriven::OperationMultipleEvalConfiguration& operationConfiguration, size_t dim,
        size_t level, size_t instances, float_t& duration, float_t& durationReference) {

    srand(static_cast<unsigned int>(time(NULL)));
    this->gridConfig.dim_ = dim;

    DataMatrix testTrainingData(instances, dim);

    for (size_t i = 0; i < instances; i++) {
        for (size_t d = 0; d < dim; d++) {
            //testTrainingData.set(i, d, unif(re));
            testTrainingData.set(i, d, fRand(0.0, 1.0));
        }
    }

    bool isRegression = true;
    LearnerLeastSquaresIdentity* learner = new LearnerLeastSquaresIdentity(isRegression, this->verbose);
    learner->setImplementation(operationConfiguration);

    LearnerLeastSquaresIdentity* learnerReference = new LearnerLeastSquaresIdentity(isRegression, this->verbose);
    SGPP::datadriven::OperationMultipleEvalConfiguration referenceOperationConfiguration(
            OperationMultipleEvalType::STREAMING, OperationMultipleEvalSubType::DEFAULT, "STREAMING");
    learnerReference->setImplementation(referenceOperationConfiguration);

    duration = learner->testRegular(this->gridConfig, testTrainingData);
    durationReference = learnerReference->testRegular(this->gridConfig, testTrainingData);
}

LearnerTiming MetaLearner::getLearnerTiming() {
    return this->myTiming;
}

LearnerTiming MetaLearner::getLearnerReferenceTiming() {
    return this->referenceTiming;
}

}
}

