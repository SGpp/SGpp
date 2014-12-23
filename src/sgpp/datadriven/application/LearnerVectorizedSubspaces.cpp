#include <datadriven/tools/LearnerVectorizedPerformanceCalculator.hpp>
#include <base/exception/factory_exception.hpp>

#include "datadriven/algorithm/DMSystemMatrixSubspaces.hpp"
#include "LearnerVectorizedSubspaces.hpp"
//#include "datadriven/DatadrivenOpFactory.hpp"
#include "base/operation/BaseOpFactory.hpp"

//TODO where to put type/subtype here

//#include "OperatorFactory.hpp"

namespace sg {
namespace datadriven {

LearnerVectorizedSubspaces::LearnerVectorizedSubspaces(sg::base::OperationMultipleEval *kernel, const bool isRegression,
		const bool verbose) :
		sg::datadriven::LearnerBase(isRegression, verbose) {
	this->dim = dim;
	this->maxLevel = maxLevel;
}

LearnerVectorizedSubspaces::~LearnerVectorizedSubspaces() {
}

sg::datadriven::DMSystemMatrixBase* LearnerVectorizedSubspaces::createDMSystem(sg::base::DataMatrix& trainDataset,
		double lambda) {
	if (this->grid_ == NULL)
		return NULL;

	//TODO: make it clear somewhere that the training data is modified
	this->adjustTrainingData(trainDataset);

	return new DMSystemMatrixSubspaces(*(this->grid_), trainDataset, lambda);
}

//adjust training data for faster kernel index calculation
//makes sure that every data tuple contains only values that are < 1 (and not <= 1)
//by selecting the last floating point represenation before 1
void LearnerVectorizedSubspaces::adjustTrainingData(sg::base::DataMatrix &trainingData) {
	size_t dim = trainingData.getNcols();
	float lastRepresentation = 1.0;
	size_t *asInteger = (size_t *) &lastRepresentation;
	*asInteger -= 1;
	sg::base::DataVector dataTuple(dim);
	for (size_t i = 0; i < trainingData.getNrows(); i++) {
		for (size_t j = 0; j < dim; j++) {
			if (trainingData.get(i, j) >= 1.0) {
				trainingData.set(i, j, lastRepresentation);
			}
		}
	}
}

void LearnerVectorizedSubspaces::postProcessing(const sg::base::DataMatrix& trainDataset,
		const sg::solver::SLESolverType& solver, const size_t numNeededIterations) {
	LearnerVectorizedPerformance currentPerf = LearnerVectorizedPerformanceCalculator::getGFlopAndGByte(*this->grid_,
			trainDataset.getNrows(), solver, numNeededIterations, sizeof(double));

	this->GFlop_ += currentPerf.GFlop_;
	this->GByte_ += currentPerf.GByte_;

	this->stepGFlop_ = currentPerf.GFlop_;
	this->stepGByte_ = currentPerf.GByte_;

	// this->GFlopsOnStep.push_back(make_pair(this->currentRefinementStep, this->stepGFlop_ / this->stepExecTime_));
	// this->GBytesOnStep.push_back(make_pair(this->currentRefinementStep, this->stepGByte_ / this->stepExecTime_));

	this->ExecTimeOnStep.push_back(std::make_pair(this->currentRefinementStep, this->stepExecTime_));

	// Calculate GFLOPS and GBytes/s and write them to console
	if (this->isVerbose_) {
		std::cout << std::endl;
		std::cout << "Current refinement GFlop/s: " << this->stepGFlop_ / this->stepExecTime_ << std::endl;
		std::cout << "Current refinement GByte/s: " << this->stepGByte_ / this->stepExecTime_ << std::endl;
		std::cout << "Overall averaged GFlop/s: " << this->GFlop_ / this->execTime_ << std::endl;
		std::cout << "Overall averaged GByte/s: " << this->GByte_ / this->execTime_ << std::endl;
		std::cout << std::endl;
	}
}

std::vector<std::pair<size_t, double> > LearnerVectorizedSubspaces::getRefinementExecTimes() {
	return this->ExecTimeOnStep;
}

sg::base::DataVector LearnerVectorizedSubspaces::predict(sg::base::DataMatrix& testDataset) {
	sg::base::DataMatrix tmpDataSet(testDataset);
	size_t originalSize = testDataset.getNrows();
	//size_t paddedSize = PaddingAssistant::padDataset(tmpDataSet);

	//cout << "originalSize: " << originalSize << endl;
	//AbstractOperationMultipleEval* MultEval = createAbstractOperationMultipleEval(*grid_, this->kernelType_, &tmpDataSet, this->dim, this->maxLevel);
	sg::base::OperationMultipleEval *MultEval = sg::op_factory::createOperationMultipleEval(*(this->grid_), tmpDataSet);
	//TODO is this a good place to implicitly do padding?
	size_t paddedSize = tmpDataSet.getNrows();
	//cout << "paddedSize: " << paddedSize << endl;

	sg::base::DataVector classesComputed(paddedSize);

	classesComputed.setAll(0.0);

	MultEval->multTranspose(*alpha_, classesComputed);
	delete MultEval;

	// removed the padded instances
	classesComputed.resize(originalSize);

	return classesComputed;
}

double LearnerVectorizedSubspaces::testRegular(const sg::base::RegularGridConfiguration& GridConfig,
		sg::base::DataMatrix& testDataset) {
	sg::base::DataMatrix tmpDataSet(testDataset);
	size_t originalSize = testDataset.getNrows();
	//size_t paddedSize = PaddingAssistant::padDataset(tmpDataSet);

	InitializeGrid(GridConfig);

	//cout << "originalSize: " << originalSize << endl;
	sg::base::OperationMultipleEval* MultEval = sg::op_factory::createOperationMultipleEval(*(this->grid_), tmpDataSet);
	//TODO don't forget to adjust changes in the padding api here, too
	size_t paddedSize = tmpDataSet.getNrows();
	//cout << "paddedSize: " << paddedSize << endl;

	sg::base::DataVector classesComputed(paddedSize);

	classesComputed.setAll(0.0);

	execTime_ = 0.0;

	sg::base::SGppStopwatch* myStopwatch = new sg::base::SGppStopwatch();
	myStopwatch->start();

	MultEval->multTranspose(*alpha_, classesComputed);
	double stopTime = myStopwatch->stop();
	this->execTime_ += stopTime;
	std::cout << "execution duration: " << this->execTime_ << std::endl;
	delete MultEval;

	// removed the padded instances
	classesComputed.resize(originalSize);

	return stopTime;
}

}
}
