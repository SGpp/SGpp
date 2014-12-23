#ifndef LEARNERVECTORIZEDSUBSPACES_HPP
#define LEARNERVECTORIZEDSUBSPACES_HPP

#include <vector>

#include "../operation/OperationMultipleEvalSubspace/CommonParameters.hpp"
#include "datadriven/application/LearnerBase.hpp"
#include "base/datatypes/DataMatrix.hpp"
#include "base/datatypes/DataVector.hpp"

namespace sg {
namespace datadriven {

/**
 * This class implements standard sparse grid regression
 * with an Identity matrix as regularization operator.
 *
 * Furthermore this Learner provides support for several
 * vectorization approaches covering GPUs, CPUs and coprocessors.
 */
class LearnerVectorizedSubspaces: public sg::datadriven::LearnerBase {
protected:
	//TODO: remove as soon as PseudoHashMap can extract this data
	size_t maxLevel;
	size_t dim;

	// vector<pair<size_t, double> > GFlopsOnStep;
	// vector<pair<size_t, double> > GBytesOnStep;
	std::vector<std::pair<size_t, double> > ExecTimeOnStep;

	virtual sg::datadriven::DMSystemMatrixBase* createDMSystem(sg::base::DataMatrix& trainDataset, double lambda);

	virtual void postProcessing(const sg::base::DataMatrix& trainDataset, const sg::solver::SLESolverType& solver,
			const size_t numNeededIterations);

	void adjustTrainingData(sg::base::DataMatrix &trainingData);

public:
	/**
	 * Constructor
	 *
	 * @param vecType selection of vectorization to employ
	 * @param isRegression set to true if a regression task should be executed
	 * @param isVerbose set to true in order to allow console output
	 */
	LearnerVectorizedSubspaces(sg::base::OperationMultipleEval *kernel, const bool isRegression,
			const bool isVerbose = true);

	//TODO allow to replace OperationMultipleEval

	/**
	 * Destructor
	 */
	virtual ~LearnerVectorizedSubspaces();

	virtual sg::base::DataVector predict(sg::base::DataMatrix& testDataset);

	double testRegular(const sg::base::RegularGridConfiguration& GridConfig, sg::base::DataMatrix& testDataset);

	std::vector<std::pair<size_t, double> > getRefinementExecTimes();
};

}

}

#endif /* LEARNERVECTORIZEDSUBSPACES_HPP */
