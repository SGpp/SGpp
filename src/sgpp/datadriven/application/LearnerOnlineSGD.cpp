#include "LearnerOnlineSGD.hpp"

namespace sg {

namespace datadriven {

LearnerOnlineSGD::LearnerOnlineSGD(
		sg::datadriven::LearnerRegularizationType& regularization,
		const bool isRegression, const bool isVerbose) :
		Learner(regularization, isRegression, isVerbose), mainTrainDataset(
		NULL), mainClasses(NULL), testTrainDataset(NULL), testClasses(
		NULL), minibatchTrainDataset(NULL), minibatchClasses(NULL), errorOnMinibatch(
		NULL), SGDCurrentIndex(0), lambda(0), gamma(0) {
}

void LearnerOnlineSGD::train(sg::base::DataMatrix& mainTrainDataset_,
		sg::base::DataVector& mainClasses_,

		sg::base::DataMatrix& testTrainDataset_,
		sg::base::DataVector& testClasses_,

		sg::base::RegularGridConfiguration& GridConfig,
		sg::datadriven::LearnerOnlineSGDRefinementConfiguration& RefineConfig,
		sg::base::AbstractRefinement& refinement,

		size_t batchSize_, double lambda_, double gamma_,

		int numRuns, std::string errorType_, std::string experimentDir) {
	using namespace sg::base;

	/*
	 * Initialization
	 */

	const int CG_IMAX = 10000;
	const double CG_EPS = 0.001;
	const double SMOOTHED_ERROR_DECLINE = 0.0001;

	if (alpha_ != NULL)
		delete alpha_;

	if (grid_ != NULL)
		delete grid_;

	if (isTrained_ == true)
		isTrained_ = false;

	InitializeGrid(GridConfig);
	if (grid_ == NULL)
		return;

	lambda = lambda_;
	gamma = gamma_;
	errorType = errorType_;

	mainTrainDataset = &mainTrainDataset_;
	mainClasses = &mainClasses_;

	testTrainDataset = &testTrainDataset_;
	testClasses = &testClasses_;

	size_t numMainData = mainTrainDataset->getNrows();
	size_t numMainDim = mainTrainDataset->getNcols();

	/*
	 * Initialize File descriptors
	 */

	std::fstream ferr0, ferr1, ferr2, fgrid, fcoor;

	ferr0.open((experimentDir + std::string("/ferr0")).c_str(),
			std::ios::out | std::ios::trunc);
	ferr1.open((experimentDir + std::string("/ferr1")).c_str(),
			std::ios::out | std::ios::trunc);
	ferr2.open((experimentDir + std::string("/ferr2")).c_str(),
			std::ios::out | std::ios::trunc);
	fgrid.open((experimentDir + std::string("/fgrid")).c_str(),
			std::ios::out | std::ios::trunc);
	fcoor.open((experimentDir + std::string("/fcoor")).c_str(),
			std::ios::out | std::ios::trunc);

	/*
	 * Initialize SGD Index order
	 */

	for (size_t i = 0; i < numMainData; i++) {
		SGDIndexOrder.push_back(i);
	}

	std::random_shuffle(SGDIndexOrder.begin(), SGDIndexOrder.end());
	SGDCurrentIndex = 0;

	/*
	 * Initialize RefinementFunctor
	 */

	RefinementFunctor *functor;

	if (RefineConfig.refinementType == "SURPLUS") {
		functor = new SurplusRefinementFunctor(alpha_,
				RefineConfig.refinementNumPoints, 0.0);

	} else if (RefineConfig.refinementType == "WEIGHTED_ERROR_MINIBATCH") {

		functor = new WeightedErrorRefinementFunctor(alpha_, grid_,
				RefineConfig.refinementNumPoints, 0.0);
		WeightedErrorRefinementFunctor* wfunctor =
				(WeightedErrorRefinementFunctor*) functor;
		wfunctor->setTrainDataset(minibatchTrainDataset);
		wfunctor->setClasses(minibatchClasses);

	} else if (RefineConfig.refinementType == "WEIGHTED_ERROR_ALL") {

		functor = new WeightedErrorRefinementFunctor(alpha_, grid_,
				RefineConfig.refinementNumPoints, 0.0);
		WeightedErrorRefinementFunctor* wfunctor =
				(WeightedErrorRefinementFunctor*) functor;
		wfunctor->setTrainDataset(mainTrainDataset);
		wfunctor->setClasses(mainClasses);

	} else if (RefineConfig.refinementType == "PERSISTENT_ERROR") {

		functor = new PersistentErrorRefinementFunctor(alpha_, grid_,
				RefineConfig.refinementNumPoints, 0.0);
		PersistentErrorRefinementFunctor* pfunctor =
				(PersistentErrorRefinementFunctor*) functor;
		pfunctor->setTrainDataset(minibatchTrainDataset);
		pfunctor->setClasses(minibatchClasses);

	} else if (RefineConfig.refinementType == "CLASSIFICATION") {
		functor = new ClassificationRefinementFunctor(alpha_, grid_,
				RefineConfig.refinementNumPoints, 0.0);
		ClassificationRefinementFunctor* cfunctor =
				(ClassificationRefinementFunctor*) functor;
		cfunctor->setTrainDataset(mainTrainDataset);
		cfunctor->setClasses(mainClasses);

	} else {
		throw base::application_exception("Invalid refinement type");
	}

	/*
	 * Perform complete run(s) through the entire training dataset
	 */

	for (int countRun = 0; countRun < numRuns; countRun++) {

		std::cout << "Run: " << countRun + 1 << std::endl;

		/*
		 * Reset Minibatch
		 */

		minibatchTrainDataset = new sg::base::DataMatrix(batchSize_,
				numMainDim);
		minibatchClasses = new sg::base::DataVector(batchSize_);
		errorOnMinibatch = new std::list<double>;

		if (RefineConfig.refinementType == "WEIGHTED_ERROR_MINIBATCH") {
			WeightedErrorRefinementFunctor* wfunctor =
					(WeightedErrorRefinementFunctor*) functor;
			wfunctor->setTrainDataset(minibatchTrainDataset);
			wfunctor->setClasses(minibatchClasses);

		} else if (RefineConfig.refinementType == "PERSISTENT_ERROR") {
			PersistentErrorRefinementFunctor* pfunctor =
					(PersistentErrorRefinementFunctor*) functor;
			pfunctor->setTrainDataset(minibatchTrainDataset);
			pfunctor->setClasses(minibatchClasses);

		}

		size_t countIterations = 0;
		while (countIterations < numMainData) {

			/*
			 * Perform SGD (until the refinement condition is true)
			 */

			/*
			 * Refinement condition: Fixed number of iterations
			 */

			if (RefineConfig.refinementCondition == "FIXED_NUMBER") {
				for (size_t countIteration = 0;
						countIteration < RefineConfig.numIterations;
						countIteration++) {
					performSGDStep();
				}
				countIterations += RefineConfig.numIterations;
			}

			/*
			 * Refinement condition: Smooth error decline
			 */

			else if (RefineConfig.refinementCondition
					== "SMOOTHED_ERROR_DECLINE") {

				if (RefineConfig.numMinibatchError == 0) {
					throw base::application_exception(
							"numMinibatchError must be > 0");
				}

				double oldErrorSum = 0;
				double oldErrorLast = 0;

				double currentError = 0;
				double ratio = 1;
				do {
					performSGDStep();
					countIterations++;

					currentError = getMinibatchError();

					if (errorOnMinibatch->size()
							>= RefineConfig.numMinibatchError) {

						// Calculate average of old minibatch errors
						for (std::list<double>::iterator it =
								errorOnMinibatch->begin();
								it != errorOnMinibatch->end(); ++it) {
							oldErrorSum += *it;
						}

						oldErrorLast = errorOnMinibatch->back();

						// Update errorOnMinibatch
						errorOnMinibatch->pop_back();
						errorOnMinibatch->push_front(currentError);

						// Update ratio
						ratio = (oldErrorLast - currentError) / oldErrorSum;
						ratio *= 1.0 / (double) RefineConfig.numMinibatchError;

						//std::cout << "Ratio: " << ratio << std::endl;

					} else {
						errorOnMinibatch->push_front(currentError);
					}
				} while (ratio > SMOOTHED_ERROR_DECLINE);
			} else {
				throw base::application_exception(
						"Invalid refinement condition");
			}

			int countTotalIterations = (int) (countIterations
					+ countRun * numMainData);

			/*
			 * Output:
			 * ferr0: Error on the current minibatch
			 * ferr1: Error on the complete training dataset
			 * ferr2: Error on the test dataset
			 * fcoor: Current coordinates of the minibatch
			 */

			ferr0 << countTotalIterations << "," << getMinibatchError()
					<< std::endl;
			ferr1 << countTotalIterations << "," << getMainError() << std::endl;
			ferr2 << countTotalIterations << "," << getTestError() << std::endl;
			fcoor << countTotalIterations << ","
					<< (*minibatchTrainDataset).toString() << std::endl;

			double percent = ((double) countTotalIterations / ((double) numMainData
					* numRuns)) * 100;
			if (percent > 100) {
				percent = 100;
			}
			printf("Accuracy: %2.2f%% (at %2.2f%%)\n", getMainError() * 100,
					percent);
			fflush(stdout);

			/*
			 * Perform one refinement operation
			 */

			// Perform the refinement
			refinement.free_refine(grid_->getStorage(), functor);
			alpha_->resizeZero(grid_->getSize());

			/*
			 * Output (fgrid): Serialized grid structure
			 */

			std::string grid_str;
			grid_->serialize(grid_str);
			fgrid << grid_str << std::endl;
		}
	}

	/*
	 * Perform CG
	 */

	std::cout << "Error before CG: " << getMainError() << std::endl;

	sg::solver::ConjugateGradients *cg = new sg::solver::ConjugateGradients(
			CG_IMAX, CG_EPS);

	sg::base::OperationMatrix *C_ = sg::op_factory::createOperationIdentity(
			*this->grid_);
	sg::datadriven::DMSystemMatrix matrix(*grid_, *mainTrainDataset, *C_,
			lambda);

	sg::base::DataVector b(alpha_->getSize());
	matrix.generateb(*mainClasses, b);

	cg->solve(matrix, *alpha_, b, true, false);

	std::cout << "Error after CG: " << getMainError() << std::endl;

	std::cout << "Error on Test Data: " << getTestError() << std::endl;

	isTrained_ = true;

	/*
	 * Close output files
	 */
	ferr0.close();
	ferr1.close();
	ferr2.close();
	fgrid.close();
	fcoor.close();
}

void LearnerOnlineSGD::performSGDStep() {
	using namespace sg::base;

	size_t numCoeff = grid_->getStorage()->size();
	size_t dim = mainTrainDataset->getNcols();

	// Get x and y pair
	DataVector x(dim);
	mainTrainDataset->getRow(SGDIndexOrder[SGDCurrentIndex], x);
	double y = mainClasses->get(SGDIndexOrder[SGDCurrentIndex]);

	// Save in minibatch
	pushMinibatch(x, y);

	// Update SGDCurrentIndex
	if (SGDCurrentIndex == SGDIndexOrder.size() - 1) {
		std::random_shuffle(SGDIndexOrder.begin(), SGDIndexOrder.end());
		SGDCurrentIndex = 0;
	} else {
		SGDCurrentIndex++;
	}

	// Calculate delta^n according to [Maier BA, 5.10]:

	// calculate error (tmp1)
	// tmp1 = (b_k^T * a^n - y_k) where
	// b_k = (phi_1(x_k) ... phi_N(x_k))
	double tmp1 = grid_->eval(*alpha_, x) - y;

	// delta^n = 2 * gamma * (b_k * tmp1 + lambda * a^n)
	DataVector delta(numCoeff);

	for (size_t i = 0; i < numCoeff; i++) {
		DataVector unit_alpha(numCoeff);
		unit_alpha.setAll(0.0);
		unit_alpha[i] = 1;

		delta[i] = grid_->eval(unit_alpha, x) * tmp1;
		delta[i] += lambda * (*alpha_)[i];
		delta[i] *= 2 * gamma;
	}

	// update alpha
	// a^{n+1} = a^n - delta^n
	for (size_t i = 0; i < numCoeff; i++) {
		(*alpha_)[i] = (*alpha_)[i] - delta[i];
	}
}

void LearnerOnlineSGD::pushMinibatch(sg::base::DataVector& x, double y) {

	size_t numMinibatchRows = minibatchTrainDataset->getNrows();
	size_t numMinibatchCols = minibatchTrainDataset->getNcols();

	sg::base::DataMatrix* newMinibatchTrainDataset = new sg::base::DataMatrix(
			numMinibatchRows, numMinibatchCols);
	sg::base::DataVector* newMinibatchClasses = new sg::base::DataVector(
			numMinibatchRows);

	newMinibatchTrainDataset->setRow(0, x);
	newMinibatchClasses->set(0, y);

	sg::base::DataVector* tmp = new sg::base::DataVector(numMinibatchCols);

	for (size_t i = 1; i < numMinibatchRows; i++) {
		minibatchTrainDataset->getRow(i - 1, *tmp);
		newMinibatchTrainDataset->setRow(i, *tmp);
		newMinibatchClasses->set(i, minibatchClasses->get(i - 1));
	}

	//delete minibatchTrainDataset;
	//delete minibatchClasses;

	minibatchTrainDataset = newMinibatchTrainDataset;
	minibatchClasses = newMinibatchClasses;
}

double LearnerOnlineSGD::getError(sg::base::DataMatrix* trainDataset,
		sg::base::DataVector* classes) {
	using namespace sg::base;

	size_t numData = trainDataset->getNrows();

	// Error vector
	DataVector *error = new DataVector(numData);
	error->setAll(0.0);

	// Values of the interpolant
	DataVector result(numData);
	sg::op_factory::createOperationMultipleEval(*grid_, trainDataset)->mult(
			*alpha_, result);

	if (errorType == "MSE") {

		for (unsigned int i = 0; i < numData; i++)
			error->set(i, result.get(i) - classes->get(i));

		// Error
		double sum = 0;
		for (unsigned int i = 0; i < numData; i++)
			sum += pow(fabs(error->get(i)), 2);

		return (sum / (double) numData);

	} else if (errorType == "ACCURACY") {

		int correct = 0;
		for (unsigned int i = 0; i < numData; i++) {
			if ((result.get(i) < 0) == (classes->get(i) < 0)) {
				correct++;
			}
		}

		return correct / (double) numData;

	} else {
		throw base::application_exception("Invalid error type");
	}

}

double LearnerOnlineSGD::getMainError() {
	return getError(mainTrainDataset, mainClasses);
}

double LearnerOnlineSGD::getTestError() {
	return getError(testTrainDataset, testClasses);
}

double LearnerOnlineSGD::getMinibatchError() {
	return getError(minibatchTrainDataset, minibatchClasses);
}

sg::base::DataVector* LearnerOnlineSGD::getAlpha() {
	if (alpha_ == NULL)
		throw base::application_exception("Not initialized");
	return alpha_;
}

LearnerOnlineSGD::~LearnerOnlineSGD() {
}
}
}
