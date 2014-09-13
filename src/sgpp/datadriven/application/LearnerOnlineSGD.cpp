#include <limits>

#include "LearnerOnlineSGD.hpp"

namespace sg {

namespace datadriven {

LearnerOnlineSGD::LearnerOnlineSGD(
		sg::datadriven::LearnerRegularizationType& regularization,
		const bool isRegression, const bool isVerbose) :
		Learner(regularization, isRegression, isVerbose), mainTrainDataset(
		NULL), mainClasses(NULL), testTrainDataset(NULL), testClasses(
		NULL), minibatchTrainDataset(NULL), minibatchClasses(NULL), errorOnMinibatch(
		NULL),  errorMB(NULL), errorTrainData(NULL), alphaAvg_(NULL), SGDCurrentIndex(0), lambda(0), gamma(0), countIterations_(0) {
}

void LearnerOnlineSGD::train(sg::base::DataMatrix& mainTrainDataset_,
		sg::base::DataVector& mainClasses_,

		sg::base::DataMatrix& testTrainDataset_,
		sg::base::DataVector& testClasses_,

		sg::base::RegularGridConfiguration& GridConfig,
		sg::datadriven::LearnerOnlineSGDRefinementConfiguration& RefineConfig,
		sg::base::AbstractRefinement& refinement,

		size_t batchSize_, double lambda_, double gamma0,

		int numRuns, std::string errorType_, std::string experimentDir) {
	using namespace sg::base;

	/*
	 * Initialization
	 */

	const int CG_IMAX = 200;
	const double CG_EPS = 0.0001;
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
	errorType = errorType_;

	mainTrainDataset = &mainTrainDataset_;
	mainClasses = &mainClasses_;

	testTrainDataset = &testTrainDataset_;
	testClasses = &testClasses_;

	minibatchSize = batchSize_;

	if (errorMB != NULL) {
		delete errorMB;
	}
	errorMB = new DataVector(minibatchSize);
	errorMB->setAll(0.0);

	if (errorTrainData != NULL) {
		delete errorTrainData;
	}
	errorTrainData = new DataVector(mainClasses->getSize());

	numMainData_ = mainTrainDataset->getNrows();
	numMainDim_ = mainTrainDataset->getNcols();

	minibatchTrainDataset = new sg::base::DataMatrix(0, numMainDim_);
	minibatchTrainDataset->addSize(minibatchSize);
	minibatchClasses = new DataVector(minibatchSize);

	DataVector * errorPerBasisFunction = new DataVector(grid_->getSize());

	alphaAvg_ = new DataVector(grid_->getSize());


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

	for (size_t i = 0; i < numMainData_; i++) {
		SGDIndexOrder.push_back(i);
	}

	// TODO: remove seed after experiments are done
	std::random_shuffle(SGDIndexOrder.begin(), SGDIndexOrder.end());
	SGDCurrentIndex = 0;

	/*
	 * Initialize RefinementFunctor
	 */

	RefinementFunctor *functor = NULL;

	if (RefineConfig.refinementType == "SURPLUS") {
		functor = new SurplusRefinementFunctor(alpha_,
				RefineConfig.refinementNumPoints, 0.0);

	} else if (RefineConfig.refinementType == "MSE") {
	    /*functor = new SurplusRefinementFunctor(alpha_,
	        RefineConfig.refinementNumPoints, 0.0);*/
	    /*functor = new SurplusRefinementFunctor(alpha_,
	            RefineConfig.refinementNumPoints, 0.0);*/
	    functor = new SurplusRefinementFunctor(errorPerBasisFunction,
	                RefineConfig.refinementNumPoints, 0.0);
	  } else if (RefineConfig.refinementType == "WEIGHTED_ERROR_MINIBATCH") {

		functor = new WeightedErrorRefinementFunctor(alpha_, grid_,
				RefineConfig.refinementNumPoints,
				-std::numeric_limits<double>::infinity());
		WeightedErrorRefinementFunctor* wfunctor =
				(WeightedErrorRefinementFunctor*) functor;
		wfunctor->setTrainDataset(minibatchTrainDataset);
		wfunctor->setClasses(minibatchClasses);
		wfunctor->setErrors(errorMB);

	} else if (RefineConfig.refinementType == "WEIGHTED_ERROR_ALL") {
		// FIXME: this case is not accounted for
		/*functor = new WeightedErrorRefinementFunctor(alpha_, grid_,
		 RefineConfig.refinementNumPoints, 0.0);
		 WeightedErrorRefinementFunctor* wfunctor =
		 (WeightedErrorRefinementFunctor*) functor;

		 wfunctor->setTrainDataset(mainTrainDataset);
		 wfunctor->setClasses(mainClasses);*/

	} else if (RefineConfig.refinementType == "PERSISTENT_ERROR") {

		functor = new PersistentErrorRefinementFunctor(alphaAvg_, grid_,
				RefineConfig.refinementNumPoints,
				-std::numeric_limits<double>::infinity());
		PersistentErrorRefinementFunctor* pfunctor =
				(PersistentErrorRefinementFunctor*) functor;
		pfunctor->setTrainDataset(minibatchTrainDataset);
		pfunctor->setClasses(minibatchClasses);
		pfunctor->setErrors(errorMB);

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

	countIterations_ = 0;
	for (int countRun = 0; countRun < numRuns; countRun++) {

		std::cout << "Run: " << countRun + 1 << std::endl;

		/*
		 * Reset Minibatch
		 */


		
		//while (countIterations_ < numMainData_) {

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
					gamma = gamma0*pow(1+gamma0*lambda_*countIterations_, -2.0/3);
					performSGDStep();
					countIterations_++;
				}
				
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
					countIterations_++;

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

			int countTotalIterations = (int) (countIterations_
					+ countRun * numMainData_);

			/*
			 * Output:
			 * ferr0: Error on the current minibatch
			 * ferr1: Error on the complete training dataset
			 * ferr2: Error on the test dataset
			 * fcoor: Current coordinates of the minibatch
			 */

			ferr0 << countIterations_ << "," << getMinibatchError()
					<< std::endl;
			ferr1 << countIterations_ << "," << getMainError() << std::endl;
			ferr2 << countIterations_ << "," << getTestError() << std::endl;
			//fcoor << countTotalIterations << ","
			//		<< (*minibatchTrainDataset).toString() << std::endl;

			double percent = ((double) countTotalIterations
					/ ((double) numMainData_ * numRuns)) * 100;
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
			//std::cout << "Indicator values: ";
			if (RefineConfig.refinementType == "PERSISTENT_ERROR" || RefineConfig.refinementType == "WEIGHTED_ERROR_MINIBATCH") {

				double accuracy = getError(minibatchTrainDataset,
						minibatchClasses, errorMB);
				if (accuracy == 1.0) {
					continue;
				}

			}

			if (RefineConfig.refinementType == "MSE") {

					  errorTrainData->sqr();
		  OperationMultipleEval* eval = sg::op_factory::createOperationMultipleEval(*grid_,
				  mainTrainDataset);
		  eval->multTranspose(*errorTrainData, *errorPerBasisFunction);
		  delete eval;
			  refinement.free_refine(grid_->getStorage(), functor);
			}
			else {
			  refinement.free_refine(grid_->getStorage(), functor);
			}
			//std::cout << std::endl;

			alpha_->resizeZero(grid_->getSize());
			alphaAvg_->resizeZero(grid_->getSize());

			/*
			 * Output (fgrid): Serialized grid structure
			 */

			std::string grid_str;
			grid_->serialize(grid_str);
			fgrid << grid_str << std::endl;
		//}
	}

	/*
	 * Perform CG
	 */

	errorType = "ACCURACY";
	std::cout << "Error before CG (" << errorType << "): " << getMainError() << std::endl;

	errorType = "MSE";
	std::cout << "Error before CG (" << errorType << "): " << getMainError() << std::endl;

	sg::solver::ConjugateGradients *cg = new sg::solver::ConjugateGradients(
			CG_IMAX, CG_EPS);

	sg::base::OperationMatrix *C_ = sg::op_factory::createOperationIdentity(
			*this->grid_);
	sg::datadriven::DMSystemMatrix matrix(*grid_, *mainTrainDataset, *C_,
			lambda);

	sg::base::DataVector b(alpha_->getSize());
	matrix.generateb(*mainClasses, b);

	cg->solve(matrix, *alpha_, b, true, false);

	errorType = "ACCURACY";
	std::cout << "Error after CG (" << errorType << "): " << getMainError() << std::endl;

	errorType = "MSE";
	std::cout << "Error after CG (" << errorType << "): " << getMainError() << std::endl;

	errorType = "ACCURACY";
	std::cout << "Error on test Data (" << errorType << "): " << getMainError() << std::endl;

	errorType = "MSE";
	std::cout << "Error on test Data (" << errorType << "): " << getMainError() << std::endl;

	isTrained_ = true;

	/*
	 * Close output files
	 */
	ferr0.close();
	ferr1.close();
	ferr2.close();
	fgrid.close();
	fcoor.close();

	delete C_;
	delete cg;
	delete functor;
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
	DataVector delta(*alpha_);
	DataVector unit_alpha(numCoeff);
	unit_alpha.setAll(0.0);

	DataVector singleAlpha(1);
	singleAlpha[0] = 1.0;


	DataMatrix dm(x.getPointer(), 1, x.getSize());
	sg::op_factory::createOperationMultipleEval(*grid_, &dm)->multTranspose(
				singleAlpha, delta);
	alpha_->mult(1-gamma*lambda);
	alpha_->axpy(-gamma*tmp1, delta);

	//double mu = 0.1; // boring exponential smoothing

	// L. Bottou exciting smoothing
	size_t t1 = (countIterations_ > numMainDim_+1) ? countIterations_ - numMainDim_ : 1;
	size_t t2 = (countIterations_ > numMainData_+1) ? countIterations_ - numMainDim_: 1;
	double mu = (t1>t2) ? static_cast<double>(t1) : static_cast<double>(t2);
	mu = 1.0/mu;


	alphaAvg_->mult(1-mu);
	alphaAvg_->axpy(mu, *alpha_);

	/*for (size_t i = 0; i < numCoeff; i++) {


		unit_alpha[i] = 1;

		delta[i] = grid_->eval(unit_alpha, x) * tmp1;
		delta[i] += lambda * (*alpha_)[i];
		delta[i] *= 2 * gamma;
		unit_alpha[i] = 0;
	}

	// update alpha
	// a^{n+1} = a^n - delta^n
	for (size_t i = 0; i < numCoeff; i++) {
		(*alpha_)[i] = (*alpha_)[i] - delta[i];
	}*/
}

void LearnerOnlineSGD::pushMinibatch(sg::base::DataVector& x, double y) {
	static size_t next_idx = 0;
	if (minibatchTrainDataset->getUnused() > 0) {
		minibatchTrainDataset->appendRow(x);
		(*minibatchClasses)[next_idx] = y;
	} else {
		minibatchTrainDataset->setRow(next_idx, x);
		(*minibatchClasses)[next_idx] = y;
	}
	next_idx = (next_idx + 1) % minibatchSize;
}

double LearnerOnlineSGD::getError(sg::base::DataMatrix* trainDataset,
		sg::base::DataVector* classes, sg::base::DataVector *error) {
	using namespace sg::base;

	size_t numData = trainDataset->getNrows();
	bool cleanup = false;
	double res = -1.0;

	// Error vector
	if (error == NULL) {
		error = new DataVector(numData);
		cleanup = true;
	}
	error->setAll(0.0);

	// Values of the interpolant
	// FIXME result is not evaluated correctly! (only -nan)
	DataVector result(numData);
	OperationMultipleEval* eval = sg::op_factory::createOperationMultipleEval(*grid_, trainDataset);
	eval->mult(*alphaAvg_, result);
	//eval->mult(*alpha_, result);
	delete eval;

	/*
	DataVector tmp_vec(trainDataset->getNcols());
	trainDataset->getRow(0, tmp_vec);
	std::cout << tmp_vec.toString() << std::endl;
	*/

	if (errorType == "MSE") {

		for (unsigned int i = 0; i < numData; i++) {
			error->set(i, classes->get(i) - result.get(i));
		}

		// Error
		double sum = 0;
		for (unsigned int i = 0; i < numData; i++) {
			sum += error->get(i) * error->get(i);
		}

		res = (sum / (double) numData);

		if (cleanup) {
			delete error;
		}

	} else if (errorType == "ACCURACY") {
		unsigned int correct = 0;
		for (unsigned int i = 0; i < numData; i++) {
			correct += (result.get(i) < 0) == (classes->get(i) < 0) ? 1 : 0;
			error->set(i, classes->get(i) - result.get(i));

		}
		res = static_cast<double>(correct) / static_cast<double>(numData);
		if (cleanup) {
			delete error;
		}

	} else {
		if (cleanup) {
			delete error;
		}

		throw base::application_exception("Invalid error type");
	}
	return res;

}


double LearnerOnlineSGD::getMainError() {
	return getError(mainTrainDataset, mainClasses, NULL);
}

double LearnerOnlineSGD::getTestError() {
	return getError(testTrainDataset, testClasses, NULL);
}

double LearnerOnlineSGD::getMinibatchError() {
	return getError(minibatchTrainDataset, minibatchClasses, NULL);
}

sg::base::DataVector* LearnerOnlineSGD::getAlpha() {
	if (alpha_ == NULL)
		throw base::application_exception("Not initialized");
	return alpha_;
}

LearnerOnlineSGD::~LearnerOnlineSGD() {
	if (alphaAvg_ != NULL){
		delete alphaAvg_;
	}
	if (minibatchTrainDataset != NULL){
		delete minibatchTrainDataset;
	}
	if (minibatchClasses != NULL){
		delete minibatchClasses;
	}
}
}
}
