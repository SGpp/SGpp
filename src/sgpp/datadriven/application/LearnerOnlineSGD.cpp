#include "LearnerOnlineSGD.hpp"

namespace sg {

namespace datadriven {

LearnerOnlineSGD::LearnerOnlineSGD(
		sg::datadriven::LearnerRegularizationType& regularization,
		const bool isRegression, const bool isVerbose) :
		Learner(regularization, isRegression, isVerbose), batch(NULL), batchSize(
				0), batchIndex(0) {
}

void LearnerOnlineSGD::train(sg::base::DataMatrix& trainDataset,
		sg::base::DataVector& classes,
		sg::base::RegularGridConfiguration& GridConfig,

		size_t numIterations, size_t batchSize_, double lambda, double gamma,

		sg::base::AbstractRefinement& refinement, int refinementCondition,
		int refinementType, int refinementNumPoints,

		int numRuns, std::ostream **outputStreams_) {
	using namespace sg::base;

	/*
	 * Constants
	 */
	const int CG_IMAX;
	const double CG_EPS;
	const double SMOOTHED_ERROR_DECLINE = 0.1;
	/*
	 * Initialization
	 */

	if (alpha_ != NULL)
		delete alpha_;

	if (grid_ != NULL)
		delete grid_;

	if (isTrained_ == true)
		isTrained_ = false;

	InitializeGrid(GridConfig);
	if (grid_ == NULL)
		return;

	batchSize = batchSize_;
	batch = new size_t[batchSize];

	outputStreams = outputStreams_;

	sg::base::RefinementFunctor *functor;

	switch (refinementType) {
	case 0:
		functor = new sg::base::SurplusRefinementFunctor(alpha_,
				refinementNumPoints, 0.0);
		break;
	default:
		throw base::application_exception("Invalid refinement type");
		break;
	}

	/*
	 * Step 1: Perform complete run(s) through the entire training dataset (on average)
	 */

	size_t numData = trainDataset.getNrows();

	for (int countRun = 0; countRun < numRuns; countRun++) {

		size_t countTotalIterations = 0;
		while (countTotalIterations < numData) {

			/*
			 * Step 1.1: Perform SGD depending on the refinement condition
			 */

			/*
			 * Condition 1.1.1: Fixed number of iterations
			 */

			if (refinementCondition == 0) {

				for (size_t countIteration = 0; countIteration < numIterations;
						countIteration++) {
					performSGDStep(trainDataset, classes, lambda, gamma);

				}

				countTotalIterations += numIterations;
			}

			/*
			 * Condition 1.1.2: Smooth error decline
			 * Smooth function: MSE
			 */

			if (refinementCondition == 1) {

				double lastError = 1;
				double currentError;
				double ratio;
				do {
					performSGDStep(trainDataset, classes, lambda, gamma);
					countTotalIterations++;

					currentError = errorOnMinibatch(trainDataset, classes);
					ratio = (lastError - currentError) / lastError;
					lastError = currentError;
				} while (ratio < SMOOTHED_ERROR_DECLINE);
			}

			/*
			 * Output 1.1.3 (ferr0): On the current minibatch
			 */

			char ferr0[20];
			int n1 = (int) (countTotalIterations + countRun * numRuns);
			double n2 = errorOnMinibatch(trainDataset, classes);
			std::sprintf(ferr0, "%d,%f", n1, n2);
			output(0, ferr0);

			/*
			 * Output 1.1.4 (ferr1): On the complete training dataset
			 */

			char ferr1[20];
			int n3 = (int) (countTotalIterations + countRun * numRuns);
			double n4 = getMSE(trainDataset, classes);
			std::sprintf(ferr1, "%d,%f", n3, n4);
			output(1, ferr1);

			/*
			 * Output 1.1.5 (ferr2): On the test dataset
			 * TODO
			 */

			/*
			 * Step 1.2: Perform one refinement operation
			 */

			refinement.free_refine(grid_->getStorage(), functor);

			alpha_->resizeZero(grid_->getSize());

			/*
			 * Output 1.2.1 (fgrid): Serialized grid structure
			 */

			std::string fgrid;
			grid_->serialize(fgrid);
			output(3, fgrid);
		}
	}

	/*
	 * Step 2: Perform CG
	 */

	sg::solver::ConjugateGradients *cg = new sg::solver::ConjugateGradients(
			CG_IMAX, CG_EPS);

	sg::base::OperationMatrix *C_ = sg::op_factory::createOperationIdentity(
			*this->grid_);
	sg::datadriven::DMSystemMatrix matrix(*grid_, trainDataset, *C_, lambda);

	sg::base::DataVector b(alpha_->getSize());
	matrix.generateb(classes, b);

	cg->solve(matrix, *alpha_, b, true, true);

	std::cout << getMSE(trainDataset, classes) << std::endl;
	isTrained_ = true;
}

void LearnerOnlineSGD::performSGDStep(sg::base::DataMatrix& trainDataset,
		sg::base::DataVector& classes, double lambda, double gamma) {
	using namespace sg::base;

	size_t numCoeff = grid_->getStorage()->size();
	size_t dim = trainDataset.getNcols();
	size_t numData = trainDataset.getNrows();

	// Get random x and y pair
	int k = getRandom((int) numData - 1);
	DataVector x(dim);
	trainDataset.getRow((size_t) k, x);
	double y = classes[k];

	pushSGDIndex(k);

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

double LearnerOnlineSGD::getMSE(sg::base::DataMatrix& trainDataset,
		sg::base::DataVector& classes) {
	using namespace sg::base;

	size_t numData = trainDataset.getNrows();

	// Error vector
	DataVector *error = new DataVector(numData);
	error->setAll(0.0);

	DataVector result(numData);
	sg::op_factory::createOperationMultipleEval(*grid_, &trainDataset)->mult(
			*alpha_, result);

	for (unsigned int i = 0; i < numData; i++)
		error->set(i, result[i] - classes[i]);

	// MSE
	double sum;
	for (unsigned int i = 0; i < numData; i++)
		sum += pow(fabs(error->get(i)), 2);

	return (sum / (double) numData);

}

int LearnerOnlineSGD::getRandom(int limit) {
	int divisor = RAND_MAX / (limit + 1);
	int r;

	do {
		r = rand() / divisor;
	} while (r > limit);

	return r;
}

void LearnerOnlineSGD::pushSGDIndex(size_t index) {
	batch[batchIndex] = index;
	batchIndex = (batchIndex + 1) % batchSize;
}

sg::base::DataVector* LearnerOnlineSGD::getAlpha() {
	if (alpha_ == NULL)
		throw base::application_exception("Not initialized");
	return alpha_;
}

double LearnerOnlineSGD::errorOnMinibatch(sg::base::DataMatrix& trainDataset,
		sg::base::DataVector& classes) {

	size_t dim = trainDataset.getNcols();

	sg::base::DataMatrix batchMatrix(batchSize, dim);
	sg::base::DataVector batchClasses(batchSize);

	for (size_t i = 0; i < batchSize; i++) {
		sg::base::DataVector row(dim);
		trainDataset.getRow(batch[i], row);
		batchMatrix.appendRow(row);
		batchClasses.append(classes[batch[i]]);
	}

	return getMSE(batchMatrix, batchClasses);
}

void LearnerOnlineSGD::output(int fd, std::string str) {
	if (isVerbose_) {
		std::cout << str << std::endl;
	}
	if (outputStreams != NULL) {
		*(outputStreams[fd]) << str << std::endl;
	}
}

LearnerOnlineSGD::~LearnerOnlineSGD() {
}
}
}
