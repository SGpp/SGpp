#include "base/grid/generation/functors/PersistentErrorRefinementFunctor.hpp"

namespace sg {
namespace base {

PersistentErrorRefinementFunctor::PersistentErrorRefinementFunctor(
		DataVector* alpha, Grid* grid, size_t refinements_num, double threshold) :
		alpha(alpha), refinements_num(refinements_num), threshold(threshold), grid(
				grid), trainDataset(
		NULL), classes(NULL), accum(NULL) {

}
PersistentErrorRefinementFunctor::~PersistentErrorRefinementFunctor() {
}

/*
 * current[j] = \sum_{i=0}^{N} (r_i + y_i) * \phi_j(x_i)
 * accum[j] = BETA * accum[j] + current[j]
 * functor value = -alpha_j * accum[j]
 */
double PersistentErrorRefinementFunctor::operator()(GridStorage* storage,
		size_t seq) {

	if (trainDataset == NULL || classes == NULL) {
		throw base::application_exception(
				"Training dataset or classes not set");
	}

	const double BETA = 0.9;
	const size_t MIN_SUPPORT = 10;

	size_t numData = trainDataset->getNrows();
	size_t dim = trainDataset->getNcols();

	// Check support
	size_t numSupport = 0;

	DataVector row(dim);

	for (size_t i = 0; i < numData; i++) {
		/* store x_i */
		trainDataset->getRow(i, row);

		/* phi_{seq}(x_i) */
		sg::base::DataVector singleAlpha(alpha->getSize());
		singleAlpha.setAll(0.0);
		singleAlpha.set(seq, 1);
		double phi = sg::op_factory::createOperationEval(*grid)->eval(
				singleAlpha, row);

		if (phi > 0) {
			numSupport++;
		}
		if (numSupport == MIN_SUPPORT) {
			break;
		}
	}

	if (numSupport != MIN_SUPPORT) {
		return -1; // threshold is 0.0
	}

	// Make sure that the error vector is as large as
	// the coefficient vector
	size_t numCoeff = alpha->getSize();

	if (accum == NULL) {
		accum = new sg::base::DataVector(numCoeff);
		accum->setAll(0.0);
	} else if (accum->getSize() != numCoeff) {
		accum->resizeZero(numCoeff);
	}

	// Calculate current
	sg::base::DataVector* current = new sg::base::DataVector(numCoeff);
	current->setAll(0.0);

	for (size_t j = 0; j < numCoeff; j++) {

		// Calculate
		// current[j] = \sum_{i=0}^{N} (r_i + y_i) * \phi_j(x_i)

		double current_j = 0;

		for (size_t i = 0; i < numData; i++) {
			/* store x_i */
			trainDataset->getRow(i, row);

			/* phi_j(x_i) */
			sg::base::DataVector singleAlpha(alpha->getSize());
			singleAlpha.setAll(0.0);
			singleAlpha.set(j, 1);
			double val1 = sg::op_factory::createOperationEval(*grid)->eval(
					singleAlpha, row);

			/* r_i */
			double val2 = classes->get(i)
					- sg::op_factory::createOperationEval(*grid)->eval(*alpha,
							row);

			/* (r_i + y_i) * \phi_j(x_i) */
			current_j += (val2 + classes->get(i)) * val1;
		}

		current->set(j, current_j);
	}

	// Accumulation
	for (size_t i = 0; i < numCoeff; i++) {
		accum->set(i, accum->get(i) * BETA + current->get(i));
	}

	//delete current;

	double func_val = -alpha->get(seq) * accum->get(seq);

	// std::cout << "Functor value (of " << seq << "): " << func_val << std::endl;
	return func_val;
}

double PersistentErrorRefinementFunctor::start() {
	return 0.0;
}

size_t PersistentErrorRefinementFunctor::getRefinementsNum() {
	return this->refinements_num;
}

double PersistentErrorRefinementFunctor::getRefinementThreshold() {
	return this->threshold;
}

void PersistentErrorRefinementFunctor::setTrainDataset(
		DataMatrix* trainDataset_) {
	trainDataset = trainDataset_;
}

void PersistentErrorRefinementFunctor::setClasses(DataVector* classes_) {
	classes = classes_;
}

}
}
