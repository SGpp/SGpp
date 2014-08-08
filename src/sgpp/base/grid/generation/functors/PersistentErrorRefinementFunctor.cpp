#include "base/grid/generation/functors/PersistentErrorRefinementFunctor.hpp"

namespace sg {
namespace base {

PersistentErrorRefinementFunctor::PersistentErrorRefinementFunctor(
		DataVector* alpha, Grid* grid, size_t refinements_num, double threshold) :
		alpha(alpha), refinements_num(refinements_num), threshold(threshold), grid(
				grid), trainDataset(
		NULL), classes(NULL), error(NULL) {

}
PersistentErrorRefinementFunctor::~PersistentErrorRefinementFunctor() {
}

double PersistentErrorRefinementFunctor::operator()(GridStorage* storage,
		size_t seq) {

	const double BETA=0.9;

	if (trainDataset == NULL || classes == NULL) {
		throw base::application_exception(
				"Training dataset or classes not set");
	}

	// Make sure that the error vector is as large as
	// the coefficient vector
	size_t numCoeff = alpha->getSize();

	if (error == NULL) {
		error = new sg::base::DataVector(numCoeff);
		error->setAll(0.0);
	} else if (error->getSize() != numCoeff) {
		error->resizeZero(numCoeff);
	}

	// Calculate the weighted error for each basis function
	sg::base::DataVector* current = new sg::base::DataVector(numCoeff);
	current->setAll(0.0);

	for (size_t i = 0; i < numCoeff; i++) {
		current->set(i, calcWeightedError(i));
	}

	// Combine the current error vector with the existing error
	// vector
	for (size_t i = 0; i < numCoeff; i++ ) {
		error->set(i, error->get(i) * BETA + current->get(i));
	}

	return error->get(seq);
}

/*
 * See WeightedErrorRefinementFunctor for a description
 */
double PersistentErrorRefinementFunctor::calcWeightedError(size_t seq) {

	if (trainDataset == NULL || classes == NULL) {
		throw base::application_exception(
				"Training dataset or classes not set");
	}

	double error = 0;

	size_t numData = trainDataset->getNrows();
	size_t dim = trainDataset->getNcols();

	for (size_t i = 0; i < numData; i++) {
		DataVector row(dim);
		trainDataset->getRow(i, row);

		/* Hack to calculate phi_j(x_i) * alpha_j */
		sg::base::DataVector singleAlpha(alpha->getSize());
		singleAlpha.setAll(0.0);
		singleAlpha.set(seq, alpha->get(seq));

		double val = sg::op_factory::createOperationEval(*grid)->eval(
				singleAlpha, row);
		double err = classes->get(seq) - val;

		error += val * err;
	}

	return error;
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
