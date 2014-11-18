#include <limits>

#include "base/grid/generation/functors/PersistentErrorRefinementFunctor.hpp"

namespace sg {
namespace base {

PersistentErrorRefinementFunctor::PersistentErrorRefinementFunctor(
		DataVector* alpha, Grid* grid, size_t refinements_num, double threshold) :
		alpha(alpha), refinements_num(refinements_num), threshold(threshold), grid(
				grid), trainDataset(
		NULL), classes(NULL), errors(NULL), accum(NULL) {

}
PersistentErrorRefinementFunctor::~PersistentErrorRefinementFunctor() {
}

double PersistentErrorRefinementFunctor::operator()(GridStorage* storage,
		size_t seq) {

	if (trainDataset == NULL || classes == NULL) {
		throw base::application_exception(
				"Training dataset or classes not set");
	}

	const double BETA = 0.1;
	const size_t MIN_SUPPORT = 5;

	size_t numData = trainDataset->getNrows();
	//size_t dim = trainDataset->getNcols();

	DataVector phi_x(numData);

	sg::base::DataVector singleAlpha(alpha->getSize());
	singleAlpha.setAll(0.0);
	singleAlpha.set(seq, 1.0);
	sg::op_factory::createOperationMultipleEval(*grid, trainDataset)->mult(
			singleAlpha, phi_x);

	if (phi_x.getNumberNonZero() < MIN_SUPPORT) {
		return start(); // threshold is 0.0
	}

	size_t numCoeff = alpha->getSize();

	if (accum == NULL) {
		accum = new sg::base::DataVector(numCoeff);
		accum->setAll(start());
	} else if (accum->getSize() != numCoeff) {
		accum->resizeZero(numCoeff);
	}

	double current_j = 0;
	double tmp = 0;
	for (size_t i = 0; i < numData; i++) {
		tmp = (errors->get(i)) * phi_x[i];
		current_j += tmp*tmp;

	}

	// Accumulation
	accum->set(seq, accum->get(seq) * (1-BETA) + BETA*current_j*fabs(alpha->get(seq)));

	double func_val = accum->get(seq);

	//std::cout << "Functor value (of " << seq << "): " << func_val << std::endl;
	return func_val;
}

double PersistentErrorRefinementFunctor::start() {
	return this->threshold;
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
void PersistentErrorRefinementFunctor::setErrors(DataVector* errors_) {
	errors = errors_;
}

}
}
