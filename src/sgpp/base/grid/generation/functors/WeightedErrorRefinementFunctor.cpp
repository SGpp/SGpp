#include "base/grid/generation/functors/WeightedErrorRefinementFunctor.hpp"

namespace sg {
namespace base {

WeightedErrorRefinementFunctor::WeightedErrorRefinementFunctor(
		DataVector* alpha, Grid* grid, size_t refinements_num, double threshold) :
		alpha(alpha), refinements_num(refinements_num), threshold(threshold), grid(
				grid), trainDataset(
		NULL), classes(NULL) {

}

WeightedErrorRefinementFunctor::~WeightedErrorRefinementFunctor() {

}

/*
 * Given alpha_j, calculates
 *
 * \sum_{i=0}^{N} abs(phi_j(x_i) * alpha_j * r_i^2)
 *
 * r_i = y_i - \sum{j=0}^{M} phi_j(x_i) * alpha_j
 */
double WeightedErrorRefinementFunctor::operator()(GridStorage* storage,
		size_t seq) {

	if (trainDataset == NULL || classes == NULL) {
		throw base::application_exception(
				"Training dataset or classes not set");
	}

	double error = 0;
	const size_t MIN_SUPPORT = 10;

	size_t numData = trainDataset->getNrows();
	size_t dim = trainDataset->getNcols();

	size_t numSupport = 0;
	for (size_t i = 0; i < numData; i++) {
		DataVector row(dim);
		trainDataset->getRow(i, row);

		/* phi_j(x_i) * alpha_j */
		sg::base::DataVector singleAlpha(alpha->getSize());
		singleAlpha.setAll(0.0);
		singleAlpha.set(seq, alpha->get(seq));
		double val1 = sg::op_factory::createOperationEval(*grid)->eval(singleAlpha, row);

		/* r_i */
		double val2 = classes->get(seq) - sg::op_factory::createOperationEval(*grid)->eval(*alpha, row);

		/* abs(phi_j(x_i) * alpha_j * r_i^2) */
		error += fabs(val1 * val2 * val2);

		if (val1 != 0) {
			numSupport++;
		}
	}

	if (numSupport < MIN_SUPPORT) {
		return -1; // under the threshold of 0.0
	}

    //std::cout << "Functor value (of " << seq << "): " << error << std::endl;
	return error;
}

double WeightedErrorRefinementFunctor::start() {
	return 0.0;
}

size_t WeightedErrorRefinementFunctor::getRefinementsNum() {
	return this->refinements_num;
}

double WeightedErrorRefinementFunctor::getRefinementThreshold() {
	return this->threshold;
}

void WeightedErrorRefinementFunctor::setTrainDataset(
    DataMatrix* trainDataset_) {
	trainDataset = trainDataset_;
}

void WeightedErrorRefinementFunctor::setClasses( DataVector* classes_) {
	classes = classes_;
}

}
}
