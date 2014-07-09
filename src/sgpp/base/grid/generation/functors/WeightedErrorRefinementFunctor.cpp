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

double WeightedErrorRefinementFunctor::operator()(GridStorage* storage,
		size_t seq) {

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

		double val = sg::op_factory::createOperationEval(*grid)->eval(singleAlpha, row);
		double err = classes->get(seq) - val;

		error += val * err * err;
	}

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

void WeightedErrorRefinementFunctor::setClasses(DataVector* classes_) {
	classes = classes_;
}

}
}
