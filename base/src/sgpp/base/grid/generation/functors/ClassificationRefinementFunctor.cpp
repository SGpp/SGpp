#include <sgpp/base/grid/generation/functors/ClassificationRefinementFunctor.hpp>

#include <sgpp/globaldef.hpp>


namespace SGPP {
namespace base {

ClassificationRefinementFunctor::ClassificationRefinementFunctor(
		DataVector* alpha, Grid* grid, size_t refinements_num, double threshold) :
		alpha(alpha), refinements_num(refinements_num), threshold(threshold),  grid(grid), trainDataset(NULL), classes(NULL) {
}

ClassificationRefinementFunctor::~ClassificationRefinementFunctor() {
}

double ClassificationRefinementFunctor::operator()(GridStorage* storage,
		size_t seq) {

	if (trainDataset == NULL || classes == NULL) {
		throw base::application_exception(
				"Training dataset or classes not set");
	}

	const size_t MIN_SUPPORT = 10;

	size_t numData = trainDataset->getNrows();
	size_t dim = trainDataset->getNcols();

	size_t numSupport = 0;
	size_t numMisclassified = 0;

	for (size_t i = 0; i < numData; i++) {
		/* store x_i */
		DataVector row(dim);
		trainDataset->getRow(i, row);

		/* phi_{seq}(x_i) */
		SGPP::base::DataVector singleAlpha(alpha->getSize());
		singleAlpha.setAll(0.0);
		singleAlpha.set(seq, 1);
		double phi = SGPP::op_factory::createOperationEval(*grid)->eval(
				singleAlpha, row);

		if (phi > 0) {
			numSupport++;

			/* evaluate grid at point x_i */
			double val = SGPP::op_factory::createOperationEval(*grid)->eval(*alpha, row);
			double y = classes->get(i);

			if ( ! ((y > 0 && val > 0) || (y < 0 && val < 0)) ) {
				numMisclassified++;
			}
		}
	}

	if (numSupport < MIN_SUPPORT) {
		return -1; // threshold is 0.0
	}

	return (double) numMisclassified;
}

double ClassificationRefinementFunctor::start() {
	return 0.0;
}

size_t ClassificationRefinementFunctor::getRefinementsNum() {
	return this->refinements_num;
}

double ClassificationRefinementFunctor::getRefinementThreshold() {
	return this->threshold;
}

void ClassificationRefinementFunctor::setTrainDataset(
		DataMatrix* trainDataset_) {
	trainDataset = trainDataset_;
}

void ClassificationRefinementFunctor::setClasses(DataVector* classes_) {
	classes = classes_;
}


}
}
