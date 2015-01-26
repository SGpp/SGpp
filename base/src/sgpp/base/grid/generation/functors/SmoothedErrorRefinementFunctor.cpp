/******************************************************************************
 * Copyright (C) 2010 Technische Universitaet Muenchen                         *
 * This file is part of the SG++ project. For conditions of distribution and   *
 * use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
 ******************************************************************************/
// @author Dirk Pflueger (pflueged@in.tum.de), Jörg Blank (blankj@in.tum.de), Alexander Heinecke (Alexander.Heinecke@mytum.de)
#include <sgpp/base/grid/generation/functors/SmoothedErrorRefinementFunctor.hpp>
#include <sgpp/base/operation/BaseOpFactory.hpp>

#include <sgpp/globaldef.hpp>


namespace SGPP {
namespace base {

SmoothedErrorRefinementFunctor::SmoothedErrorRefinementFunctor(
		DataVector* alpha, Grid* grid, size_t refinements_num, double threshold) :
		alpha(alpha), refinements_num(refinements_num), threshold(threshold), grid(
				grid), trainDataset(
		NULL), classes(NULL) {

}

SmoothedErrorRefinementFunctor::~SmoothedErrorRefinementFunctor() {

}

double SmoothedErrorRefinementFunctor::operator()(GridStorage* storage,
		size_t seq) {

	if (trainDataset == NULL || classes == NULL) {
		throw base::application_exception(
				"Training dataset or classes not set");
	}

	std::vector<double> errors;

	size_t numData = trainDataset->getNrows();
	size_t dim = trainDataset->getNcols();

	DataVector unit(alpha->getSize());
	unit.setAll(0.0);
	unit[seq] = 1.0;

	for (size_t i = 0; i < numData; i++) {
		DataVector row(dim);
		trainDataset->getRow(i, row);

		// check if this particular training data is in the support
		double tmp = SGPP::op_factory::createOperationEval(*grid)->eval(unit,
				row);
		if (tmp == 0) {
			continue;
		}

		// if it is inside the support,
		// calculate error and add it to
		// the vector
		double err = classes->get(i)
				- SGPP::op_factory::createOperationEval(*grid)->eval(*alpha, row);
		errors.push_back(err);
	}

	// Uniform weight
	double SmoothedErrors = 0;
	size_t numErrors = errors.size();

	if (numErrors == 0) {
		return 0.0;
	}

	for (int i = 0; i < numErrors; i++) {
		SmoothedErrors += errors[i] * (1.0 / numErrors);
	}

	return SmoothedErrors;
}

double SmoothedErrorRefinementFunctor::start() {
	return 0.0;
}

size_t SmoothedErrorRefinementFunctor::getRefinementsNum() {
	return this->refinements_num;
}

double SmoothedErrorRefinementFunctor::getRefinementThreshold() {
	return this->threshold;
}

void SmoothedErrorRefinementFunctor::setTrainDataset(
		DataMatrix* trainDataset_) {
	trainDataset = trainDataset_;
}

void SmoothedErrorRefinementFunctor::setClasses(DataVector* classes_) {
	classes = classes_;
}

}
}
