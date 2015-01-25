/* ****************************************************************************
 * Copyright (C) 2014 Universitaet Stuttgart                                   *
 * This file is part of the SG++ project. For conditions of distribution and   *
 * use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
 **************************************************************************** */
// @author Andreas Doerr, Marcel Schneider, Matthias Moegerle
#include "SampleGenerator.hpp"

using namespace sg::base;

namespace sg {
namespace mcm {

void SampleGenerator::getSamples(DataMatrix& samples) {

	// Number of columns has to correspond to the number of dimensions
	if (samples.getNcols() != dimensions)
		return;

	// generate one sample for every row of the given DataMatrix
	for (size_t i = 0; i < samples.getNrows(); i++) {
		DataVector dv(dimensions);
		getSample(dv);
		samples.setRow(i, dv);
	}
}

size_t SampleGenerator::getDimensions() {
	return dimensions;
}

void SampleGenerator::setDimensions(size_t dimensions) {
	this->dimensions = dimensions;
}

}
}
