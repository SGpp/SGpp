// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at 
// sgpp.sparsegrids.org

#include <sgpp/quadrature/sample/SampleGenerator.hpp>

using namespace SGPP::base;

#include <sgpp/globaldef.hpp>


namespace SGPP {
namespace quadrature {

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