/*
 * SampleProvider.hpp
 *
 *  Created on: Feb 8, 2016
 *      Author: perun
 */

#ifndef SAMPLEPROVIDER_HPP_
#define SAMPLEPROVIDER_HPP_

#include <memory>

#include <sgpp/datadriven/tools/Dataset.hpp>

#include <sgpp/globaldef.hpp>


namespace SGPP {
namespace datadriven {

class SampleProvider {
public:
	SampleProvider(){};
	virtual ~SampleProvider(){};

	/**
	 * Selects a certain number of samples
	 * @param how_many number of samples to return
	 * @param dataset pointer to the returned dataset
	 */
	virtual Dataset nextSamples(int how_many) = 0;

	/**
	 * Returns all samples
	 * @param dataset pointer to the returned dataset
	 */
	virtual Dataset allSamples() = 0;

};
}
}

#endif /* SAMPLEPROVIDER_HPP_ */
