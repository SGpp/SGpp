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
#include <sgpp/datadriven/datamining/DataMiningConfiguration.hpp>

#include <sgpp/globaldef.hpp>

namespace sgpp {
namespace datadriven {

class SampleProvider {
 public:
  SampleProvider(datadriven::DataMiningConfiguration& config){};
  virtual ~SampleProvider(){};

  /**
   * Selects a certain number of samples
   * @param how_many number of samples to return
   * @return dataset
   */
  virtual Dataset nextSamples(int how_many) = 0;

  /**
   * Returns all samples
   * @return dataset
   */
  virtual Dataset allSamples() = 0;

  /**
   * Returns the dimensionality of the data source
   * @return dimensionality
   */
  size_t getDim() { return dim; }

 protected:
  size_t dim;
};
}
}

#endif /* SAMPLEPROVIDER_HPP_ */
