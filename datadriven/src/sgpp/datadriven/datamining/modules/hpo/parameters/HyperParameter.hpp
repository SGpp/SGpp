/*
 * HyperParameter.h
 *
 *  Created on: 24.01.2018
 *      Author: Eric Koepke
 */

#ifndef DATADRIVEN_SRC_SGPP_DATADRIVEN_DATAMINING_MODULES_HPO_HYPERPARAMETER_H_
#define DATADRIVEN_SRC_SGPP_DATADRIVEN_DATAMINING_MODULES_HPO_HYPERPARAMETER_H_

#include <sgpp/datadriven/datamining/modules/hpo/harmonica/ConfigurationBit.hpp>

#include <list>

namespace sgpp {
namespace datadriven {

/**
 * Class to represent a hyperparameter
 */
class HyperParameter {
public:
  /**
   * Default constructor for implicit use by data structures
   */
  HyperParameter() :bits(), nBits(0), name(){}
  /**
   * Normal constructor
   * @param nBits number of bits for representation in harmonica
   * @param name name of the hyperparameter
   */
	HyperParameter(int nBits, std::string& name)
          :bits(), nBits(nBits), name(name){}
  /**
   * Default Destructor
   */
	virtual ~HyperParameter() = default;

  /**
   * Connects parameter to bit representation in harmonica
   * @param configBits
   */
	void makeConfigBits(std::vector<ConfigurationBit*>& configBits);

  /**
   * sets value according to the associated bits
   */
	virtual void setHarmonica() = 0;

protected:
  /**
   * associated configuration bits for harmonica
   */
  std::vector<ConfigurationBit> bits;
  /**
   * number of bits for harmonica
   */
  int nBits;
  /**
   * name of the hyperparameter
   */
  std::string name;
};

} /* namespace datadriven */
} /* namespace sgpp */

#endif /* DATADRIVEN_SRC_SGPP_DATADRIVEN_DATAMINING_MODULES_HPO_HYPERPARAMETER_H_ */
