/*
 * CategoricalParameter.hpp
 *
 *  Created on: Jan 25, 2018
 *      Author: polarbart
 */

#ifndef CategoricalParameter_HPP_
#define CategoricalParameter_HPP_

#include "HyperParameter.hpp"
#include "ContinuousParameter.hpp"

namespace sgpp {
namespace datadriven {


class ExponentialParameter: public ContinuousParameter {
public:
  ExponentialParameter(std::string& name, double min, double max, int nBits) :
          ContinuousParameter(name, min, max, nBits) {}
  //~ExponentialParameter();

  double getValue() override;

};
} /* namespace datadriven */
} /* namespace sgpp */
#endif /* CategoricalParameter_HPP_ */
