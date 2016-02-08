/*
 * LinearLoadModel.hpp
 *
 *  Created on: Oct 9, 2013
 *      Author: heenemo
 */

#ifndef LINEARLOADMODEL_HPP_
#define LINEARLOADMODEL_HPP_

#include <cmath>
#include <iostream>
#include "sgpp/distributedcombigrid/utils/LevelVector.hpp"
#include "sgpp/distributedcombigrid/utils/Types.hpp"
#include "sgpp/distributedcombigrid/loadmodel/LoadModel.hpp"

namespace combigrid {

class LinearLoadModel: public LoadModel {
 public:
  LinearLoadModel() {
  }

  inline real eval(const LevelVector& l) const;

  virtual ~LinearLoadModel() {
  }
};

inline real LinearLoadModel::eval(const LevelVector& l) const {
  // todo does not hold in general case
  real ret(1.0);

  for (size_t i = 0; i < l.size(); ++i) {
    ret *= std::pow(real(2.0), static_cast<real>(l[i]));
  }

  return ret;
}

} /* namespace combigrid */
#endif /* LINEARLOADMODEL_HPP_ */
