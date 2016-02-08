/*
 * LoadModel.hpp
 *
 *  Created on: Oct 9, 2013
 *      Author: heenemo
 */

#ifndef LOADMODEL_HPP_
#define LOADMODEL_HPP_

#include "sgpp/distributedcombigrid/utils/Types.hpp"
#include "sgpp/distributedcombigrid/utils/LevelVector.hpp"

namespace combigrid {

class LoadModel {
 public:
  LoadModel() {
  }
  ;

  virtual real eval(const LevelVector& l) const = 0;

  virtual ~LoadModel() {
  }
  ;
};

} /* namespace combigrid */
#endif /* LOADMODEL_HPP_ */
