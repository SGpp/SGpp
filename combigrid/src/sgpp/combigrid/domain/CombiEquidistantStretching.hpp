/*
 * CombiEquidistantStretching.hpp
 *
 *  Created on: 15 Sep 2014
 *      Author: kenny
 */

#ifndef COMBIEQUIDISTANTSTRETCHING_HPP_
#define COMBIEQUIDISTANTSTRETCHING_HPP_
/**
 *
 * Create an equidistant  grid..
 *
 *
 */

#include <sgpp/combigrid/domain/AbstractStretchingMaker.hpp>

namespace combigrid {

class CombiEquidistantStretching : public AbstractStretchingMaker {
 public:
  //  CombiEquidistantStretching(){}

  ~CombiEquidistantStretching() {}

  /**
   *
   * @param level integer specifying the current grid level, the
   * corresponding nr of points is 2**level + 1
   * @param min the left boundary of the interval
   * @param max the right boundary of the interval
   * @param stretching the output vector of pre-computed grid points
   * @param jacobian the evaluated jacobian at all points of the stretching,
   * taking into consideration size of the interval and underlying tranformations
   */
  void get1DStretching(int level, double min, double max,
                       std::vector<double>* stretching,
                       std::vector<double>* jacobian) const;

  Stretching getStretchingType() const { return EQUIDISTANT; }
};
}

#endif /* COMBIEQUIDISTANTSTRETCHING_HPP_ */
