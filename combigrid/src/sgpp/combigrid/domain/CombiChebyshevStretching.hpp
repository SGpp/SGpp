/*
 * CombiChebyshevStretching.hpp
 *
 *  Created on: 28 Jul 2014
 *      Author: kenny
 *
 */

#ifndef COMBICHEBYSHEVSTRETCHING_HPP_
#define COMBICHEBYSHEVSTRETCHING_HPP_

/**
 *   Implements a simple coordinate transform :
 *   T : [-1;1] -> [a;b]
 *
 *   y(x)  = 0.5*(b+a) - 0.5*(b-a)cos(Pi*x) , or in the discrete form:
 *   y_j = 0.5*(b+a) - 0.5*(b-a)*cos(Pi*j/(N-1)); , where N is the nr of points
 *on the discrete grid...
 *  The class is called ChebyshevStretching, well because it results in a
 *non-equidistant or "stretched" grid, and Chebyshev, because
 *  the points:
 *  cos(Pi*j/(N)) are the extrema of Chebyshev polynomial T_N(x), defined
 *beautifully as  T_N(x) = cos(N*arcos(x))
 **/

#include <sgpp/combigrid/domain/AbstractStretchingMaker.hpp>

namespace combigrid {

class CombiChebyshevStretching : public AbstractStretchingMaker {
 public:
  CombiChebyshevStretching() : AbstractStretchingMaker() { ; }

  virtual ~CombiChebyshevStretching() { ; }
  /**
   * @param level - integer specifying the current grid level . the
   *corresponding nr of points is 2^level + 1
   * @param min - the left boundary of the interval
   * @param max - the right boundary of the interval
   * @param stretching - the output vector of pre-computed grid points...
   * @param jacobian - the evaluated jacobian at all points of the stretching ,
   *taking into consideration
   * size of the interval and underlying tranformations.
   *
   */
  void get1DStretching(int level, double min, double max,
                       std::vector<double>& stretching,
                       std::vector<double>& jacobian) const;

  Stretching getStretchingType() const { return CHEBYSHEV; }
};
}

#endif /* COMBICHEBYSHEVSTRETCHING_HPP_ */
