/*
 * CombiChebyshevStretching.cpp
 *
 *  Created on: 28 Jul 2014
 *      Author: kenny
 */

#include <sgpp/combigrid/domain/CombiChebyshevStretching.hpp>
#include <sgpp/globaldef.hpp>

/**
 *   Implements a simple coordinate transform :
 *   T : [-1;1] -> [a;b]
 *   y(x)  = 0.5*(b+a) - 0.5*(b-a)cos(Pi*x) , or in the discrete form:
 *   y_j = 0.5*(b+a) - 0.5*(b-a)*cos(Pi*j/(N-1)); , where N is the nr of points
 *on the discrete grid...
 *  The class is called ChebyshevStretching, well because it results in a
 *non-equidistant or "stretched" grid, and Chebyshev, because
 *  the points:
 *  cos(Pi*j/(N)) are the extrema of Chebyshev polynomial T_N(x), defined
 *beautifully as  T_N(x) = cos(N*arcos(x))
 **/

void combigrid::CombiChebyshevStretching::get1DStretching(
    int level, double min, double max, std::vector<double>* stretching,
    std::vector<double>* jacobian) const {
  int nr_pts = powerOfTwo[level] + 1;
  // resize the stretching vector...
  stretching->resize(nr_pts);
  jacobian->resize(nr_pts);

  double factor = M_PI / (nr_pts - 1.0);

  TRANSFORMATION_TYPE type = FINITE;

  if (min == n_INF && max < p_INF) type = SEMI_INFINITE_NINF;

  if (min > n_INF && max == p_INF) type = SEMI_INFINITE_PINF;

  if (min == n_INF && max == p_INF) type = INFINITE;

  (*stretching)[0] = min;
  (*jacobian)[0] = transformationJacobian(min, max, -1.0, type);
  int loop_len = nr_pts - 1;

  for (int i = 1; i < loop_len; i++) {
    // transform the stretching with respect to the integration interval...n
    // save the transformation jacobian...
    (*stretching)[i] = transforminterval(min, max, -cos(i * factor), type);
    (*jacobian)[i] = transformationJacobian(min, max, -cos(i * factor), type);
  }

  (*stretching)[loop_len] = max;  // seems to be correct ...
  (*jacobian)[loop_len] = transformationJacobian(min, max, +1.0, type);
}
