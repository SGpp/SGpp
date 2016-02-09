/*
 * CombiEquidistantStretching.cpp
 *
 *  Created on: 15 Sep 2014
 *      Author: kenny
 */
#include <sgpp/combigrid/domain/CombiEquidistantStretching.hpp>

/**
 *
 *
 *  simple 1D equidistant grid on the interval (min;max)
 **/

void combigrid::CombiEquidistantStretching::get1DStretching(
    int level, double min, double max, std::vector<double>& stretching,
    std::vector<double>& jacobian) const {
  int nr_pts = powerOfTwo[level] + 1;
  // resize the stretching vector...
  stretching.resize(nr_pts);
  jacobian.resize(nr_pts);

  // initialize the uniform grid on the interval to -1;1 and depending on
  // the type of the interval (-inf,a] [a,inf) (-inf; inf) or [a;b]
  // transform the points back to [a;b]
  // at the same time store this transformations's jacobian
  // i.e. dx(t)/dt where x(t) : [-1;1] -> (a,b)

  TRANSFORMATION_TYPE type = FINITE;

  if (min == n_INF && max < p_INF) type = SEMI_INFINITE_NINF;

  if (min > n_INF && max == p_INF) type = SEMI_INFINITE_PINF;

  if (min == n_INF && max == p_INF) type = INFINITE;

  double dt = (2.0) / double(nr_pts - 1);
  // notice that dt is the spacing on the [-1;1] interval...
  //
  stretching[0] = min;
  jacobian[0] = transformationJacobian(min, max, -1.0, type);
  int loop_len = nr_pts - 1;

  for (int i = 1; i < loop_len; i += 1) {
    stretching[i] = transforminterval(min, max, -1.0 + i * dt, type);
    jacobian[i] = transformationJacobian(min, max, -1.0 + i * dt, type);
  }

  stretching[loop_len] = max;  // seems to be correct ...
  jacobian[loop_len] = transformationJacobian(min, max, +1.0, type);
}
