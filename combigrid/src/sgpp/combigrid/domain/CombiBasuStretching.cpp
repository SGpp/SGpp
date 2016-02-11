// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <sgpp/combigrid/domain/CombiBasuStretching.hpp>

#include <vector>

void combigrid::CombiBasuStretching::get1DStretching(int level, double min, double max,
                                                     std::vector<double>* stretching,
                                                     std::vector<double>* jacobian) const {
  int N = powerOfTwo[level];  // N is the number of subintervals... N+1 the
                              // number of points!
  jacobian->clear();
  jacobian->resize(N + 1, 1.0);
  stretching->clear();

  TRANSFORMATION_TYPE type = FINITE;

  if (min == n_INF && max < p_INF) type = SEMI_INFINITE_NINF;

  if (min > n_INF && max == p_INF) type = SEMI_INFINITE_PINF;

  if (min == n_INF && max == p_INF) type = INFINITE;

  // multiplication factor for the abscissas which takes values 1/-1 depending
  // on the type of semi-infinite intervals we are looking at
  double factor = 1.0;
  // term translating the absissas from 0;inf to arbitrary integration interval
  double add = min;

  switch (type) {
    case FINITE:
      COMBIGRID_OUT_WRN(
          "Basu stretching over finite domains not supported! Returning zero "
          "array!!",
          __FILE__, __LINE__)
      factor = 0.0;
      return;

    case SEMI_INFINITE_NINF:
      factor = -1.0;
      add = max;
      break;

    case SEMI_INFINITE_PINF:
      factor = +1.0;
      add = min;
      break;

    case INFINITE:
      factor = -1.0;

      // take negative abscissas
      for (int s = N / 2; s >= 0; s--) {
        double sec_ = 1 / cos(M_PI * s / (N));
        double t_ = log(sec_ * sec_);
        stretching->push_back(factor * t_);
      }

      // take positive abscissas
      factor = 1.0;

      for (int s = 1; s <= N / 2; s++) {
        double sec_ = 1 / cos(M_PI * s / (N));
        double t_ = log(sec_ * sec_);
        stretching->push_back(factor * t_);
      }

      return;
      break;

    default:
      COMBIGRID_OUT_WRN(
          "Basu stretching over finite domains not supported! Returning zero "
          "array!!",
          __FILE__, __LINE__)
      factor = 0.0;
      return;
  }

  for (int s = 0; s <= N; s++) {
    double sec_ = 1 / cos(M_PI * s / (2 * N));
    double t_ = log(sec_ * sec_);
    stretching->push_back(add + factor * t_);
  }
}
