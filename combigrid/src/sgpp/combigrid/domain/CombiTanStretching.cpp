// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <sgpp/combigrid/domain/CombiTanStretching.hpp>
#include <math.h>

#include <vector>

void combigrid::TanStretching::get1DStretching(int level, double min, double max,
                                               std::vector<double>* stretching,
                                               std::vector<double>* jacobian) const {
  int nrPoints = combigrid::powerOfTwo[level] + 1;
  stretching->resize(nrPoints);
  std::vector<double> tmpPoints(nrPoints);

  for (int ii = 0; ii < nrPoints; ii++)
    tmpPoints[ii] = static_cast<double>(2 * ii - combigrid::powerOfTwo[level]) /
                    static_cast<double>(combigrid::powerOfTwo[level]);

  for (int ii = 0; ii < nrPoints; ii++)
    // pi/2 = 1.57079632679490
    (*stretching)[ii] = tan((1.57079632679490 - intFact_) * tmpPoints[ii]);

  // do the scaling
  for (int ii = 0; ii < nrPoints; ii++)
    (*stretching)[ii] =
        min + (max - min) * 0.5 * (1 + ((*stretching)[ii] / (*stretching)[nrPoints - 1]));

  // PeTz
  // TODO(???): CALCULATE the tan stretching's JACOBIAN!!!
}
