// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at 
// sgpp.sparsegrids.org

#include "CombiUniformStretching.hpp"


void combigrid::UniformStretching::get1DStretching(
  int level , double min, double max,
  std::vector<double>& stretching) const {

  int nrPoints = combigrid::powerOfTwo[level] + 1;
  stretching.resize( nrPoints , 0.0 );

  // set each point just uniformly
  for (int i = 0 ; i < nrPoints ; i++) {
    stretching[i] = min + ((double)i) * (max - min) / (nrPoints - 1);
  }
}