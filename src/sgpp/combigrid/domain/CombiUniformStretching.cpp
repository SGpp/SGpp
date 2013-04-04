/* ****************************************************************************
* Copyright (C) 2011 Technische Universitaet Muenchen                         *
* This file is part of the SG++ project. For conditions of distribution and   *
* use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
**************************************************************************** */
// @author Janos Benk (benk@in.tum.de)

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
