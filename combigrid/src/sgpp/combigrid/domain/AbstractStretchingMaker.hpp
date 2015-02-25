// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#ifndef ABSTRACTSTRETCHINGMAKER_HPP_
#define ABSTRACTSTRETCHINGMAKER_HPP_

#include <sgpp/combigrid/utils/combigrid_ultils.hpp>

namespace combigrid {
  /** class to create stretching in 1D*/
  class AbstractStretchingMaker {
    public:
      /**
       * @param level [IN] level of the array
       * @param min [IN] minimum value of the domain
       * @param max [IN] maximum value of the domain
       * @param stretching [OUT] 2^level + 1 elements */
      virtual void get1DStretching(
        int level , double min, double max,
        std::vector<double>& stretching) const = 0;

  };
}

#endif /* ABSTRACTSTRETCHINGMAKER_HPP_ */