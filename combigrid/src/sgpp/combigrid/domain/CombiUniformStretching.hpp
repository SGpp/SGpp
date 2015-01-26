// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at 
// sgpp.sparsegrids.org

#ifndef COMBIUNIFORMSTRETCHING_HPP_
#define COMBIUNIFORMSTRETCHING_HPP_

#include <sgpp/combigrid/domain/AbstractStretchingMaker.hpp>

namespace combigrid {

  /** uniform stretching used only for testing purposes */
  class UniformStretching : public AbstractStretchingMaker {
    public:

      UniformStretching(): AbstractStretchingMaker() {
        ;
      }

      virtual ~UniformStretching() {
        ;
      }

      void get1DStretching(
        int level , double min, double max,
        std::vector<double>& stretching) const;

  };

}

#endif /* COMBIUNIFORMSTRETCHING_HPP_ */