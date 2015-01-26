/* ****************************************************************************
* Copyright (C) 2011 Technische Universitaet Muenchen                         *
* This file is part of the SG++ project. For conditions of distribution and   *
* use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
**************************************************************************** */
// @author Janos Benk (benk@in.tum.de)

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
