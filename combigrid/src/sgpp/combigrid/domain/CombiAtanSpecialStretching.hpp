/* ****************************************************************************
* Copyright (C) 2011 Technische Universitaet Muenchen                         *
* This file is part of the SG++ project. For conditions of distribution and   *
* use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
**************************************************************************** */
// @author Janos Benk (benk@in.tum.de)

#ifndef COMBIATANSPECIALSTRETCHING_HPP_
#define COMBIATANSPECIALSTRETCHING_HPP_

#include <sgpp/combigrid/domain/AbstractStretchingMaker.hpp>

namespace combigrid {

  /** Stretching formula with the following matlab formula: <br>
   *  L = 6; <br>
   *  N = -1:(1/2^L):1; <br>
   *  Fakt = 5.3; <br>
   *  V = tan(((pi/2) - 1/Fakt)*N); <br>
   *  V1 = 3*atan( (N).^5+0.1*(N) ); <br>
   *  V = V ./ max(V); <br>
   *  V1 = V1 ./ max(V1); <br>
  */
  class AtanSpecialStretching : public AbstractStretchingMaker {
    public:

      AtanSpecialStretching(): AbstractStretchingMaker() {
        ;
      }

      virtual ~AtanSpecialStretching() {
        ;
      }

      void get1DStretching(
        int level , double min, double max,
        std::vector<double>& stretching) const;
  };

}

#endif /* COMBIATANSPECIALSTRETCHING_HPP_ */
