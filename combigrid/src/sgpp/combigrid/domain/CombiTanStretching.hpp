/* ****************************************************************************
* Copyright (C) 2011 Technische Universitaet Muenchen                         *
* This file is part of the SG++ project. For conditions of distribution and   *
* use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
**************************************************************************** */
// @author Janos Benk (benk@in.tum.de)

#ifndef COMBITANSTRETCHING_HPP_
#define COMBITANSTRETCHING_HPP_

#include "combigrid/domain/AbstractStretchingMaker.hpp"

namespace combigrid {

  /** The stretching function of this class is: <br>
   *octave:28> intFact = 1/4; <br>
    octave:29> x=(-pi/2+intFact):2^-L:(pi/2-intFact); <br>
    octave:30> y = tan(x); <br>
   * */
  class TanStretching : public AbstractStretchingMaker {
    public:
      /** Ctor
       * @param intFact must be smaller than one*/
      TanStretching(double intFact = 1.0 / 7.0 ): AbstractStretchingMaker() , intFact_(intFact) {
        if (intFact_ > 1.5) intFact_ = 1.0 / 1.5;

        if (intFact_ < 0.01) intFact_ = 1.0 / 10.0;
      }

      virtual ~TanStretching() {
        ;
      }

      void get1DStretching(
        int level , double min, double max,
        std::vector<double>& stretching) const;

    private:

      /** internal factor for the formula */
      double intFact_;

  };

}

#endif /* COMBITANSTRETCHING_HPP_ */
