/* ****************************************************************************
* Copyright (C) 2009 Technische Universitaet Muenchen                         *
* This file is part of the SG++ project. For conditions of distribution and   *
* use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
**************************************************************************** */
// @author Jörg Blank (blankj@in.tum.de), Alexander Heinecke (Alexander.Heinecke@mytum.de), Dirk Pflueger (pflueged@in.tum.de)

#ifndef MODIFIED_WAVELET_BASE_HPP
#define MODIFIED_WAVELET_BASE_HPP

#include <cmath>
#include <sgpp/base/basis/Basis.hpp>


namespace sg {
  namespace base {

    /**
     * Modified Wavelet basis functions.
     * @todo (pflueged) Concernc also most other basis functions: Why static cast in int and not in unsigned int?? If index is unsigned int (which is the case) then for level 32, indices get lost/corrupted.
     */
    template<class LT, class IT>
    class ModifiedWaveletBasis: public Basis<LT, IT> {
      public:
        double eval(LT level, IT index, double p) {
          //std::cout << "Level " << level <<" Index "<<index<<" Point "<<p<<" BasisValue ";
          if (level == 1) {
            //std::cout<<1<<std::endl;
            return 1;
          }

          double sup = 1 << level;
          sup = 1 / sup;

          if (static_cast<int>(index) == static_cast<int>((1 << level) - 1) && p > (1 - 1.5602 * sup) ) {
            //std::cout<<0.5014+1.3803*(0.5602+p/sup-index)<<std::endl;
            return 0.5014 + 1.3803 * (0.5602 + p / sup - index);
          }

          if ( index == 1 && p < 1.5602 * sup ) {
            //std::cout<<0.5014+1.3803*(0.5602-p/sup+1)<<std::endl;
            return 0.5014 + 1.3803 * (0.5602 - p / sup + 1);
          }

          double t = pow(p / sup - index, 2);
          //std::cout<<(1-t)*exp(-0.9986*t)<<std::endl;
          return (1 - t) * exp(-0.9986 * t);

        }

    };
    // default type-def (unsigned int for level and index)
    typedef ModifiedWaveletBasis<unsigned int, unsigned int> SModWaveletBase;
  }
}

#endif /* MODIFIED_WAVELET_BASE_HPP */
