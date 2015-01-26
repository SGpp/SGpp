/* ****************************************************************************
* Copyright (C) 2011 Technische Universitaet Muenchen                         *
* This file is part of the SG++ project. For conditions of distribution and   *
* use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
**************************************************************************** */
// @author Janos Benk (benk@in.tum.de)

#ifndef PROLONGATIONRESTRICTION_HPP_
#define PROLONGATIONRESTRICTION_HPP_

#include "combigrid.hpp"
#include <sgpp/combigrid/combigrid/AbstractCombiGrid.hpp>

namespace combigrid {

  /** class which has two static methods for restrictions and prolongation of the multigrid method. <br>
   * */
  class ProlongationRestriction {
    public:

      /** empty Ctor*/
      ProlongationRestriction() {
        ;
      }

      /** empty Dtor*/
      virtual ~ProlongationRestriction() {
        ;
      }

      /** prolongation with the specified coefficients
         *
       * @param fgFine fine grid to be restricted
         * @param vectFine vector containing the elements of the fine full grid
         * @param coefFine coefficient with which the values from the fine grid should be weighted
         * @param fgCoarse resulting coarse grid
         * @param vectCoarse vector containing the elements of the coarse grid
         * @param coefCoarse coefficient with which the elements from the coarse grid shall be weigthed
         * @param nrSpace number of unknowns per node
         */
      static void prolongation(const FullGridD* fgFine ,
                               std::vector<double>& vectFine ,
                               double coefFine ,
                               const FullGridD* fgCoarse ,
                               const std::vector<double>& vectCoarse ,
                               double coefCoarse ,
                               int nrSpace) ;

      /** restriction with the specified coefficients
         *
       * @param fgFine fine grid which shall be restricted
       * @param vectFine vector containing the elements of the fine full grid
       * @param coefFine coefficient with which the values from the fine grid should be weighted
       * @param fgCoarse resulting coarse grid
       * @param vectCoarse vector containing the elements of the coarse grid
       * @param coefCoarse coefficient with which the elements from the coarse grid shall be weigthed
       * @param nrSpace number of unknowns per node
       */
      static void restriction( const FullGridD* fgFine ,
                               const std::vector<double>& vectFine ,
                               double coefFine ,
                               const FullGridD* fgCoarse ,
                               std::vector<double>& vectCoarse ,
                               double coefCoarse ,
                               int nrSpace) ;

  };

}

#endif /* PROLONGATIONRESTRICTION_HPP_ */
