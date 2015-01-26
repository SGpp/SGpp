/* ****************************************************************************
* Copyright (C) 2011 Technische Universitaet Muenchen                         *
* This file is part of the SG++ project. For conditions of distribution and   *
* use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
**************************************************************************** */
// @author Janos Benk (benk@in.tum.de)

#ifndef RunPoisson_HPP_
#define RunPoisson_HPP_

#include "combigrid.hpp"
#include <sgpp/combigrid/utils/combigrid_ultils.hpp>
#include <sgpp/combigrid/combigrid/AbstractCombiGrid.hpp>
#include <sgpp/combigrid/domain/CombiGridDomain.hpp>
#include <sgpp/combigrid/multigridFG/operators/PoissonOperator.hpp>

namespace combigrid {

  /** class to run the test cases for tikhonov regularization <br>
   * The input */
  class RunPoisson {
    public:

      /** empty Ctor*/
      RunPoisson() {
        ;
      }

      /** empty Dtor*/
      virtual ~RunPoisson() {
        ;
      }

      /** static function to run the Poisson problem on a full grid
       * @param domain the domain for the full grid
       * @param levels  level vector
       * @param sigma  diffusion coefficients in each direction
       * @param startValue the value of all unknowns at the beginning
       * @param callbackRHS  right hand side call back
       * */
      static FullGridD* computeFGPoisson(
        GridDomain& domain,
        const std::vector<int>& levels,
        const std::vector<double>& sigma ,
        double startValue ,
        const CallBackRHS* callbackRHS );

      /** static function to run the Poisson problem on a full grid, which the user already created
       * @param fg  the full grid
       * @param sigma  diffusion coefficients in each direction
       * @param unknowns  the initial value of unkowns
       * @param callbackRHS right hand side call back
       * */
      static void computeFGPoisson_FG(
        FullGridD* fg ,
        const std::vector<double>& sigma ,
        std::vector<double>& unknowns ,
        const CallBackRHS* callbackRHS );

    private :

  };

}

#endif /* RunPoisson_HPP_ */
