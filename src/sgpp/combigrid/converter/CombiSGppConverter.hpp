/* ****************************************************************************
* Copyright (C) 2011 Technische Universitaet Muenchen                         *
* This file is part of the SG++ project. For conditions of distribution and   *
* use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
**************************************************************************** */
// @author Janos Benk (benk@in.tum.de)

#ifndef COMBISGPPCONVERTER_HPP_
#define COMBISGPPCONVERTER_HPP_

// ---- the combi includes --------
#include "combigrid/utils/combigrid_ultils.hpp"
#include "combigrid/fullgrid/CombiFullGrid.hpp"
#include "combigrid/combigridkernel/CombiGridKernel.hpp"
#include "combigrid/combigrid/AbstractCombiGrid.hpp"

// ---- the SGpp includes --------
#include "base/datatypes/DataVector.hpp"
#include "base/grid/Grid.hpp"
#include "base/grid/GridStorage.hpp"

namespace combigrid {

  /** calss which converts a combi grid to a SGpp sparse grid and back <br>
   * */
  class CombiSGppConverter {

    public:

      /** empty Ctor*/
      CombiSGppConverter() {
        ;
      }

      /** creates a SGpp grid out of a combi grid
       * @param storage the grid storage
       * @param combikernel contains all the information to create a SGpp grid*/
      static void createSGpp( sg::base::GridStorage* storage , const CombiGridKernelD* combikernel );

      /** add the value of the FG to the SGpp grid. We assume that the array has been already resized
       * @param fg the full grid
       * @param coef the multiplicator coefficient (the combination coefficient)*
       * @param storage SGpp grid storage
       * @param alpha coefficient vector */
      static void FullGridToSGpp(const FullGridD* fg , double coef , sg::base::GridStorage* storage , sg::base::DataVector* alpha);

      /** similar as the previous method, but it updates the min and max vector <br>
       * these vectors can be used to estimate the goodnes of the combination of the FGs */
      static void FullGridToSGpp(const FullGridD* fg , double coef , sg::base::GridStorage* storage ,
                                 sg::base::DataVector* alpha , sg::base::DataVector* minAlpha , sg::base::DataVector* maxAlpha );

      /** Sets the values of the FG from the SGpp grid values, this method is used by the decomposition
       * of the SGpp grid into the combi grid
       * @param storage SGpp grid storage
       * @param alpha coefficient vector of the SGpp
       * @param fg the full grid whos values will be set */
      static void SGppToFullGrid( sg::base::GridStorage* storage , sg::base::DataVector* alpha , FullGridD* fg );

  };
}


#endif /* COMBISGPPCONVERTER_HPP_ */
