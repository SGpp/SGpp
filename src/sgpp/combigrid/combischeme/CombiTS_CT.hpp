/* ****************************************************************************
* Copyright (C) 2011 Technische Universitaet Muenchen                         *
* This file is part of the SG++ project. For conditions of distribution and   *
* use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
**************************************************************************** */
// @author Janos Benk (benk@in.tum.de)

#ifndef COMBITS_CT_HPP_
#define COMBITS_CT_HPP_

#include "combigrid/utils/combigrid_ultils.hpp"
#include "CombiSchemeBasis.hpp"

namespace combigrid {

  /** class of the two scale combination scheme (square root CT) <br>*/
  class TS_CT : public CombiSchemeBasis {

    public:

      /** Ctor
       * @param dim dimension of the scheme
       * @param level global level */
      TS_CT( int dim , int level );

      /** Ctor
       * @param dim
       * @param levels the level vector for the dimension adaptive case */
      TS_CT( int dim , const std::vector<int>& levels );

      /** Ctor for cases when in specific dimensions no combi should be done
       * @param dim
       * @param levels the level vector for the dimension adaptive case
       * @param makeCombiInDimension */
      TS_CT( int dim , const std::vector<int>& levels ,
             const std::vector<bool>& makeCombiInDimension );


      /** Ctor for manual steared TS scheme where the user specifies the higher and the lower levels
       * @param minlevels the min levels
       * @param maxlevels the max levels  */
      TS_CT( const std::vector<int>& minlevels ,
             const std::vector<int>& maxlevels  );

    private:

  };
}

#endif /* COMBITS_CT_HPP_ */
