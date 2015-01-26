/* ****************************************************************************
* Copyright (C) 2011 Technische Universitaet Muenchen                         *
* This file is part of the SG++ project. For conditions of distribution and   *
* use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
**************************************************************************** */
// @author Janos Benk (benk@in.tum.de)

#ifndef COMBIS_CT_HPP_
#define COMBIS_CT_HPP_

#include "combigrid/utils/combigrid_ultils.hpp"
#include "CombiSchemeBasis.hpp"

namespace combigrid {

  /** class to model the standard combination technique. <br>
   * For grids with boundary there will be only the trapezoid grid.
   * The linear grid with boundary has too many boundary points
   * and is not often used in practical applications.
   *  */
  class S_CT : public CombiSchemeBasis {

    public:


      /** simple Ctor with the same level */
      S_CT( int dim , int level );

      /** simple Ctor with the same level and level of truncation*/
      S_CT( int dim , int level , int level_truncation);

      /** S-CT Ctor with dimension adaptive levels */
      S_CT( int dim , const std::vector< int >& levels );

      /** S-CT Ctor with dimension homogenious levels but with dimension adaptive truncation */
      S_CT( int dim , int level , const std::vector< int >& levels_trunk ) ;

      /** S-CT Ctor with dimension adaptive levels and with dimension adaptive levels and truncation */
      S_CT( int dim , const std::vector< int >& levels , const std::vector< int >& levels_trunk );


    protected:


      /** the init function which is general for all types of combi schemes
       * */
      void init( std::vector< int >& levels , std::vector< double >& ratio_ , std::vector< int >& l_user_ );


      /** recursive function for the combination scheme
       * */
      void getTrapezoidsums(std::vector<int>& v , size_t dim , int sum , std::vector< double >& ratio_ ,  std::vector< int >& l_user_);
  };
}


#endif /* COMBIS_CT_HPP_ */
