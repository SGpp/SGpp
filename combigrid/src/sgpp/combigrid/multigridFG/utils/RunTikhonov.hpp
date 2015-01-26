/* ****************************************************************************
* Copyright (C) 2011 Technische Universitaet Muenchen                         *
* This file is part of the SG++ project. For conditions of distribution and   *
* use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
**************************************************************************** */
// @author Janos Benk (benk@in.tum.de)

#ifndef RUNTIKHONOV_HPP_
#define RUNTIKHONOV_HPP_

#include "combigrid.hpp"
#include <sgpp/combigrid/utils/combigrid_ultils.hpp>
#include <sgpp/combigrid/combigrid/AbstractCombiGrid.hpp>
#include <sgpp/combigrid/domain/CombiGridDomain.hpp>

namespace combigrid {

  /** class to run the test cases for tikhonov regularization <br>
   * The input */
  class RunTikhonov {
    public:

      /** empty Ctor*/
      RunTikhonov() {
        ;
      }

      /** empty Dtor*/
      virtual ~RunTikhonov() {
        ;
      }

      /** static function to run the Tikhonov reguralization on a full grid
       * @param domain the domain for the full grid on which the regularization will be done
       * @param levels level vector
       * @param lambda the lambda factor
       * @param Xfile  file path with the X coordinates
       * @param YFile  file with the Y values
       * @return a full grid with the solution */
      static FullGridD* computeFGTikhonov(
        GridDomain& domain,
        const std::vector<int>& levels,
        double lambda,
        const std::string& Xfile ,
        const std::string& YFile );

      /** static function to run the Tikhonov reguralization on a full grid, which grid the user already created
       * @param fg [IN/OUT] input full grid , must have a valid domain
       * @param unknowns [IN/OUT] vector of unknowns
       * @param lambda [IN] the lambda factor
       * @param XCoords  [IN] X points/coordinates
       * @param YPoint  [IN] Y values */
      static void computeFGTikhonov_FG(
        FullGridD* fg,
        std::vector<double>& unknowns ,
        double lambda,
        std::vector<double>& XCoords,
        std::vector<double>& YPoint);

      /** function to compute the Tikhonov problem on a combination grid
       * @param combiG [IN/OUT] input combination grid, each full grid should have a valid domain
       * @param lambda [IN] the lambda factor
       * @param XCoords  [IN] X points/coordinates
       * @param YPoint  [IN] Y values */
      static void computeTikhonov_CG(
        AbstractCombiGrid* combiG,
        double lambda,
        std::vector<double>& XCoords,
        std::vector<double>& YPoint);

      /** function to compute the Tikhonov problem on a full grid, and finds automatically the
       * optimal lambda coefficient
       * @param fg [IN/OUT] input combination grid, each full grid should have a valid domain
       * @param lambda [OUT] the lambda factor
       * @param XCoords  [IN] X points/coordinates
       * @param YPoint  [IN] Y values */
      static void computeTikhonov_FG_crossvalidation(
        FullGridD* fg,
        double& lambda,
        std::vector<double>& XCoords,
        std::vector<double>& YPoint);

      /** function to compute the Tikhonov problem on a combination grid
       * @param combiG [IN/OUT] input combination grid, each full grid should have a valid domain
       * @param lambda [OUT] the lambda factor
       * @param XCoords  [IN] X points/coordinates
       * @param YPoint  [IN] Y values */
      static void computeTikhonov_CG_crossvalidation(
        AbstractCombiGrid* combiG,
        double& lambda,
        std::vector<double>& XCoords,
        std::vector<double>& YPoint);


      /** read in the input values for the regularization*/
      static void readInInput(const std::string& XfileS , const std::string& YFileS ,
                              int& dimensions ,  int& nrPoints ,
                              std::vector<double>& XCoords, std::vector<double>& YPoint);

    private :

      /** measure error between the full grid and the points*/
      static double measureErrorFG(FullGridD* fg,
                                   int dim ,
                                   std::vector<double>& XCoords,
                                   std::vector<double>& YPoint);

      /** measure error between the combination grid and the points*/
      static double measureErrorCG(AbstractCombiGrid* combiG,
                                   int dim ,
                                   std::vector<double>& XCoords,
                                   std::vector<double>& YPoint);
  };

}

#endif /* RUNTIKHONOV_HPP_ */
