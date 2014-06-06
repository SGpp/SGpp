/* ****************************************************************************
* Copyright (C) 2011 Technische Universitaet Muenchen                         *
* This file is part of the SG++ project. For conditions of distribution and   *
* use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
**************************************************************************** */
// @author Janos Benk (benk@in.tum.de)

#ifndef TESTSERIALS_CTGRID_HPP_
#define TESTSERIALS_CTGRID_HPP_

#include "combigrid.hpp"

namespace combigrid {

  /** Class to thes the all the full grids functionalities */
  class TestSerialS_CTGrid {
    public:

      /** empty Ctor */
      TestSerialS_CTGrid();

      static void test_all_cases() {
        COMBIGRID_OUT_LEVEL3( 4 , "TestSerialS_CTGrid ... ");
        // 2D
        test1();
        test12();

        // 3D test
        test2();
        test21();
      }

      /** test function for the grid, this should be exactly interpolated */
      static double testfunction2D( double x1, double x2) {
        return (2 + 3.0 * x1 + 0.5 * x2);
      }

      static double testfunction2D1( double x1, double x2) {
        return ::exp(0.25 * x1 + 0.5 * x2);
      }

      static double testfunction3D( double x1, double x2, double x3) {
        return (2 + 3.0 * x1 + 0.5 * x2 + 1.5 * x3);
      }

      static double testfunction3D1( double x1, double x2, double x3) {
        return (0.01 + 1.0 * x1 + 0.5 * x2 + 1.5 * x3) * (0.01 + 2.0 * x1 + 0.5 * x2 + 1.5 * x3);
      }

      static void test1() {
        std::vector<double> coords( 2 , 0.0);
        CombiSchemeBasis* combischeme = new S_CT( 2 , 6 );
        AbstractCombiGrid* combigrid = new SerialCombiGrid(combischeme);
        AbstractCombiGrid* combigrid1 = new SerialCombiGrid(combischeme, false );

        // creates the full grids in the combination grid
        combigrid->createFullGrids();
        combigrid1->createFullGrids();

        // for each full grid set the values
        for ( int fg = 0 ; fg < combigrid->getNrFullGrid() ; fg++) {
          FullGridD* fgp = combigrid->getFullGrid(fg);

          // for this full grid set the values
          for (int i = 0 ; i < fgp->getNrElements() ; i++) {
            fgp->getCoords( i , coords );
            fgp->getElementVector()[i] = testfunction2D( coords[0] , coords[1] );
          }
        }

        for ( int fg = 0 ; fg < combigrid1->getNrFullGrid() ; fg++) {
          FullGridD* fgp = combigrid1->getFullGrid(fg);

          // for this full grid set the values
          for (int i = 0 ; i < fgp->getNrElements() ; i++) {
            fgp->getCoords( i , coords );
            fgp->getElementVector()[i] = testfunction2D( coords[0] , coords[1] );
          }
        }

        coords[0] = 1.0 / 3.0;
        coords[1] = 1.0 / 3.0;
        COMBIGRID_ERROR_TEST_EQUAL( combigrid->eval(coords) , testfunction2D( coords[0] , coords[1] ) , 1e-10 , "" );
        COMBIGRID_ERROR_TEST_EQUAL( combigrid1->eval(coords) , testfunction2D( coords[0] , coords[1] ) , 1e-10 , "" );

        coords[0] = 0.0001;
        coords[1] = 1.0 / 3.0;
        COMBIGRID_ERROR_TEST_EQUAL( combigrid->eval(coords) , testfunction2D( coords[0] , coords[1] ) , 1e-10 , "" );
        COMBIGRID_ERROR_TEST_EQUAL( combigrid1->eval(coords) , testfunction2D( coords[0] , coords[1] ) , 1e-6 , "" );

        coords[0] = 1.0 / 3.0;
        coords[1] = 0.0001;
        COMBIGRID_ERROR_TEST_EQUAL( combigrid->eval(coords) , testfunction2D( coords[0] , coords[1] ) , 1e-10 , "" );
        COMBIGRID_ERROR_TEST_EQUAL( combigrid1->eval(coords) , testfunction2D( coords[0] , coords[1] ) , 1e-10 , "" );

        coords[0] = 0.00001;
        coords[1] = 0.00001;
        COMBIGRID_ERROR_TEST_EQUAL( combigrid->eval(coords) , testfunction2D( coords[0] , coords[1] ) , 1e-10 , "" );
        COMBIGRID_ERROR_TEST_EQUAL( combigrid1->eval(coords) , testfunction2D( coords[0] , coords[1] ) , 1e-10 , "" );

        coords[0] = 1.0;
        coords[1] = 1.0;
        COMBIGRID_ERROR_TEST_EQUAL( combigrid->eval(coords) , testfunction2D( coords[0] , coords[1] ) , 1e-10 , "" );
        COMBIGRID_ERROR_TEST_EQUAL( combigrid1->eval(coords) , testfunction2D( coords[0] , coords[1] ) , 1e-10 , "" );

        coords[0] = 0.999;
        coords[1] = 0.95;
        COMBIGRID_ERROR_TEST_EQUAL( combigrid->eval(coords) , testfunction2D( coords[0] , coords[1] ) , 1e-10 , "" );
        COMBIGRID_ERROR_TEST_EQUAL( combigrid1->eval(coords) , testfunction2D( coords[0] , coords[1] ) , 1e-10 , "" );

        coords[0] = 0.999;
        coords[1] = 0.001;
        COMBIGRID_ERROR_TEST_EQUAL( combigrid->eval(coords) , testfunction2D( coords[0] , coords[1] ) , 1e-10 , "" );
        COMBIGRID_ERROR_TEST_EQUAL( combigrid1->eval(coords) , testfunction2D( coords[0] , coords[1] ) , 1e-10 , "" );

        delete combigrid;
        delete combigrid1;
        delete combischeme;
      }


      static void test12() {
        std::vector<double> coords( 2 , 0.0);
        CombiSchemeBasis* combischeme = new S_CT( 2 , 9 , 3);
        AbstractCombiGrid* combigrid = new SerialCombiGrid(combischeme);
        AbstractCombiGrid* combigrid1 = new SerialCombiGrid(combischeme, false );

        // creates the full grids in the combination grid
        combigrid->createFullGrids();
        combigrid1->createFullGrids();

        // for each full grid set the values
        for ( int fg = 0 ; fg < combigrid->getNrFullGrid() ; fg++) {
          FullGridD* fgp = combigrid->getFullGrid(fg);

          // for this full grid set the values
          for (int i = 0 ; i < fgp->getNrElements() ; i++) {
            fgp->getCoords( i , coords );
            fgp->getElementVector()[i] = testfunction2D1( coords[0] , coords[1] );
          }
        }

        for ( int fg = 0 ; fg < combigrid1->getNrFullGrid() ; fg++) {
          FullGridD* fgp = combigrid1->getFullGrid(fg);

          // for this full grid set the values
          for (int i = 0 ; i < fgp->getNrElements() ; i++) {
            fgp->getCoords( i , coords );
            fgp->getElementVector()[i] = testfunction2D1( coords[0] , coords[1] );
          }
        }

        //COMBIGRID_ERROR_TEST_EQUAL( combigrid->getNrFullGrid() , 3.0 , 1e-10 , "" );
        //COMBIGRID_ERROR_TEST_EQUAL( combigrid1->getNrFullGrid() , 3.0 , 1e-10 , "" );

        coords[0] = 1.0 / 3.0;
        coords[1] = 1.0 / 3.0;
        COMBIGRID_ERROR_TEST_EQUAL( combigrid->eval(coords) , testfunction2D1( coords[0] , coords[1] ) , 1e-5 , "" );
        COMBIGRID_ERROR_TEST_EQUAL( combigrid1->eval(coords) , testfunction2D1( coords[0] , coords[1] ) , 1e-5 , "" );

        coords[0] = 0.0001;
        coords[1] = 1.0 / 3.0;
        COMBIGRID_ERROR_TEST_EQUAL( combigrid->eval(coords) , testfunction2D1( coords[0] , coords[1] ) , 1e-5 , "" );
        COMBIGRID_ERROR_TEST_EQUAL( combigrid1->eval(coords) , testfunction2D1( coords[0] , coords[1] ) , 1e-5 , "" );

        coords[0] = 1.0 / 3.0;
        coords[1] = 0.0001;
        COMBIGRID_ERROR_TEST_EQUAL( combigrid->eval(coords) , testfunction2D1( coords[0] , coords[1] ) , 1e-5 , "" );
        COMBIGRID_ERROR_TEST_EQUAL( combigrid1->eval(coords) , testfunction2D1( coords[0] , coords[1] ) , 1e-5 , "" );

        coords[0] = 0.00001;
        coords[1] = 0.00001;
        COMBIGRID_ERROR_TEST_EQUAL( combigrid->eval(coords) , testfunction2D1( coords[0] , coords[1] ) , 1e-5 , "" );
        COMBIGRID_ERROR_TEST_EQUAL( combigrid1->eval(coords) , testfunction2D1( coords[0] , coords[1] ) , 1e-5 , "" );

        coords[0] = 1.0;
        coords[1] = 1.0;
        COMBIGRID_ERROR_TEST_EQUAL( combigrid->eval(coords) , testfunction2D1( coords[0] , coords[1] ) , 1e-5 , "" );
        COMBIGRID_ERROR_TEST_EQUAL( combigrid1->eval(coords) , testfunction2D1( coords[0] , coords[1] ) , 1e-5 , "" );

        coords[0] = 0.999;
        coords[1] = 0.95;
        COMBIGRID_ERROR_TEST_EQUAL( combigrid->eval(coords) , testfunction2D1( coords[0] , coords[1] ) , 1e-5 , "" );
        COMBIGRID_ERROR_TEST_EQUAL( combigrid1->eval(coords) , testfunction2D1( coords[0] , coords[1] ) , 1e-5 , "" );

        coords[0] = 0.999;
        coords[1] = 0.001;
        COMBIGRID_ERROR_TEST_EQUAL( combigrid->eval(coords) , testfunction2D1( coords[0] , coords[1] ) , 1e-5 , "" );
        COMBIGRID_ERROR_TEST_EQUAL( combigrid1->eval(coords) , testfunction2D1( coords[0] , coords[1] ) , 1e-5 , "" );

        delete combigrid;
        delete combigrid1;
        delete combischeme;
      }



      static void test2() {

        std::vector<double> coords( 3 , 0.0);
        std::vector<int> levels( 3 , 6);
        CombiSchemeBasis* combischeme = new S_CT( 3 , 6 , 1);
        AbstractCombiGrid* combigrid = new SerialCombiGrid(combischeme);
        AbstractCombiGrid* combigrid1 = new SerialCombiGrid(combischeme, false );

        // creates the full grids in the combination grid
        combigrid->createFullGrids();
        combigrid1->createFullGrids();

        // for each full grid set the values
        for ( int fg = 0 ; fg < combigrid->getNrFullGrid() ; fg++) {
          FullGridD* fgp = combigrid->getFullGrid(fg);

          // for this full grid set the values
          for (int i = 0 ; i < fgp->getNrElements() ; i++) {
            fgp->getCoords( i , coords );
            fgp->getElementVector()[i] = testfunction3D( coords[0] , coords[1] , coords[2]);
          }
        }

        for ( int fg = 0 ; fg < combigrid1->getNrFullGrid() ; fg++) {
          FullGridD* fgp = combigrid1->getFullGrid(fg);

          // for this full grid set the values
          for (int i = 0 ; i < fgp->getNrElements() ; i++) {
            fgp->getCoords( i , coords );
            fgp->getElementVector()[i] = testfunction3D( coords[0] , coords[1] , coords[2]);
          }
        }

        coords[0] = 1.0 / 3.0;
        coords[1] = 1.0 / 3.0;
        coords[2] = 1.0 / 3.0;
        COMBIGRID_ERROR_TEST_EQUAL( combigrid->eval(coords) , testfunction3D( coords[0] , coords[1] , coords[2] ) , 1e-10 , "" );
        COMBIGRID_ERROR_TEST_EQUAL( combigrid1->eval(coords) , testfunction3D( coords[0] , coords[1] , coords[2] ) , 1e-10 , "" );

        coords[0] = 0.0001;
        coords[1] = 1.0 / 3.0;
        coords[2] = 1.0 / 3.0;
        COMBIGRID_ERROR_TEST_EQUAL( combigrid->eval(coords) , testfunction3D( coords[0] , coords[1] , coords[2] ) , 1e-10 , "" );
        COMBIGRID_ERROR_TEST_EQUAL( combigrid1->eval(coords) , testfunction3D( coords[0] , coords[1] , coords[2] ) , 1e-10 , "" );

        coords[0] = 1.0 / 3.0;
        coords[1] = 0.0001;
        coords[2] = 0.0001;
        COMBIGRID_ERROR_TEST_EQUAL( combigrid->eval(coords) , testfunction3D( coords[0] , coords[1] , coords[2] ) , 1e-10 , "" );
        COMBIGRID_ERROR_TEST_EQUAL( combigrid1->eval(coords) , testfunction3D( coords[0] , coords[1] , coords[2] ) , 1e-10 , "" );

        coords[0] = 0.0;
        coords[1] = 0.0;
        coords[2] = 0.0;
        COMBIGRID_ERROR_TEST_EQUAL( combigrid->eval(coords) , testfunction3D( coords[0] , coords[1] , coords[2] ) , 1e-10 , "" );
        COMBIGRID_ERROR_TEST_EQUAL( combigrid1->eval(coords) , testfunction3D( coords[0] , coords[1] , coords[2] ) , 1e-10 , "" );

        coords[0] = 1.0;
        coords[1] = 1.0;
        coords[2] = 1.0;
        COMBIGRID_ERROR_TEST_EQUAL( combigrid->eval(coords) , testfunction3D( coords[0] , coords[1] , coords[2] ) , 1e-10 , "" );
        COMBIGRID_ERROR_TEST_EQUAL( combigrid1->eval(coords) , testfunction3D( coords[0] , coords[1] , coords[2] ) , 1e-10 , "" );

        coords[0] = 0.999;
        coords[1] = 0.95;
        coords[2] = 0.95;
        COMBIGRID_ERROR_TEST_EQUAL( combigrid->eval(coords) , testfunction3D( coords[0] , coords[1] , coords[2] ) , 1e-10 , "" );
        COMBIGRID_ERROR_TEST_EQUAL( combigrid1->eval(coords) , testfunction3D( coords[0] , coords[1] , coords[2] ) , 1e-10 , "" );

        coords[0] = 0.999;
        coords[1] = 0.001;
        coords[2] = 0.999;
        COMBIGRID_ERROR_TEST_EQUAL( combigrid->eval(coords) , testfunction3D( coords[0] , coords[1] , coords[2] ) , 1e-10 , "" );
        COMBIGRID_ERROR_TEST_EQUAL( combigrid1->eval(coords) , testfunction3D( coords[0] , coords[1] , coords[2] ) , 1e-10 , "" );

        delete combigrid;
        delete combigrid1;
        delete combischeme;
      }



      static void test21() {

        std::vector<double> coords( 3 , 0.0);
        std::vector<int> levels( 3 , 9);
        std::vector<int> levels_trunc( 3 , 3);
        levels[0] = 10;
        levels[1] = 9;
        levels[2] = 10;
        levels_trunc[0] = 1;
        levels_trunc[1] = 3;
        levels_trunc[2] = 1;
        CombiSchemeBasis* combischeme = new S_CT( 3 , levels , levels_trunc );
        AbstractCombiGrid* combigrid = new SerialCombiGrid(combischeme);
        AbstractCombiGrid* combigrid1 = new SerialCombiGrid(combischeme, false );

        // creates the full grids in the combination grid
        combigrid->createFullGrids();
        combigrid1->createFullGrids();

        // for each full grid set the values
        for ( int fg = 0 ; fg < combigrid->getNrFullGrid() ; fg++) {
          FullGridD* fgp = combigrid->getFullGrid(fg);

          // for this full grid set the values
          for (int i = 0 ; i < fgp->getNrElements() ; i++) {
            fgp->getCoords( i , coords );
            fgp->getElementVector()[i] = testfunction3D1( coords[0] , coords[1] , coords[2]);
          }
        }

        for ( int fg = 0 ; fg < combigrid1->getNrFullGrid() ; fg++) {
          FullGridD* fgp = combigrid1->getFullGrid(fg);

          // for this full grid set the values
          for (int i = 0 ; i < fgp->getNrElements() ; i++) {
            fgp->getCoords( i , coords );
            fgp->getElementVector()[i] = testfunction3D1( coords[0] , coords[1] , coords[2]);
          }
        }

        coords[0] = 1.0 / 3.0;
        coords[1] = 1.0 / 3.0;
        coords[2] = 1.0 / 3.0;
        COMBIGRID_ERROR_TEST_EQUAL( combigrid->eval(coords) , testfunction3D1( coords[0] , coords[1] , coords[2] ) , 1e-4 , "" );
        COMBIGRID_ERROR_TEST_EQUAL( combigrid1->eval(coords) , testfunction3D1( coords[0] , coords[1] , coords[2] ) , 5e-4 , "" );

        coords[0] = 0.0001;
        coords[1] = 1.0 / 3.0;
        coords[2] = 1.0 / 3.0;
        COMBIGRID_ERROR_TEST_EQUAL( combigrid->eval(coords) , testfunction3D1( coords[0] , coords[1] , coords[2] ) , 1e-4 , "" );
        COMBIGRID_ERROR_TEST_EQUAL( combigrid1->eval(coords) , testfunction3D1( coords[0] , coords[1] , coords[2] ) , 1e-3 , "" );

        coords[0] = 1.0 / 3.0;
        coords[1] = 0.0001;
        coords[2] = 0.0001;
        COMBIGRID_ERROR_TEST_EQUAL( combigrid->eval(coords) , testfunction3D1( coords[0] , coords[1] , coords[2] ) , 1e-4 , "" );
        COMBIGRID_ERROR_TEST_EQUAL( combigrid1->eval(coords) , testfunction3D1( coords[0] , coords[1] , coords[2] ) , 1e-3 , "" );

        coords[0] = 0.0;
        coords[1] = 0.0;
        coords[2] = 0.0;
        COMBIGRID_ERROR_TEST_EQUAL( combigrid->eval(coords) , testfunction3D1( coords[0] , coords[1] , coords[2] ) , 1e-4 , "" );
        COMBIGRID_ERROR_TEST_EQUAL( combigrid1->eval(coords) , testfunction3D1( coords[0] , coords[1] , coords[2] ) , 1e-3 , "" );

        coords[0] = 1.0;
        coords[1] = 1.0;
        coords[2] = 1.0;
        COMBIGRID_ERROR_TEST_EQUAL( combigrid->eval(coords) , testfunction3D1( coords[0] , coords[1] , coords[2] ) , 1e-4 , "" );
        COMBIGRID_ERROR_TEST_EQUAL( combigrid1->eval(coords) , testfunction3D1( coords[0] , coords[1] , coords[2] ) , 5e-3 , "" );

        coords[0] = 0.999;
        coords[1] = 0.95;
        coords[2] = 0.95;
        COMBIGRID_ERROR_TEST_EQUAL( combigrid->eval(coords) , testfunction3D1( coords[0] , coords[1] , coords[2] ) , 1e-4 , "" );
        COMBIGRID_ERROR_TEST_EQUAL( combigrid1->eval(coords) , testfunction3D1( coords[0] , coords[1] , coords[2] ) , 5e-3 , "" );

        coords[0] = 0.999;
        coords[1] = 0.001;
        coords[2] = 0.999;
        COMBIGRID_ERROR_TEST_EQUAL( combigrid->eval(coords) , testfunction3D1( coords[0] , coords[1] , coords[2] ) , 1e-4 , "" );
        COMBIGRID_ERROR_TEST_EQUAL( combigrid1->eval(coords) , testfunction3D1( coords[0] , coords[1] , coords[2] ) , 5e-3 , "" );

        delete combigrid;
        delete combigrid1;
        delete combischeme;

      }

  };
}

#endif /* TESTSERIALS_CTGRID_HPP_ */
