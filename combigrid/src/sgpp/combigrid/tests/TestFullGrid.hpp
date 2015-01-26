// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at 
// sgpp.sparsegrids.org

#ifndef TESTFULLGRID_HPP_
#define TESTFULLGRID_HPP_

#include "combigrid.hpp"

namespace combigrid {

  /** Class to thes the all the full grids functionalities */
  class TestFullGrid {
    public:

      TestFullGrid();

      static void test_all_cases() {
        COMBIGRID_OUT_LEVEL3( 4 , "Testing Full SGPP::base::Grid ... ");
        // test here all the functionalities of the full grid
        test1();

        // 3D test
        test2();
      }

      /** test function for the grid, this should be exactly interpolated */
      static double testfunction2D( double x1, double x2) {
        return (2 + 3.0 * x1 + 0.5 * x2);
      }

      static double testfunction3D( double x1, double x2, double x3) {
        return (2 + 3.0 * x1 + 0.5 * x2 + 1.5 * x3);
      }

      static void test1() {
        // test the full grid
        FullGrid<double>* fg = new FullGrid<double>(2, 3);
        FullGrid<double>* fg_wg = new FullGrid<double>(2, 3, false);

        fg->createFullGrid();
        fg_wg->createFullGrid();

        std::vector<double> coords(2, 0.0);

        // ---
        for (int i = 0 ; i < fg->getNrElements() ; i++) {
          fg->getCoords( i , coords );
          fg->getElementVector()[i] = testfunction2D( coords[0] , coords[1] );
        }

        // ---
        for (int i = 0 ; i < fg_wg->getNrElements() ; i++) {
          fg_wg->getCoords( i , coords );
          fg_wg->getElementVector()[i] = testfunction2D( coords[0] , coords[1] );
        }

        coords[0] = 1.0 / 3.0;
        coords[1] = 1.0 / 3.0;
        COMBIGRID_ERROR_TEST_EQUAL( fg->eval(coords) , testfunction2D( coords[0] , coords[1] ) , 1e-10 , "" );
        COMBIGRID_ERROR_TEST_EQUAL( fg_wg->eval(coords) , testfunction2D( coords[0] , coords[1] ) , 1e-10 , "" );

        coords[0] = 0.0001;
        coords[1] = 1.0 / 3.0;
        COMBIGRID_ERROR_TEST_EQUAL( fg->eval(coords) , testfunction2D( coords[0] , coords[1] ) , 1e-10 , "" );
        COMBIGRID_ERROR_TEST_EQUAL( fg_wg->eval(coords) , testfunction2D( coords[0] , coords[1] ) , 1e-10 , "" );

        coords[0] = 1.0 / 3.0;
        coords[1] = 0.0001;
        COMBIGRID_ERROR_TEST_EQUAL( fg->eval(coords) , testfunction2D( coords[0] , coords[1] ) , 1e-10 , "" );
        COMBIGRID_ERROR_TEST_EQUAL( fg_wg->eval(coords) , testfunction2D( coords[0] , coords[1] ) , 1e-10 , "" );

        coords[0] = 0.00001;
        coords[1] = 0.00001;
        COMBIGRID_ERROR_TEST_EQUAL( fg->eval(coords) , testfunction2D( coords[0] , coords[1] ) , 1e-10 , "" );
        COMBIGRID_ERROR_TEST_EQUAL( fg_wg->eval(coords) , testfunction2D( coords[0] , coords[1] ) , 1e-10 , "" );

        coords[0] = 1.0;
        coords[1] = 1.0;
        COMBIGRID_ERROR_TEST_EQUAL( fg->eval(coords) , testfunction2D( coords[0] , coords[1] ) , 1e-10 , "" );
        COMBIGRID_ERROR_TEST_EQUAL( fg_wg->eval(coords) , testfunction2D( coords[0] , coords[1] ) , 1e-10 , "" );

        coords[0] = 0.999;
        coords[1] = 0.95;
        COMBIGRID_ERROR_TEST_EQUAL( fg->eval(coords) , testfunction2D( coords[0] , coords[1] ) , 1e-10 , "" );
        COMBIGRID_ERROR_TEST_EQUAL( fg_wg->eval(coords) , testfunction2D( coords[0] , coords[1] ) , 1e-10 , "" );

        coords[0] = 0.999;
        coords[1] = 0.001;
        COMBIGRID_ERROR_TEST_EQUAL( fg->eval(coords) , testfunction2D( coords[0] , coords[1] ) , 1e-10 , "" );
        COMBIGRID_ERROR_TEST_EQUAL( fg_wg->eval(coords) , testfunction2D( coords[0] , coords[1] ) , 1e-10 , "" );

        // delete the full grid
        delete fg;
        delete fg_wg;
      }


      static void test2() {
        // test the full grid

        std::vector<int> levels(3, 5);
        std::vector<bool> hasBdrPoints(3, true);
        levels[1] = 6;
        levels[2] = 4;
        hasBdrPoints[1] = false;
        FullGrid<double>* fg = new FullGrid<double>( 3 , levels , hasBdrPoints );
        levels[1] = 4;
        levels[2] = 6;
        hasBdrPoints[1] = true;
        hasBdrPoints[0] = false;
        hasBdrPoints[2] = false;
        FullGrid<double>* fg_wg = new FullGrid<double>( 3 , levels , hasBdrPoints );

        fg->createFullGrid();
        fg_wg->createFullGrid();

        std::vector<double> coords(3, 0.0);

        // ---
        for (int i = 0 ; i < fg->getNrElements() ; i++) {
          fg->getCoords( i , coords );
          fg->getElementVector()[i] = testfunction3D( coords[0] , coords[1] , coords[2] );
        }

        // ---
        //COMBIGRID_OUT_LEVEL3( 4 , " Set values for FG without boundary ");
        for (int i = 0 ; i < fg_wg->getNrElements() ; i++) {
          fg_wg->getCoords( i , coords );
          fg_wg->getElementVector()[i] = testfunction3D( coords[0] , coords[1] , coords[2] );
        }

        coords[0] = 1.0 / 3.0;
        coords[1] = 1.0 / 3.0;
        coords[2] = 1.0 / 3.0;
        COMBIGRID_ERROR_TEST_EQUAL( fg->eval(coords) , testfunction3D( coords[0] , coords[1] , coords[2] ) , 1e-10 , "" );
        COMBIGRID_ERROR_TEST_EQUAL( fg_wg->eval(coords) , testfunction3D( coords[0] , coords[1] , coords[2] ) , 1e-10 , "" );

        coords[0] = 0.0001;
        coords[1] = 1.0 / 3.0;
        coords[2] = 1.0 / 3.0;
        COMBIGRID_ERROR_TEST_EQUAL( fg->eval(coords) , testfunction3D( coords[0] , coords[1] , coords[2] ) , 1e-10 , "" );
        COMBIGRID_ERROR_TEST_EQUAL( fg_wg->eval(coords) , testfunction3D( coords[0] , coords[1] , coords[2] ) , 1e-10 , "" );

        coords[0] = 1.0 / 3.0;
        coords[1] = 0.0001;
        coords[2] = 0.0001;
        COMBIGRID_ERROR_TEST_EQUAL( fg->eval(coords) , testfunction3D( coords[0] , coords[1] , coords[2] ) , 1e-10 , "" );
        COMBIGRID_ERROR_TEST_EQUAL( fg_wg->eval(coords) , testfunction3D( coords[0] , coords[1] , coords[2] ) , 1e-10 , "" );

        coords[0] = 0.0;
        coords[1] = 0.0;
        coords[2] = 0.0;
        COMBIGRID_ERROR_TEST_EQUAL( fg->eval(coords) , testfunction3D( coords[0] , coords[1] , coords[2] ) , 1e-10 , "" );
        COMBIGRID_ERROR_TEST_EQUAL( fg_wg->eval(coords) , testfunction3D( coords[0] , coords[1] , coords[2] ) , 1e-10 , "" );

        coords[0] = 1.0;
        coords[1] = 1.0;
        coords[2] = 1.0;
        COMBIGRID_ERROR_TEST_EQUAL( fg->eval(coords) , testfunction3D( coords[0] , coords[1] , coords[2] ) , 1e-10 , "" );
        COMBIGRID_ERROR_TEST_EQUAL( fg_wg->eval(coords) , testfunction3D( coords[0] , coords[1] , coords[2] ) , 1e-10 , "" );

        coords[0] = 0.999;
        coords[1] = 0.95;
        coords[2] = 0.95;
        COMBIGRID_ERROR_TEST_EQUAL( fg->eval(coords) , testfunction3D( coords[0] , coords[1] , coords[2] ) , 1e-10 , "" );
        COMBIGRID_ERROR_TEST_EQUAL( fg_wg->eval(coords) , testfunction3D( coords[0] , coords[1] , coords[2] ) , 1e-10 , "" );

        coords[0] = 0.999;
        coords[1] = 0.001;
        coords[2] = 0.999;
        COMBIGRID_ERROR_TEST_EQUAL( fg->eval(coords) , testfunction3D( coords[0] , coords[1] , coords[2] ) , 1e-10 , "" );
        COMBIGRID_ERROR_TEST_EQUAL( fg_wg->eval(coords) , testfunction3D( coords[0] , coords[1] , coords[2] ) , 1e-10 , "" );

        // delete the full grid
        delete fg;
        delete fg_wg;
      }

  };
}

#endif /* TESTFULLGRID_HPP_ */