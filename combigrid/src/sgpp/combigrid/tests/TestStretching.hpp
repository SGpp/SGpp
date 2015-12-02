// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#ifndef TESTSTRETCHING_HPP_
#define TESTSTRETCHING_HPP_

#include <combigrid.hpp>

namespace combigrid {

  /** Class to thes the all the full grids functionalities */
  class TestStretching {
    public:

      TestStretching() {
        ;
      }

      static void test_all_cases() {
        COMBIGRID_OUT_LEVEL3( 4 , "Testing SGPP::base::Stretching ... ");
        // test here all the functionalities of the full grid
        // 2D test
        test1();
        test12();
        test13();
        test14();

        // 3D test
        test2();
      }

      /** test function for the grid, this should be exactly interpolated */
      static double testfunction2D( double x1, double x2) {
        //COMBIGRID_OUT_LEVEL3( 4 , "testfunction2D x1=" << x1 << " , x2=" << x2 << " ret: " << (2+3.0*x1+0.5*x2));
        return (2 + 3.0 * x1 + 0.5 * x2);
      }

      static double testfunction3D( double x1, double x2, double x3) {
        return (2.0 + 3.0 * x1 + 0.5 * x2 + 1.5 * x3);
      }

      static void test1() {
        // test the full grid

        std::vector<double> min(2);
        min[0] = -1;
        min[1] = 2;
        std::vector<double> max(2);
        max[0] = 1;
        max[1] = 3;
        std::vector<int> levels(2);
        levels[0] = 6;
        levels[1] = 7;
        //std::vector<int> levels(2); levels[0] = 3; levels[1] = 3;
        double funcVal = 0.0 , gridVal = 0.0;

        AbstractStretchingMaker* stretchingMaker = new UniformStretching();
        GridDomain* domain = new GridDomain( 2 , levels , min , max , (*stretchingMaker) );

        FullGrid<double>* fg = new FullGrid<double>(2, 5);
        FullGrid<double>* fg_wg = new FullGrid<double>(2, levels, false);
        //FullGrid<double>* fg = new FullGrid<double>(2,2);
        //FullGrid<double>* fg_wg = new FullGrid<double>(2,2,false);

        fg->createFullGrid();
        fg_wg->createFullGrid();
        fg->setDomain( domain );
        fg_wg->setDomain( domain );

        std::vector<double> coords(2, 0.0);
        std::vector<double> coords_tmp(2, 0.0);

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

        coords[0] = min[0] + (max[0] - min[0]) * 1.0 / 3.0;
        coords[1] = min[1] + (max[1] - min[1]) * 1.0 / 3.0;
        coords_tmp = coords;
        funcVal = testfunction2D( coords[0] , coords[1] );
        gridVal = fg->eval(coords_tmp);
        COMBIGRID_ERROR_TEST_EQUAL( gridVal , funcVal , 1e-10 , "" );
        coords_tmp = coords;
        gridVal = fg_wg->eval(coords_tmp);
        COMBIGRID_ERROR_TEST_EQUAL( gridVal , funcVal , 1e-10 , "" );
        coords_tmp = coords;

        coords[0] = min[0] + (max[0] - min[0]) * 0.0001;
        coords[1] = min[1] + (max[1] - min[1]) * 1.0 / 3.0;
        coords_tmp = coords;
        funcVal = testfunction2D( coords[0] , coords[1] );
        gridVal = fg->eval(coords_tmp);
        COMBIGRID_ERROR_TEST_EQUAL( gridVal , funcVal , 1e-10 , "" );
        coords_tmp = coords;
        gridVal = fg_wg->eval(coords_tmp);
        COMBIGRID_ERROR_TEST_EQUAL( gridVal , funcVal , 1e-10 , "" );
        coords_tmp = coords;

        coords[0] = min[0] + (max[0] - min[0]) * 1.0 / 3.0;
        coords[1] = min[1] + (max[1] - min[1]) * 0.0001;
        coords_tmp = coords;
        funcVal = testfunction2D( coords[0] , coords[1] );
        gridVal = fg->eval(coords_tmp);
        COMBIGRID_ERROR_TEST_EQUAL( gridVal , funcVal , 1e-10 , "" );
        coords_tmp = coords;
        gridVal = fg_wg->eval(coords_tmp);
        COMBIGRID_ERROR_TEST_EQUAL( gridVal , funcVal , 1e-10 , "" );
        coords_tmp = coords;

        coords[0] = min[0] + (max[0] - min[0]) * 0.00001;
        coords[1] = min[1] + (max[1] - min[1]) * 0.00001;
        coords_tmp = coords;
        funcVal = testfunction2D( coords[0] , coords[1] );
        gridVal = fg->eval(coords_tmp);
        COMBIGRID_ERROR_TEST_EQUAL( gridVal , funcVal , 1e-10 , "" );
        coords_tmp = coords;
        gridVal = fg_wg->eval(coords_tmp);
        COMBIGRID_ERROR_TEST_EQUAL( gridVal , funcVal , 1e-10 , "" );
        coords_tmp = coords;

        coords[0] = min[0] + (max[0] - min[0]) * 1.0;
        coords[1] = min[1] + (max[1] - min[1]) * 1.0;
        coords_tmp = coords;
        funcVal = testfunction2D( coords[0] , coords[1] );
        gridVal = fg->eval(coords_tmp);
        COMBIGRID_ERROR_TEST_EQUAL( gridVal , funcVal , 1e-10 , "" );
        coords_tmp = coords;
        gridVal = fg_wg->eval(coords_tmp);
        COMBIGRID_ERROR_TEST_EQUAL( gridVal , funcVal , 1e-10 , "" );
        coords_tmp = coords;

        coords[0] = min[0] + (max[0] - min[0]) * 0.999;
        coords[1] = min[1] + (max[1] - min[1]) * 0.95;
        coords_tmp = coords;
        funcVal = testfunction2D( coords[0] , coords[1] );
        gridVal = fg->eval(coords_tmp);
        COMBIGRID_ERROR_TEST_EQUAL( gridVal , funcVal , 1e-10 , "" );
        coords_tmp = coords;
        gridVal = fg_wg->eval(coords_tmp);
        COMBIGRID_ERROR_TEST_EQUAL( gridVal , funcVal , 1e-10 , "" );
        coords_tmp = coords;

        coords[0] = min[0] + (max[0] - min[0]) * 0.999;
        coords[1] = min[1] + (max[1] - min[1]) * 0.001;
        coords_tmp = coords;
        funcVal = testfunction2D( coords[0] , coords[1] );
        gridVal = fg->eval(coords_tmp);
        COMBIGRID_ERROR_TEST_EQUAL( gridVal , funcVal , 1e-10 , "" );
        coords_tmp = coords;
        gridVal = fg_wg->eval(coords_tmp);
        COMBIGRID_ERROR_TEST_EQUAL( gridVal , funcVal , 1e-10 , "" );
        coords_tmp = coords;

        // delete the full grid
        delete stretchingMaker;
        delete domain;
        delete fg;
        delete fg_wg;
      }


      static void test12() {
        // test the full grid

        std::vector<double> min(2);
        min[0] = -3;
        min[1] = 2;
        std::vector<double> max(2);
        max[0] = 1;
        max[1] = 4;
        std::vector<int> levels(2);
        levels[0] = 6;
        levels[1] = 7;
        double funcVal = 0.0 , gridVal = 0.0;

        AbstractStretchingMaker* stretchingMaker = new AtanSpecialStretching();
        GridDomain* domain = new GridDomain( 2 , levels , min , max , (*stretchingMaker) );

        FullGrid<double>* fg = new FullGrid<double>(2, 5);
        FullGrid<double>* fg_wg = new FullGrid<double>(2, 5, false);

        fg->createFullGrid();
        fg_wg->createFullGrid();
        fg->setDomain( domain );
        fg_wg->setDomain( domain );

        std::vector<double> coords(2, 0.0);
        std::vector<double> coords_tmp(2, 0.0);

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

        coords[0] = min[0] + (max[0] - min[0]) * 1.0 / 3.0;
        coords[1] = min[1] + (max[1] - min[1]) * 1.0 / 3.0;
        coords_tmp = coords;
        funcVal = testfunction2D( coords[0] , coords[1] );
        gridVal = fg->eval(coords_tmp);
        COMBIGRID_ERROR_TEST_EQUAL( gridVal , funcVal , 1e-10 , "" );
        coords_tmp = coords;
        gridVal = fg_wg->eval(coords_tmp);
        COMBIGRID_ERROR_TEST_EQUAL( gridVal , funcVal , 1e-10 , "" );
        coords_tmp = coords;

        coords[0] = min[0] + (max[0] - min[0]) * 0.0001;
        coords[1] = min[1] + (max[1] - min[1]) * 1.0 / 3.0;
        coords_tmp = coords;
        funcVal = testfunction2D( coords[0] , coords[1] );
        gridVal = fg->eval(coords_tmp);
        COMBIGRID_ERROR_TEST_EQUAL( gridVal , funcVal , 1e-10 , "" );
        coords_tmp = coords;
        gridVal = fg_wg->eval(coords_tmp);
        COMBIGRID_ERROR_TEST_EQUAL( gridVal , funcVal , 1e-10 , "" );
        coords_tmp = coords;

        coords[0] = min[0] + (max[0] - min[0]) * 1.0 / 3.0;
        coords[1] = min[1] + (max[1] - min[1]) * 0.0001;
        coords_tmp = coords;
        funcVal = testfunction2D( coords[0] , coords[1] );
        gridVal = fg->eval(coords_tmp);
        COMBIGRID_ERROR_TEST_EQUAL( gridVal , funcVal , 1e-10 , "" );
        coords_tmp = coords;
        gridVal = fg_wg->eval(coords_tmp);
        COMBIGRID_ERROR_TEST_EQUAL( gridVal , funcVal , 1e-10 , "" );
        coords_tmp = coords;

        coords[0] = min[0] + (max[0] - min[0]) * 0.00001;
        coords[1] = min[1] + (max[1] - min[1]) * 0.00001;
        coords_tmp = coords;
        funcVal = testfunction2D( coords[0] , coords[1] );
        gridVal = fg->eval(coords_tmp);
        COMBIGRID_ERROR_TEST_EQUAL( gridVal , funcVal , 1e-10 , "" );
        coords_tmp = coords;
        gridVal = fg_wg->eval(coords_tmp);
        COMBIGRID_ERROR_TEST_EQUAL( gridVal , funcVal , 1e-10 , "" );
        coords_tmp = coords;

        coords[0] = min[0] + (max[0] - min[0]) * 0.999;
        coords[1] = min[1] + (max[1] - min[1]) * 0.95;
        coords_tmp = coords;
        funcVal = testfunction2D( coords[0] , coords[1] );
        gridVal = fg->eval(coords_tmp);
        COMBIGRID_ERROR_TEST_EQUAL( gridVal , funcVal , 1e-10 , "" );
        coords_tmp = coords;
        gridVal = fg_wg->eval(coords_tmp);
        COMBIGRID_ERROR_TEST_EQUAL( gridVal , funcVal , 1e-10 , "" );
        coords_tmp = coords;

        coords[0] = min[0] + (max[0] - min[0]) * 0.999;
        coords[1] = min[1] + (max[1] - min[1]) * 0.001;
        coords_tmp = coords;
        funcVal = testfunction2D( coords[0] , coords[1] );
        gridVal = fg->eval(coords_tmp);
        COMBIGRID_ERROR_TEST_EQUAL( gridVal , funcVal , 1e-10 , "" );
        coords_tmp = coords;
        gridVal = fg_wg->eval(coords_tmp);
        COMBIGRID_ERROR_TEST_EQUAL( gridVal , funcVal , 1e-10 , "" );
        coords_tmp = coords;

        coords[0] = min[0] + (max[0] - min[0]) * 1.0;
        coords[1] = min[1] + (max[1] - min[1]) * 1.0;
        coords_tmp = coords;
        funcVal = testfunction2D( coords[0] , coords[1] );
        gridVal = fg->eval(coords_tmp);
        COMBIGRID_ERROR_TEST_EQUAL( gridVal , funcVal , 1e-10 , "" );
        coords_tmp = coords;
        gridVal = fg_wg->eval(coords_tmp);
        COMBIGRID_ERROR_TEST_EQUAL( gridVal , funcVal , 1e-10 , "" );
        coords_tmp = coords;

        // delete the full grid
        delete stretchingMaker;
        delete domain;
        delete fg;
        delete fg_wg;
      }


      static void test13() {
        // test the full grid

        std::vector<double> min(2);
        min[0] = -2;
        min[1] = 1;
        std::vector<double> max(2);
        max[0] = 1;
        max[1] = 3;
        std::vector<int> levels(2);
        levels[0] = 8;
        levels[1] = 7;
        double funcVal = 0.0 , gridVal = 0.0;

        AbstractStretchingMaker* stretchingMaker = new TanStretching();
        GridDomain* domain = new GridDomain( 2 , levels , min , max , (*stretchingMaker) );

        FullGrid<double>* fg = new FullGrid<double>(2, 5);
        FullGrid<double>* fg_wg = new FullGrid<double>(2, levels, false);

        fg->createFullGrid();
        fg_wg->createFullGrid();
        fg->setDomain( domain );
        fg_wg->setDomain( domain );

        std::vector<double> coords(2, 0.0);
        std::vector<double> coords_tmp(2, 0.0);

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

        coords[0] = min[0] + (max[0] - min[0]) * 1.0 / 3.0;
        coords[1] = min[1] + (max[1] - min[1]) * 1.0 / 3.0;
        coords_tmp = coords;
        funcVal = testfunction2D( coords[0] , coords[1] );
        gridVal = fg->eval(coords_tmp);
        COMBIGRID_ERROR_TEST_EQUAL( gridVal , funcVal , 1e-10 , "" );
        coords_tmp = coords;
        gridVal = fg_wg->eval(coords_tmp);
        COMBIGRID_ERROR_TEST_EQUAL( gridVal , funcVal , 1e-10 , "" );
        coords_tmp = coords;

        coords[0] = min[0] + (max[0] - min[0]) * 0.0001;
        coords[1] = min[1] + (max[1] - min[1]) * 1.0 / 3.0;
        coords_tmp = coords;
        funcVal = testfunction2D( coords[0] , coords[1] );
        gridVal = fg->eval(coords_tmp);
        COMBIGRID_ERROR_TEST_EQUAL( gridVal , funcVal , 1e-10 , "" );
        coords_tmp = coords;
        gridVal = fg_wg->eval(coords_tmp);
        COMBIGRID_ERROR_TEST_EQUAL( gridVal , funcVal , 1e-10 , "" );
        coords_tmp = coords;

        coords[0] = min[0] + (max[0] - min[0]) * 1.0 / 3.0;
        coords[1] = min[1] + (max[1] - min[1]) * 0.0001;
        coords_tmp = coords;
        funcVal = testfunction2D( coords[0] , coords[1] );
        gridVal = fg->eval(coords_tmp);
        COMBIGRID_ERROR_TEST_EQUAL( gridVal , funcVal , 1e-10 , "" );
        coords_tmp = coords;
        gridVal = fg_wg->eval(coords_tmp);
        COMBIGRID_ERROR_TEST_EQUAL( gridVal , funcVal , 1e-10 , "" );
        coords_tmp = coords;

        coords[0] = min[0] + (max[0] - min[0]) * 0.00001;
        coords[1] = min[1] + (max[1] - min[1]) * 0.00001;
        coords_tmp = coords;
        funcVal = testfunction2D( coords[0] , coords[1] );
        gridVal = fg->eval(coords_tmp);
        COMBIGRID_ERROR_TEST_EQUAL( gridVal , funcVal , 1e-10 , "" );
        coords_tmp = coords;
        gridVal = fg_wg->eval(coords_tmp);
        COMBIGRID_ERROR_TEST_EQUAL( gridVal , funcVal , 1e-10 , "" );
        coords_tmp = coords;

        coords[0] = min[0] + (max[0] - min[0]) * 1.0;
        coords[1] = min[1] + (max[1] - min[1]) * 1.0;
        coords_tmp = coords;
        funcVal = testfunction2D( coords[0] , coords[1] );
        gridVal = fg->eval(coords_tmp);
        COMBIGRID_ERROR_TEST_EQUAL( gridVal , funcVal , 1e-10 , "" );
        coords_tmp = coords;
        gridVal = fg_wg->eval(coords_tmp);
        COMBIGRID_ERROR_TEST_EQUAL( gridVal , funcVal , 1e-10 , "" );
        coords_tmp = coords;

        coords[0] = min[0] + (max[0] - min[0]) * 0.999;
        coords[1] = min[1] + (max[1] - min[1]) * 0.95;
        coords_tmp = coords;
        funcVal = testfunction2D( coords[0] , coords[1] );
        gridVal = fg->eval(coords_tmp);
        COMBIGRID_ERROR_TEST_EQUAL( gridVal , funcVal , 1e-10 , "" );
        coords_tmp = coords;
        gridVal = fg_wg->eval(coords_tmp);
        COMBIGRID_ERROR_TEST_EQUAL( gridVal , funcVal , 1e-10 , "" );
        coords_tmp = coords;

        coords[0] = min[0] + (max[0] - min[0]) * 0.999;
        coords[1] = min[1] + (max[1] - min[1]) * 0.001;
        coords_tmp = coords;
        funcVal = testfunction2D( coords[0] , coords[1] );
        gridVal = fg->eval(coords_tmp);
        COMBIGRID_ERROR_TEST_EQUAL( gridVal , funcVal , 1e-10 , "" );
        coords_tmp = coords;
        gridVal = fg_wg->eval(coords_tmp);
        COMBIGRID_ERROR_TEST_EQUAL( gridVal , funcVal , 1e-10 , "" );
        coords_tmp = coords;

        // delete the full grid
        delete stretchingMaker;
        delete domain;
        delete fg;
        delete fg_wg;
      }


      static void test14() {
        // test the full grid

        std::vector<double> min(2);
        min[0] = -2;
        min[1] = 1;
        std::vector<double> max(2);
        max[0] = 1;
        max[1] = 3;
        std::vector<int> levels(2);
        levels[0] = 8;
        levels[1] = 7;
        double funcVal = 0.0 , gridVal = 0.0;

        GridDomain* domain = new GridDomain( 2 , min , max );

        FullGrid<double>* fg = new FullGrid<double>(2, 5);
        FullGrid<double>* fg_wg = new FullGrid<double>(2, levels, false);

        fg->createFullGrid();
        fg_wg->createFullGrid();
        fg->setDomain( domain );
        fg_wg->setDomain( domain );

        std::vector<double> coords(2, 0.0);
        std::vector<double> coords_tmp(2, 0.0);

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

        coords[0] = min[0] + (max[0] - min[0]) * 1.0 / 3.0;
        coords[1] = min[1] + (max[1] - min[1]) * 1.0 / 3.0;
        coords_tmp = coords;
        funcVal = testfunction2D( coords[0] , coords[1] );
        gridVal = fg->eval(coords_tmp);
        COMBIGRID_ERROR_TEST_EQUAL( gridVal , funcVal , 1e-10 , "" );
        coords_tmp = coords;
        gridVal = fg_wg->eval(coords_tmp);
        COMBIGRID_ERROR_TEST_EQUAL( gridVal , funcVal , 1e-10 , "" );
        coords_tmp = coords;

        coords[0] = min[0] + (max[0] - min[0]) * 0.0001;
        coords[1] = min[1] + (max[1] - min[1]) * 1.0 / 3.0;
        coords_tmp = coords;
        funcVal = testfunction2D( coords[0] , coords[1] );
        gridVal = fg->eval(coords_tmp);
        COMBIGRID_ERROR_TEST_EQUAL( gridVal , funcVal , 1e-10 , "" );
        coords_tmp = coords;
        gridVal = fg_wg->eval(coords_tmp);
        COMBIGRID_ERROR_TEST_EQUAL( gridVal , funcVal , 1e-10 , "" );
        coords_tmp = coords;

        coords[0] = min[0] + (max[0] - min[0]) * 1.0 / 3.0;
        coords[1] = min[1] + (max[1] - min[1]) * 0.0001;
        coords_tmp = coords;
        funcVal = testfunction2D( coords[0] , coords[1] );
        gridVal = fg->eval(coords_tmp);
        COMBIGRID_ERROR_TEST_EQUAL( gridVal , funcVal , 1e-10 , "" );
        coords_tmp = coords;
        gridVal = fg_wg->eval(coords_tmp);
        COMBIGRID_ERROR_TEST_EQUAL( gridVal , funcVal , 1e-10 , "" );
        coords_tmp = coords;

        coords[0] = min[0] + (max[0] - min[0]) * 0.00001;
        coords[1] = min[1] + (max[1] - min[1]) * 0.00001;
        coords_tmp = coords;
        funcVal = testfunction2D( coords[0] , coords[1] );
        gridVal = fg->eval(coords_tmp);
        COMBIGRID_ERROR_TEST_EQUAL( gridVal , funcVal , 1e-10 , "" );
        coords_tmp = coords;
        gridVal = fg_wg->eval(coords_tmp);
        COMBIGRID_ERROR_TEST_EQUAL( gridVal , funcVal , 1e-10 , "" );
        coords_tmp = coords;

        coords[0] = min[0] + (max[0] - min[0]) * 1.0;
        coords[1] = min[1] + (max[1] - min[1]) * 1.0;
        coords_tmp = coords;
        funcVal = testfunction2D( coords[0] , coords[1] );
        gridVal = fg->eval(coords_tmp);
        COMBIGRID_ERROR_TEST_EQUAL( gridVal , funcVal , 1e-10 , "" );
        coords_tmp = coords;
        gridVal = fg_wg->eval(coords_tmp);
        COMBIGRID_ERROR_TEST_EQUAL( gridVal , funcVal , 1e-10 , "" );
        coords_tmp = coords;

        coords[0] = min[0] + (max[0] - min[0]) * 0.999;
        coords[1] = min[1] + (max[1] - min[1]) * 0.95;
        coords_tmp = coords;
        funcVal = testfunction2D( coords[0] , coords[1] );
        gridVal = fg->eval(coords_tmp);
        COMBIGRID_ERROR_TEST_EQUAL( gridVal , funcVal , 1e-10 , "" );
        coords_tmp = coords;
        gridVal = fg_wg->eval(coords_tmp);
        COMBIGRID_ERROR_TEST_EQUAL( gridVal , funcVal , 1e-10 , "" );
        coords_tmp = coords;

        coords[0] = min[0] + (max[0] - min[0]) * 0.999;
        coords[1] = min[1] + (max[1] - min[1]) * 0.001;
        coords_tmp = coords;
        funcVal = testfunction2D( coords[0] , coords[1] );
        gridVal = fg->eval(coords_tmp);
        COMBIGRID_ERROR_TEST_EQUAL( gridVal , funcVal , 1e-10 , "" );
        coords_tmp = coords;
        gridVal = fg_wg->eval(coords_tmp);
        COMBIGRID_ERROR_TEST_EQUAL( gridVal , funcVal , 1e-10 , "" );
        coords_tmp = coords;

        // delete the full grid
        delete domain;
        delete fg;
        delete fg_wg;
      }

      static void test2() {
        // test the full grid

        std::vector<double> min(3);
        min[0] = -1;
        min[1] = 3;
        min[2] = -2;
        std::vector<double> max(3);
        max[0] = 2;
        max[1] = 5;
        max[2] = 2;
        std::vector<int> levels(3, 9);
        levels[0] = 7;
        levels[1] = 8;
        levels[2] = 9;
        double funcVal = 0.0 , gridVal = 0.0;

        AbstractStretchingMaker* stretchingMaker = new AtanSpecialStretching(); //new UniformStretching(); //new TanStretching(0.25);
        GridDomain* domain = new GridDomain( 3 , levels , min , max , (*stretchingMaker) );

        std::vector<bool> hasBdrPoints(3, true);
        levels[1] = 7;
        levels[1] = 4;
        levels[2] = 7;
        hasBdrPoints[1] = false;
        FullGrid<double>* fg = new FullGrid<double>( 3 , levels , hasBdrPoints );
        levels[1] = 3;
        levels[1] = 6;
        levels[2] = 8;
        hasBdrPoints[0] = false;
        hasBdrPoints[1] = false;
        hasBdrPoints[2] = false;
        FullGrid<double>* fg_wg = new FullGrid<double>( 3 , levels , hasBdrPoints );

        fg->createFullGrid();
        fg_wg->createFullGrid();
        fg->setDomain( domain );
        fg_wg->setDomain( domain );

        std::vector<double> coords(3, 0.0);
        std::vector<double> coords_tmp(3, 0.0);

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

        coords[0] = min[0] + (max[0] - min[0]) * 1.0 / 3.0;
        coords[1] = min[1] + (max[1] - min[1]) * 1.0 / 3.0;
        coords[2] = min[2] + (max[2] - min[2]) * 1.0 / 3.0;
        coords_tmp = coords;
        funcVal =  testfunction3D( coords[0] , coords[1] , coords[2] ) ;
        gridVal = fg->eval(coords_tmp);
        COMBIGRID_ERROR_TEST_EQUAL( gridVal , funcVal , 1e-10 , "" );
        coords_tmp = coords;
        gridVal = fg_wg->eval(coords_tmp);
        COMBIGRID_ERROR_TEST_EQUAL( gridVal , funcVal , 1e-10 , "" );
        coords_tmp = coords;

        coords[0] = min[0] + (max[0] - min[0]) * 0.0001;
        coords[1] = min[1] + (max[1] - min[1]) * 1.0 / 3.0;
        coords[2] = min[2] + (max[2] - min[2]) * 1.0 / 3.0;
        coords_tmp = coords;
        funcVal =  testfunction3D( coords[0] , coords[1] , coords[2] ) ;
        gridVal = fg->eval(coords_tmp);
        COMBIGRID_ERROR_TEST_EQUAL( gridVal , funcVal , 1e-10 , "" );
        coords_tmp = coords;
        gridVal = fg_wg->eval(coords_tmp);
        COMBIGRID_ERROR_TEST_EQUAL( gridVal , funcVal , 1e-10 , "" );
        coords_tmp = coords;

        coords[0] = min[0] + (max[0] - min[0]) * 1.0 / 3.0;
        coords[1] = min[1] + (max[1] - min[1]) * 0.0001;
        coords[2] = min[2] + (max[2] - min[2]) * 0.0001;
        coords_tmp = coords;
        funcVal =  testfunction3D( coords[0] , coords[1] , coords[2] ) ;
        gridVal = fg->eval(coords_tmp);
        COMBIGRID_ERROR_TEST_EQUAL( gridVal , funcVal , 1e-10 , "" );
        coords_tmp = coords;
        gridVal = fg_wg->eval(coords_tmp);
        COMBIGRID_ERROR_TEST_EQUAL( gridVal , funcVal , 1e-10 , "" );
        coords_tmp = coords;

        coords[0] = min[0] + (max[0] - min[0]) * 0.0;
        coords[1] = min[1] + (max[1] - min[1]) * 0.0;
        coords[2] = min[2] + (max[2] - min[2]) * 0.0;
        coords_tmp = coords;
        funcVal =  testfunction3D( coords[0] , coords[1] , coords[2] ) ;
        gridVal = fg->eval(coords_tmp);
        COMBIGRID_ERROR_TEST_EQUAL( gridVal , funcVal , 1e-10 , "" );
        coords_tmp = coords;
        gridVal = fg_wg->eval(coords_tmp);
        COMBIGRID_ERROR_TEST_EQUAL( gridVal , funcVal , 1e-10 , "" );
        coords_tmp = coords;

        coords[0] = min[0] + (max[0] - min[0]) * 1.0;
        coords[1] = min[1] + (max[1] - min[1]) * 1.0;
        coords[2] = min[2] + (max[2] - min[2]) * 1.0;
        coords_tmp = coords;
        funcVal =  testfunction3D( coords[0] , coords[1] , coords[2] ) ;
        gridVal = fg->eval(coords_tmp);
        COMBIGRID_ERROR_TEST_EQUAL( gridVal , funcVal , 1e-10 , "" );
        coords_tmp = coords;
        gridVal = fg_wg->eval(coords_tmp);
        COMBIGRID_ERROR_TEST_EQUAL( gridVal , funcVal , 1e-10 , "" );
        coords_tmp = coords;

        coords[0] = min[0] + (max[0] - min[0]) * 0.999;
        coords[1] = min[1] + (max[1] - min[1]) * 0.95;
        coords[2] = min[2] + (max[2] - min[2]) * 0.95;
        coords_tmp = coords;
        funcVal = testfunction3D( coords[0] , coords[1] , coords[2] );
        gridVal = fg->eval(coords_tmp);
        COMBIGRID_ERROR_TEST_EQUAL( gridVal , funcVal , 1e-10 , "" );
        coords_tmp = coords;
        gridVal = fg_wg->eval(coords_tmp);
        COMBIGRID_ERROR_TEST_EQUAL( gridVal , funcVal , 1e-10 , "" );
        coords_tmp = coords;

        coords[0] = min[0] + (max[0] - min[0]) * 0.999;
        coords[1] = min[1] + (max[1] - min[1]) * 0.001;
        coords[2] = min[2] + (max[2] - min[2]) * 0.999;
        coords_tmp = coords;
        funcVal =  testfunction3D( coords[0] , coords[1] , coords[2] );
        gridVal = fg->eval(coords_tmp);
        COMBIGRID_ERROR_TEST_EQUAL( gridVal , funcVal , 1e-10 , "" );
        coords_tmp = coords;
        gridVal = fg_wg->eval(coords_tmp);
        COMBIGRID_ERROR_TEST_EQUAL( gridVal , funcVal , 1e-10 , "" );
        coords_tmp = coords;

        // delete the full grid
        delete stretchingMaker;
        delete domain;
        delete fg;
        delete fg_wg;
      }
  };
}

#endif /* TESTSTRETCHING_HPP_ */