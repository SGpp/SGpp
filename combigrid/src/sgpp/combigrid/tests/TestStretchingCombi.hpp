// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#ifndef TESTTRETCHINGCOMBI_HPP_
#define TESTTRETCHINGCOMBI_HPP_

#include <combigrid.hpp>

namespace combigrid {

/** Class to thes the all the full grids functionalities */
class TestStretchingCombi {
 public:

  /** empty Ctor */
  TestStretchingCombi();

  static void test_all_cases() {
    COMBIGRID_OUT_LEVEL3( 4 , "TestStretching Combi ... ");
    // 2D
    test1();
    test12();

    // 3D test
    test2();
    test21();
  }

  /** test function for the grid, this should be exactly interpolated */
  static double testfunction2D( double x1, double x2) {
    //COMBIGRID_OUT_LEVEL3( 4 , "testfunction2D x1=" << x1 << " , x2=" << x2 << " ret: " << (2+3.0*x1+0.5*x2));
    return (2 + 3.0 * x1 + 0.5 * x2);
  }

  static double testfunction2D1( double x1, double x2) {
    return ::exp(0.25 * x1 + 0.5 * x2);
  }

  static double testfunction3D( double x1, double x2, double x3) {
    //COMBIGRID_OUT_LEVEL3( 4 , "testfunction3D x1=" << x1 << " , x2=" << x2 << " , x3=" << x3 << " ret: " << (2+3.0*x1+0.5*x2 + 1.5*x3));
    return (2 + 3.0 * x1 + 0.5 * x2 + 1.5 * x3);
  }

  static double testfunction3D1( double x1, double x2, double x3) {
    return (0.01 + 1.0 * x1 + 0.5 * x2 + 1.5 * x3) * (0.01 + 2.0 * x1 + 0.5 * x2 +
           1.5 * x3);
  }

  static void test1() {
    std::vector<double> coords( 2 , 0.0);
    std::vector<double> coords_tmp(2, 0.0);
    std::vector<int> minlevels(2);
    minlevels[0] = 4;
    minlevels[1] = 4;
    std::vector<int> maxlevels(2);
    maxlevels[0] = 6;
    maxlevels[1] = 6;
    CombiSchemeBasis* combischeme = new TS_CT( minlevels, maxlevels);
    AbstractCombiGrid* combigrid = new SerialCombiGrid(combischeme);
    AbstractCombiGrid* combigrid1 = new SerialCombiGrid(combischeme, false );

    std::vector<double> min(2);
    min[0] = -2;
    min[1] = 2;
    std::vector<double> max(2);
    max[0] = 2;
    max[1] = 5;
    std::vector<int> levels(2, 6);
    double funcVal = 0.0 , gridVal = 0.0;

    AbstractStretchingMaker* stretchingMaker = new
    AtanSpecialStretching(); //new UniformStretching(); //new AtanSpecialStretching();
    GridDomain* domain = new GridDomain( 2 , levels , min , max ,
                                         (*stretchingMaker) );

    // creates the full grids in the combination grid
    combigrid->createFullGrids();
    combigrid1->createFullGrids();

    // for each full grid set the values
    for ( int fg = 0 ; fg < combigrid->getNrFullGrid() ; fg++) {
      FullGridD* fgp = combigrid->getFullGrid(fg);
      fgp->setDomain(domain);

      // for this full grid set the values
      for (int i = 0 ; i < fgp->getNrElements() ; i++) {
        fgp->getCoords( i , coords );
        fgp->getElementVector()[i] = testfunction2D( coords[0] , coords[1] );
      }

      fgp->setDomain(NULL);
    }

    for ( int fg = 0 ; fg < combigrid1->getNrFullGrid() ; fg++) {
      FullGridD* fgp = combigrid1->getFullGrid(fg);
      // for this full grid set the values
      fgp->setDomain(domain);

      for (int i = 0 ; i < fgp->getNrElements() ; i++) {
        fgp->getCoords( i , coords );
        fgp->getElementVector()[i] = testfunction2D( coords[0] , coords[1] );
      }

      fgp->setDomain(NULL);
    }

    combigrid->setDomain(domain);
    combigrid1->setDomain(domain);

    coords[0] = min[0] + (max[0] - min[0]) * 1.0 / 3.0;
    coords[1] = min[1] + (max[1] - min[1]) * 1.0 / 3.0;
    coords_tmp = coords;
    funcVal = testfunction2D( coords[0] , coords[1] );
    gridVal = combigrid->eval(coords_tmp);
    COMBIGRID_ERROR_TEST_EQUAL( gridVal , funcVal , 1e-10 , "" );
    coords_tmp = coords;
    gridVal = combigrid1->eval(coords_tmp);
    COMBIGRID_ERROR_TEST_EQUAL( gridVal , funcVal , 1e-10 , "" );
    coords_tmp = coords;

    coords[0] = min[0] + (max[0] - min[0]) * 0.4001;
    coords[1] = min[1] + (max[1] - min[1]) * 1.0 / 3.0;
    coords_tmp = coords;
    funcVal = testfunction2D( coords[0] , coords[1] );
    gridVal = combigrid->eval(coords_tmp);
    COMBIGRID_ERROR_TEST_EQUAL( gridVal , funcVal , 1e-10 , "" );
    coords_tmp = coords;
    gridVal = combigrid1->eval(coords_tmp);
    COMBIGRID_ERROR_TEST_EQUAL( gridVal , funcVal , 1e-10 , "" );
    coords_tmp = coords;

    coords[0] = min[0] + (max[0] - min[0]) * 1.0 / 3.0;
    coords[1] = min[1] + (max[1] - min[1]) * 0.6001;
    coords_tmp = coords;
    funcVal = testfunction2D( coords[0] , coords[1] );
    gridVal = combigrid->eval(coords_tmp);
    COMBIGRID_ERROR_TEST_EQUAL( gridVal , funcVal , 1e-10 , "" );
    coords_tmp = coords;
    gridVal = combigrid1->eval(coords_tmp);
    COMBIGRID_ERROR_TEST_EQUAL( gridVal , funcVal , 1e-10 , "" );
    coords_tmp = coords;

    coords[0] = min[0] + (max[0] - min[0]) * 0.61001;
    coords[1] = min[1] + (max[1] - min[1]) * 0.72001;
    coords_tmp = coords;
    funcVal = testfunction2D( coords[0] , coords[1] );
    gridVal = combigrid->eval(coords_tmp);
    COMBIGRID_ERROR_TEST_EQUAL( gridVal , funcVal , 1e-10 , "" );
    coords_tmp = coords;
    gridVal = combigrid1->eval(coords_tmp);
    COMBIGRID_ERROR_TEST_EQUAL( gridVal , funcVal , 1e-10 , "" );
    coords_tmp = coords;

    coords[0] = min[0] + (max[0] - min[0]) * 0.730;
    coords[1] = min[1] + (max[1] - min[1]) * 0.67;
    coords_tmp = coords;
    funcVal = testfunction2D( coords[0] , coords[1] );
    gridVal = combigrid->eval(coords_tmp);
    COMBIGRID_ERROR_TEST_EQUAL( gridVal , funcVal , 1e-10 , "" );
    coords_tmp = coords;
    gridVal = combigrid1->eval(coords_tmp);
    COMBIGRID_ERROR_TEST_EQUAL( gridVal , funcVal , 1e-10 , "" );
    coords_tmp = coords;

    coords[0] = min[0] + (max[0] - min[0]) * 0.4999;
    coords[1] = min[1] + (max[1] - min[1]) * 0.495;
    coords_tmp = coords;
    funcVal = testfunction2D( coords[0] , coords[1] );
    gridVal = combigrid->eval(coords_tmp);
    COMBIGRID_ERROR_TEST_EQUAL( gridVal , funcVal , 1e-10 , "" );
    coords_tmp = coords;
    gridVal = combigrid1->eval(coords_tmp);
    COMBIGRID_ERROR_TEST_EQUAL( gridVal , funcVal , 1e-10 , "" );
    coords_tmp = coords;

    coords[0] = min[0] + (max[0] - min[0]) * 0.5999;
    coords[1] = min[1] + (max[1] - min[1]) * 0.7201;
    coords_tmp = coords;
    funcVal = testfunction2D( coords[0] , coords[1] );
    gridVal = combigrid->eval(coords_tmp);
    COMBIGRID_ERROR_TEST_EQUAL( gridVal , funcVal , 1e-10 , "" );
    coords_tmp = coords;
    gridVal = combigrid1->eval(coords_tmp);
    COMBIGRID_ERROR_TEST_EQUAL( gridVal , funcVal , 1e-10 , "" );
    coords_tmp = coords;

    delete stretchingMaker;
    delete domain;
    delete combigrid;
    delete combigrid1;
    delete combischeme;
  }


  static void test12() {
    std::vector<double> coords( 2 , 0.0);
    std::vector<double> coords_tmp(2, 0.0);
    CombiSchemeBasis* combischeme = new S_CT( 2 , 10 , 3 );
    AbstractCombiGrid* combigrid = new SerialCombiGrid(combischeme , true );
    AbstractCombiGrid* combigrid1 = new SerialCombiGrid(combischeme , false );

    std::vector<double> min(2);
    min[0] = -4;
    min[1] = 3;
    std::vector<double> max(2);
    max[0] = -2;
    max[1] = 5;
    std::vector<int> levels(2, 10);
    double funcVal = 0.0 , gridVal = 0.0;

    AbstractStretchingMaker* stretchingMaker = new AtanSpecialStretching();
    GridDomain* domain = new GridDomain( 2 , levels , min , max ,
                                         (*stretchingMaker) );

    // creates the full grids in the combination grid
    combigrid->createFullGrids();
    combigrid1->createFullGrids();

    // set for all full grids the domain
    combigrid->setDomainAllFG(domain);
    combigrid1->setDomainAllFG(domain);

    // for each full grid set the values
    for ( int fg = 0 ; fg < combigrid->getNrFullGrid() ; fg++) {
      FullGridD* fgp = combigrid->getFullGrid(fg);

      for (int i = 0 ; i < fgp->getNrElements() ; i++) {
        fgp->getCoords( i , coords );
        fgp->getElementVector()[i] = testfunction2D1( coords[0] , coords[1] );
      }
    }

    // for each full grid set the values
    for ( int fg = 0 ; fg < combigrid1->getNrFullGrid() ; fg++) {
      FullGridD* fgp = combigrid1->getFullGrid(fg);

      for (int i = 0 ; i < fgp->getNrElements() ; i++) {
        fgp->getCoords( i , coords );
        fgp->getElementVector()[i] = testfunction2D1( coords[0] , coords[1] );
      }
    }

    coords[0] = min[0] + (max[0] - min[0]) * 1.0 / 3.0;
    coords[1] = min[1] + (max[1] - min[1]) * 1.0 / 3.0;
    coords_tmp = coords;
    funcVal = testfunction2D1( coords[0] , coords[1] );
    gridVal = combigrid->eval(coords_tmp);
    COMBIGRID_ERROR_TEST_EQUAL(gridVal , funcVal  , 1e-4 , "" );
    coords_tmp = coords;
    gridVal = combigrid1->eval(coords_tmp);
    COMBIGRID_ERROR_TEST_EQUAL(gridVal , funcVal  , 1e-4 , "" );
    coords_tmp = coords;

    coords[0] = min[0] + (max[0] - min[0]) * 0.1001;
    coords[1] = min[1] + (max[1] - min[1]) * 1.0 / 3.0;
    coords_tmp = coords;
    funcVal = testfunction2D1( coords[0] , coords[1] );
    gridVal = combigrid->eval(coords_tmp);
    COMBIGRID_ERROR_TEST_EQUAL( gridVal , funcVal  , 1e-4 , "" );
    coords_tmp = coords;
    gridVal = combigrid1->eval(coords_tmp);
    COMBIGRID_ERROR_TEST_EQUAL(gridVal , funcVal  , 1e-4 , "" );
    coords_tmp = coords;

    coords[0] = min[0] + (max[0] - min[0]) * 0.9001;
    coords[1] = min[1] + (max[1] - min[1]) * 2.0 / 3.0;
    coords_tmp = coords;
    funcVal = testfunction2D1( coords[0] , coords[1] );
    gridVal = combigrid->eval(coords_tmp);
    COMBIGRID_ERROR_TEST_EQUAL( gridVal , funcVal  , 1e-4 , "" );
    coords_tmp = coords;
    gridVal = combigrid1->eval(coords_tmp);
    COMBIGRID_ERROR_TEST_EQUAL(gridVal , funcVal  , 1e-4 , "" );
    coords_tmp = coords;

    coords[0] = min[0] + (max[0] - min[0]) * 0.0001;
    coords[1] = min[1] + (max[1] - min[1]) * 1.0 / 3.0;
    coords_tmp = coords;
    funcVal = testfunction2D1( coords[0] , coords[1] );
    gridVal = combigrid->eval(coords_tmp);
    COMBIGRID_ERROR_TEST_EQUAL( gridVal , funcVal  , 1e-4 , "" );
    coords_tmp = coords;

    coords[0] = min[0] + (max[0] - min[0]) * 1.0 / 3.0;
    coords[1] = min[1] + (max[1] - min[1]) * 0.0001;
    coords_tmp = coords;
    funcVal = testfunction2D1( coords[0] , coords[1] );
    gridVal = combigrid->eval(coords_tmp);
    COMBIGRID_ERROR_TEST_EQUAL( gridVal , funcVal  , 1e-4 , "" );
    coords_tmp = coords;

    coords[0] = min[0] + (max[0] - min[0]) * 0.00001;
    coords[1] = min[1] + (max[1] - min[1]) * 0.00001;
    coords_tmp = coords;
    funcVal = testfunction2D1( coords[0] , coords[1] );
    gridVal = combigrid->eval(coords_tmp);
    COMBIGRID_ERROR_TEST_EQUAL( gridVal , funcVal  , 1e-4 , "" );
    coords_tmp = coords;

    coords[0] = min[0] + (max[0] - min[0]) * 1.0;
    coords[1] = min[1] + (max[1] - min[1]) * 1.0;
    coords_tmp = coords;
    funcVal = testfunction2D1( coords[0] , coords[1] );
    gridVal = combigrid->eval(coords_tmp);
    COMBIGRID_ERROR_TEST_EQUAL( gridVal , funcVal  , 1e-4 , "" );
    coords_tmp = coords;

    coords[0] = min[0] + (max[0] - min[0]) * 0.999;
    coords[1] = min[1] + (max[1] - min[1]) * 0.95;
    coords_tmp = coords;
    funcVal = testfunction2D1( coords[0] , coords[1] );
    gridVal = combigrid->eval(coords_tmp);
    COMBIGRID_ERROR_TEST_EQUAL( gridVal , funcVal  , 1e-4 , "" );
    coords_tmp = coords;

    coords[0] = min[0] + (max[0] - min[0]) * 0.999;
    coords[1] = min[1] + (max[1] - min[1]) * 0.001;
    coords_tmp = coords;
    funcVal = testfunction2D1( coords[0] , coords[1] );
    gridVal = combigrid->eval(coords_tmp);
    COMBIGRID_ERROR_TEST_EQUAL( gridVal , funcVal  , 1e-4 , "" );
    coords_tmp = coords;

    delete stretchingMaker;
    delete domain;
    delete combigrid;
    delete combischeme;
  }



  static void test2() {

    std::vector<double> coords(3, 0.0);
    std::vector<double> coords_tmp(3, 0.0);
    std::vector<int> levels( 3 , 6);
    CombiSchemeBasis* combischeme = new S_CT( 3 , 6 , 2);
    //std::vector<int> minlevels(3); minlevels[0]=2; minlevels[1]=2; minlevels[2]=2;
    //std::vector<int> maxlevels(3); maxlevels[0]=3; maxlevels[1]=3; maxlevels[2]=3;
    //CombiSchemeBasis* combischeme = new TS_CT( minlevels, maxlevels);
    AbstractCombiGrid* combigrid = new SerialCombiGrid(combischeme, true );
    AbstractCombiGrid* combigrid1 = new SerialCombiGrid(combischeme, false );

    std::vector<double> min(3);
    min[0] = -1;
    min[1] = 0.1;
    min[2] = -0.3;
    std::vector<double> max(3);
    max[0] = 0.1;
    max[1] = 0.3;
    max[2] = 0.8;
    std::vector<int> levels_d(3, 7);
    double funcVal = 0.0 , gridVal = 0.0;
    AbstractStretchingMaker* stretchingMaker = new TanStretching(
      0.3); //new UniformStretching(); //new TanStretching(0.7);
    GridDomain* domain = new GridDomain( 3 , levels_d , min , max ,
                                         (*stretchingMaker) );

    // creates the full grids in the combination grid
    combigrid->createFullGrids();
    combigrid1->createFullGrids();

    // set for all full grids the domain
    combigrid->setDomainAllFG(domain);
    combigrid1->setDomainAllFG(domain);

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

    coords[0] = min[0] + (max[0] - min[0]) * 1.0 / 3.0;
    coords[1] = min[1] + (max[1] - min[1]) * 1.0 / 3.0;
    coords[2] = min[2] + (max[2] - min[2]) * 1.0 / 3.0;
    coords_tmp = coords;
    funcVal = testfunction3D( coords[0] , coords[1] , coords[2]);
    gridVal = combigrid->eval(coords_tmp);
    COMBIGRID_ERROR_TEST_EQUAL(  gridVal , funcVal , 1e-10 , "" );
    coords_tmp = coords;
    gridVal = combigrid1->eval(coords_tmp);
    COMBIGRID_ERROR_TEST_EQUAL(  gridVal , funcVal , 1e-10 , "" );
    coords_tmp = coords;

    coords[0] = min[0] + (max[0] - min[0]) * 0.1001;
    coords[1] = min[1] + (max[1] - min[1]) * 1.0 / 3.0;
    coords[2] = min[2] + (max[2] - min[2]) * 1.1 / 3.0;
    coords_tmp = coords;
    funcVal = testfunction3D( coords[0] , coords[1] , coords[2]);
    gridVal = combigrid->eval(coords_tmp);
    COMBIGRID_ERROR_TEST_EQUAL(  gridVal , funcVal , 1e-10 , "" );
    coords_tmp = coords;
    gridVal = combigrid1->eval(coords_tmp);
    COMBIGRID_ERROR_TEST_EQUAL(  gridVal , funcVal , 1e-10 , "" );
    coords_tmp = coords;

    coords[0] = min[0] + (max[0] - min[0]) * 1.1 / 3.0;
    coords[1] = min[1] + (max[1] - min[1]) * 0.2301;
    coords[2] = min[2] + (max[2] - min[2]) * 0.3001;
    coords_tmp = coords;
    funcVal = testfunction3D( coords[0] , coords[1] , coords[2]);
    gridVal = combigrid->eval(coords_tmp);
    COMBIGRID_ERROR_TEST_EQUAL(  gridVal , funcVal , 1e-10 , "" );
    coords_tmp = coords;
    gridVal = combigrid1->eval(coords_tmp);
    COMBIGRID_ERROR_TEST_EQUAL(  gridVal , funcVal , 1e-10 , "" );
    coords_tmp = coords;

    coords[0] = min[0] + (max[0] - min[0]) * 0.999;
    coords[1] = min[1] + (max[1] - min[1]) * 0.95;
    coords[2] = min[2] + (max[2] - min[2]) * 0.95;
    coords_tmp = coords;
    funcVal = testfunction3D( coords[0] , coords[1] , coords[2]);
    gridVal = combigrid->eval(coords_tmp);
    COMBIGRID_ERROR_TEST_EQUAL(  gridVal , funcVal , 1e-10 , "" );
    coords_tmp = coords;
    gridVal = combigrid1->eval(coords_tmp);
    COMBIGRID_ERROR_TEST_EQUAL(  gridVal , funcVal , 1e-10 , "" );
    coords_tmp = coords;

    coords[0] = min[0] + (max[0] - min[0]) * 0.0;
    coords[1] = min[1] + (max[1] - min[1]) * 0.0;
    coords[2] = min[2] + (max[2] - min[2]) * 0.0;
    coords_tmp = coords;
    funcVal = testfunction3D( coords[0] , coords[1] , coords[2]);
    gridVal = combigrid->eval(coords_tmp);
    COMBIGRID_ERROR_TEST_EQUAL(  gridVal , funcVal , 1e-10 , "" );
    coords_tmp = coords;
    gridVal = combigrid1->eval(coords_tmp);
    COMBIGRID_ERROR_TEST_EQUAL(  gridVal , funcVal , 1e-10 , "" );
    coords_tmp = coords;

    coords[0] = min[0] + (max[0] - min[0]) * 1.0;
    coords[1] = min[1] + (max[1] - min[1]) * 1.0;
    coords[2] = min[2] + (max[2] - min[2]) * 1.0;
    coords_tmp = coords;
    funcVal = testfunction3D( coords[0] , coords[1] , coords[2]);
    gridVal = combigrid->eval(coords_tmp);
    COMBIGRID_ERROR_TEST_EQUAL(  gridVal , funcVal , 1e-10 , "" );
    coords_tmp = coords;
    gridVal = combigrid1->eval(coords_tmp);
    COMBIGRID_ERROR_TEST_EQUAL(  gridVal , funcVal , 1e-10 , "" );
    coords_tmp = coords;

    coords[0] = min[0] + (max[0] - min[0]) * 0.999;
    coords[1] = min[1] + (max[1] - min[1]) * 0.001;
    coords[2] = min[2] + (max[2] - min[2]) * 0.999;
    coords_tmp = coords;
    funcVal = testfunction3D( coords[0] , coords[1] , coords[2]);
    gridVal = combigrid->eval(coords_tmp);
    COMBIGRID_ERROR_TEST_EQUAL(  gridVal , funcVal , 1e-10 , "" );
    coords_tmp = coords;
    gridVal = combigrid1->eval(coords_tmp);
    COMBIGRID_ERROR_TEST_EQUAL(  gridVal , funcVal , 1e-10 , "" );
    coords_tmp = coords;

    delete stretchingMaker;
    delete domain;
    delete combigrid;
    delete combigrid1;
    delete combischeme;
  }



  static void test21() {

    std::vector<double> coords( 3 , 0.0);
    std::vector<double> coords_tmp( 3 , 0.0);
    std::vector<int> levels( 3 , 9);
    std::vector<int> levels_trunc( 3 , 3);
    levels[0] = 7;
    levels[1] = 9;
    levels[2] = 8;
    levels_trunc[0] = 3;
    levels_trunc[1] = 2;
    levels_trunc[2] = 1;
    CombiSchemeBasis* combischeme = new S_CT( 3 , levels , levels_trunc );
    AbstractCombiGrid* combigrid = new SerialCombiGrid(combischeme);
    AbstractCombiGrid* combigrid1 = new SerialCombiGrid(combischeme, false );

    std::vector<double> min(3);
    min[0] = 0.1;
    min[1] = 0.3;
    min[2] = -0.9;
    std::vector<double> max(3);
    max[0] = 1.3;
    max[1] = 1.1;
    max[2] = -0.1;
    std::vector<int> levels_d(3, 9);
    double funcVal = 0.0 , gridVal = 0.0;
    AbstractStretchingMaker* stretchingMaker = new TanStretching(
      0.7); //new UniformStretching();
    GridDomain* domain = new GridDomain( 3 , levels_d , min , max ,
                                         (*stretchingMaker) );

    // creates the full grids in the combination grid
    combigrid->createFullGrids();
    combigrid1->createFullGrids();

    // set for all full grids the domain
    combigrid->setDomainAllFG(domain);
    combigrid1->setDomainAllFG(domain);

    // for each full grid set the values
    for ( int fg = 0 ; fg < combigrid->getNrFullGrid() ; fg++) {
      FullGridD* fgp = combigrid->getFullGrid(fg);

      // for this full grid set the values
      for (int i = 0 ; i < fgp->getNrElements() ; i++) {
        fgp->getCoords( i , coords );
        fgp->getElementVector()[i] = testfunction3D1( coords[0] , coords[1] ,
                                     coords[2]);
      }
    }

    for ( int fg = 0 ; fg < combigrid1->getNrFullGrid() ; fg++) {
      FullGridD* fgp = combigrid1->getFullGrid(fg);

      // for this full grid set the values
      for (int i = 0 ; i < fgp->getNrElements() ; i++) {
        fgp->getCoords( i , coords );
        fgp->getElementVector()[i] = testfunction3D1( coords[0] , coords[1] ,
                                     coords[2]);
      }
    }

    coords[0] = min[0] + (max[0] - min[0]) * 1.0 / 3.0;
    coords[1] = min[1] + (max[1] - min[1]) * 2.0 / 3.0;
    coords[2] = min[2] + (max[2] - min[2]) * 1.6 / 3.0;
    coords_tmp = coords;
    funcVal = testfunction3D1( coords[0] , coords[1] , coords[2]);
    gridVal = combigrid->eval(coords_tmp);
    COMBIGRID_ERROR_TEST_EQUAL( gridVal , funcVal , 1e-4 , "" );
    coords_tmp = coords;
    gridVal = combigrid1->eval(coords_tmp);
    COMBIGRID_ERROR_TEST_EQUAL( gridVal , funcVal , 1e-4 , "" );
    coords_tmp = coords;

    coords[0] = min[0] + (max[0] - min[0]) * 0.1001;
    coords[1] = min[1] + (max[1] - min[1]) * 1.0 / 3.0;
    coords[2] = min[2] + (max[2] - min[2]) * 1.0 / 3.0;
    coords_tmp = coords;
    funcVal = testfunction3D1( coords[0] , coords[1] , coords[2]);
    gridVal = combigrid->eval(coords_tmp);
    COMBIGRID_ERROR_TEST_EQUAL( gridVal , funcVal , 1e-4 , "" );
    coords_tmp = coords;
    gridVal = combigrid1->eval(coords_tmp);
    COMBIGRID_ERROR_TEST_EQUAL( gridVal , funcVal , 1e-4 , "" );
    coords_tmp = coords;

    coords[0] = min[0] + (max[0] - min[0]) * 1.0 / 3.0;
    coords[1] = min[1] + (max[1] - min[1]) * 0.2101;
    coords[2] = min[2] + (max[2] - min[2]) * 0.2301;
    coords_tmp = coords;
    funcVal = testfunction3D1( coords[0] , coords[1] , coords[2]);
    gridVal = combigrid->eval(coords_tmp);
    COMBIGRID_ERROR_TEST_EQUAL( gridVal , funcVal , 1e-4 , "" );
    coords_tmp = coords;
    gridVal = combigrid1->eval(coords_tmp);
    COMBIGRID_ERROR_TEST_EQUAL( gridVal , funcVal , 1e-4 , "" );
    coords_tmp = coords;

    coords[0] = min[0] + (max[0] - min[0]) * 0.0;
    coords[1] = min[1] + (max[1] - min[1]) * 0.0;
    coords[2] = min[2] + (max[2] - min[2]) * 0.0;
    coords_tmp = coords;
    funcVal = testfunction3D1( coords[0] , coords[1] , coords[2]);
    gridVal = combigrid->eval(coords_tmp);
    COMBIGRID_ERROR_TEST_EQUAL( gridVal , funcVal , 1e-4 , "" );
    coords_tmp = coords;
    //gridVal = combigrid1->eval(coords_tmp); COMBIGRID_ERROR_TEST_EQUAL( gridVal , funcVal , 1e-4 , "" );coords_tmp = coords;

    coords[0] = min[0] + (max[0] - min[0]) * 0.7999;
    coords[1] = min[1] + (max[1] - min[1]) * 0.1201;
    coords[2] = min[2] + (max[2] - min[2]) * 0.899;
    coords_tmp = coords;
    funcVal = testfunction3D1( coords[0] , coords[1] , coords[2]);
    gridVal = combigrid->eval(coords_tmp);
    COMBIGRID_ERROR_TEST_EQUAL( gridVal , funcVal , 5e-4 , "" );
    coords_tmp = coords;
    gridVal = combigrid1->eval(coords_tmp);
    COMBIGRID_ERROR_TEST_EQUAL( gridVal , funcVal , 5e-4 , "" );
    coords_tmp = coords;

    coords[0] = min[0] + (max[0] - min[0]) * 0.6099;
    coords[1] = min[1] + (max[1] - min[1]) * 0.1981;
    coords[2] = min[2] + (max[2] - min[2]) * 0.83;
    coords_tmp = coords;
    funcVal = testfunction3D1( coords[0] , coords[1] , coords[2]);
    gridVal = combigrid->eval(coords_tmp);
    COMBIGRID_ERROR_TEST_EQUAL( gridVal , funcVal , 5e-4 , "" );
    coords_tmp = coords;
    gridVal = combigrid1->eval(coords_tmp);
    COMBIGRID_ERROR_TEST_EQUAL( gridVal , funcVal , 5e-4 , "" );
    coords_tmp = coords;

    coords[0] = min[0] + (max[0] - min[0]) * 0.999;
    coords[1] = min[1] + (max[1] - min[1]) * 0.95;
    coords[2] = min[2] + (max[2] - min[2]) * 0.95;
    coords_tmp = coords;
    funcVal = testfunction3D1( coords[0] , coords[1] , coords[2]);
    gridVal = combigrid->eval(coords_tmp);
    COMBIGRID_ERROR_TEST_EQUAL( gridVal , funcVal , 5e-4 , "" );
    coords_tmp = coords;
    //gridVal = combigrid1->eval(coords_tmp); COMBIGRID_ERROR_TEST_EQUAL( gridVal , funcVal , 1e-4 , "" );coords_tmp = coords;

    coords[0] = min[0] + (max[0] - min[0]) * 1.0;
    coords[1] = min[1] + (max[1] - min[1]) * 1.0;
    coords[2] = min[2] + (max[2] - min[2]) * 1.0;
    coords_tmp = coords;
    funcVal = testfunction3D1( coords[0] , coords[1] , coords[2]);
    gridVal = combigrid->eval(coords_tmp);
    COMBIGRID_ERROR_TEST_EQUAL( gridVal , funcVal , 5e-4 , "" );
    coords_tmp = coords;
    //gridVal = combigrid1->eval(coords_tmp); COMBIGRID_ERROR_TEST_EQUAL( gridVal , funcVal , 1e-4 , "" );coords_tmp = coords;

    delete stretchingMaker;
    delete domain;
    delete combigrid;
    delete combigrid1;
    delete combischeme;

  }

};
}

#endif /* TESTTRETCHINGCOMBI_HPP_ */