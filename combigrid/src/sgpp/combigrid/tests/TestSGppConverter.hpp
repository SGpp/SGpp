// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#ifndef TESTSGPPCONVERTER_HPP_
#define TESTSGPPCONVERTER_HPP_

#include "combigrid.hpp"

namespace combigrid {

  /** Class to thes the all the full grids functionalities */
  class TestSGppConverter {
    public:

      TestSGppConverter();

      static void test_all_cases() {
        COMBIGRID_OUT_LEVEL3( 4 , "SGppConverter ... ");
        // 2D test
        test1();

        // 3D test
        test2();
        test21(); // S-CT
        test22(); // S-CT without boundary
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
        CombiSchemeBasis* combischeme = new TS_CT( 2 , 6 );
        AbstractCombiGrid* combigrid = new SerialCombiGrid(combischeme);
        SGPP::base::GridStorage* gridstorageSGpp;
        SGPP::base::DataVector* alphas , *minAlpha , *maxAlpha ;

        // creates the full grids in the combination grid
        combigrid->createFullGrids();

        // for each full grid set the values
        for ( int fg = 0 ; fg < combigrid->getNrFullGrid() ; fg++) {
          FullGridD* fgp = combigrid->getFullGrid(fg);

          // for this full grid set the values
          for (int i = 0 ; i < fgp->getNrElements() ; i++) {
            fgp->getCoords( i , coords );
            fgp->getElementVector()[i] = testfunction2D( coords[0] , coords[1] );
          }
        }

        // test the number of subspaces
        COMBIGRID_ERROR_TEST_EQUAL( combigrid->getNrFullGrid() , 3.0 , 1e-10 , "" );

        // create the grid storage
        gridstorageSGpp = combigrid->createSGppGridStorage();

        // create the coefficient vectors for the SGpp grid
        alphas = new SGPP::base::DataVector(gridstorageSGpp->size());
        minAlpha = new SGPP::base::DataVector(gridstorageSGpp->size());
        maxAlpha = new SGPP::base::DataVector(gridstorageSGpp->size());

        // recompose full grids -> SGpp
        combigrid->reCompose( gridstorageSGpp , alphas , minAlpha , maxAlpha );
        combigrid->reCompose( gridstorageSGpp , alphas );

        double minA, maxA, minMin , maxMax;
        maxMax = (*maxAlpha)[0];
        minMin = (*minAlpha)[0];
        minA = maxA = (*alphas)[0];

        for (int ind = 1 ; ind < (int)gridstorageSGpp->size() ; ind++) {
          maxMax = ( (*maxAlpha)[ind] > maxMax) ? (*maxAlpha)[ind] : maxMax;
          minMin = ( (*minAlpha)[ind] < minMin) ? (*minAlpha)[ind] : minMin;
          minA = ( (*minAlpha)[ind] < minA) ? (*minAlpha)[ind] : minA;
          maxA = ( (*minAlpha)[ind] > maxA) ? (*minAlpha)[ind] : maxA;
        }

        coords[0] = 0.0;
        coords[1] = 0.0;
        COMBIGRID_ERROR_TEST_EQUAL( minMin , testfunction2D( coords[0] , coords[1] ) , 1e-10 , "" );
        COMBIGRID_ERROR_TEST_EQUAL( minA , testfunction2D( coords[0] , coords[1] ) , 1e-10 , "" );
        coords[0] = 1.0;
        coords[1] = 1.0;
        COMBIGRID_ERROR_TEST_EQUAL( maxMax , testfunction2D( coords[0] , coords[1] ) , 1e-10 , "" );
        COMBIGRID_ERROR_TEST_EQUAL( maxA , testfunction2D( coords[0] , coords[1] ) , 1e-10 , "" );

        // for each full grid set the values
        for ( int fg = 0 ; fg < combigrid->getNrFullGrid() ; fg++) {
          FullGridD* fgp = combigrid->getFullGrid(fg);

          for (int i = 0 ; i < fgp->getNrElements() ; i++) {
            // set all the full grid values with dummy values
            // so that we can test if the decomposition really works
            fgp->getElementVector()[i] = 912.34;
          }
        }

        // --- decompose SGpp -> full grids
        combigrid->deCompose( gridstorageSGpp , alphas );

        // test if the combi grid has the correct values after the decomposition
        coords[0] = 1.0 / 3.0;
        coords[1] = 1.0 / 3.0;
        COMBIGRID_ERROR_TEST_EQUAL( combigrid->eval(coords) , testfunction2D( coords[0] , coords[1] ) , 1e-10 , "" );
        coords[0] = 0.0001;
        coords[1] = 1.0 / 3.0;
        COMBIGRID_ERROR_TEST_EQUAL( combigrid->eval(coords) , testfunction2D( coords[0] , coords[1] ) , 1e-10 , "" );
        coords[0] = 1.0 / 3.0;
        coords[1] = 0.0001;
        COMBIGRID_ERROR_TEST_EQUAL( combigrid->eval(coords) , testfunction2D( coords[0] , coords[1] ) , 1e-10 , "" );
        coords[0] = 0.00001;
        coords[1] = 0.00001;
        COMBIGRID_ERROR_TEST_EQUAL( combigrid->eval(coords) , testfunction2D( coords[0] , coords[1] ) , 1e-10 , "" );
        coords[0] = 1.0;
        coords[1] = 1.0;
        COMBIGRID_ERROR_TEST_EQUAL( combigrid->eval(coords) , testfunction2D( coords[0] , coords[1] ) , 1e-10 , "" );
        coords[0] = 0.999;
        coords[1] = 0.95;
        COMBIGRID_ERROR_TEST_EQUAL( combigrid->eval(coords) , testfunction2D( coords[0] , coords[1] ) , 1e-10 , "" );
        coords[0] = 0.999;
        coords[1] = 0.001;
        COMBIGRID_ERROR_TEST_EQUAL( combigrid->eval(coords) , testfunction2D( coords[0] , coords[1] ) , 1e-10 , "" );

        // delete the allocated variables
        delete alphas;
        delete minAlpha;
        delete maxAlpha;
        delete gridstorageSGpp;
        delete combigrid;
        delete combischeme;
      }


      static void test2() {

        std::vector<double> coords( 3 , 0.0);
        std::vector<int> levels( 3 , 7);
        std::vector<bool> makeCombi( 3 , true);
        CombiSchemeBasis* combischeme = new TS_CT( 3 , levels , makeCombi );
        AbstractCombiGrid* combigrid = new SerialCombiGrid(combischeme);
        SGPP::base::GridStorage* gridstorageSGpp;
        SGPP::base::DataVector* alphas , *minAlpha , *maxAlpha ;

        // creates the full grids in the combination grid
        combigrid->createFullGrids();

        // for each full grid set the values
        for ( int fg = 0 ; fg < combigrid->getNrFullGrid() ; fg++) {
          FullGridD* fgp = combigrid->getFullGrid(fg);

          // for this full grid set the values
          for (int i = 0 ; i < fgp->getNrElements() ; i++) {
            fgp->getCoords( i , coords );
            fgp->getElementVector()[i] = testfunction3D1( coords[0] , coords[1] , coords[2]);
          }
        }

        COMBIGRID_ERROR_TEST_EQUAL( combigrid->getNrFullGrid() , 4.0 , 1e-10 , "" );

        // create the grid storage
        gridstorageSGpp = combigrid->createSGppGridStorage();

        // create the coefficient vectors for the SGpp grid
        alphas = new SGPP::base::DataVector(gridstorageSGpp->size());
        minAlpha = new SGPP::base::DataVector(gridstorageSGpp->size());
        maxAlpha = new SGPP::base::DataVector(gridstorageSGpp->size());

        // recompose full grids -> SGpp
        combigrid->reCompose( gridstorageSGpp , alphas , minAlpha , maxAlpha );
        combigrid->reCompose( gridstorageSGpp , alphas );

        double minA, maxA, minMin , maxMax;
        maxMax = (*maxAlpha)[0];
        minMin = (*minAlpha)[0];
        minA = maxA = (*alphas)[0];

        for (int ind = 1 ; ind < (int)gridstorageSGpp->size() ; ind++) {
          maxMax = ( (*maxAlpha)[ind] > maxMax) ? (*maxAlpha)[ind] : maxMax;
          minMin = ( (*minAlpha)[ind] < minMin) ? (*minAlpha)[ind] : minMin;
          minA = ( (*minAlpha)[ind] < minA) ? (*minAlpha)[ind] : minA;
          maxA = ( (*minAlpha)[ind] > maxA) ? (*minAlpha)[ind] : maxA;
        }

        coords[0] = 0.0;
        coords[1] = 0.0;
        coords[2] = 0.0;
        COMBIGRID_ERROR_TEST_EQUAL( minMin , testfunction3D1( coords[0] , coords[1] , coords[2]) , 1e-10 , "" );
        COMBIGRID_ERROR_TEST_EQUAL( minA , testfunction3D1( coords[0] , coords[1] , coords[2]) , 1e-10 , "" );
        coords[0] = 1.0;
        coords[1] = 1.0;
        coords[2] = 1.0;
        COMBIGRID_ERROR_TEST_EQUAL( maxMax , testfunction3D1( coords[0] , coords[1] , coords[2]) , 1e-10 , "" );
        COMBIGRID_ERROR_TEST_EQUAL( maxA , testfunction3D1( coords[0] , coords[1] , coords[2]) , 1e-10 , "" );

        // for each full grid set the values
        for ( int fg = 0 ; fg < combigrid->getNrFullGrid() ; fg++) {
          FullGridD* fgp = combigrid->getFullGrid(fg);

          for (int i = 0 ; i < fgp->getNrElements() ; i++) {
            // set all the full grid values with dummy values
            // so that we can test if the decomposition really works
            fgp->getElementVector()[i] = 912.34;
          }
        }

        // --- decompose SGpp -> full grids
        combigrid->deCompose( gridstorageSGpp , alphas );

        coords[0] = 1.0 / 3.0;
        coords[1] = 1.0 / 3.0;
        coords[2] = 1.0 / 3.0;
        COMBIGRID_ERROR_TEST_EQUAL( combigrid->eval(coords) , testfunction3D1( coords[0] , coords[1] , coords[2] ) , 1e-4 , "" );
        coords[0] = 0.0001;
        coords[1] = 1.0 / 3.0;
        coords[2] = 1.0 / 3.0;
        COMBIGRID_ERROR_TEST_EQUAL( combigrid->eval(coords) , testfunction3D1( coords[0] , coords[1] , coords[2] ) , 1e-4 , "" );
        coords[0] = 1.0 / 3.0;
        coords[1] = 0.0001;
        coords[2] = 0.0001;
        COMBIGRID_ERROR_TEST_EQUAL( combigrid->eval(coords) , testfunction3D1( coords[0] , coords[1] , coords[2] ) , 1e-4 , "" );
        coords[0] = 0.0;
        coords[1] = 0.0;
        coords[2] = 0.0;
        COMBIGRID_ERROR_TEST_EQUAL( combigrid->eval(coords) , testfunction3D1( coords[0] , coords[1] , coords[2] ) , 1e-4 , "" );
        coords[0] = 1.0;
        coords[1] = 1.0;
        coords[2] = 1.0;
        COMBIGRID_ERROR_TEST_EQUAL( combigrid->eval(coords) , testfunction3D1( coords[0] , coords[1] , coords[2] ) , 1e-4 , "" );
        coords[0] = 0.999;
        coords[1] = 0.95;
        coords[2] = 0.95;
        COMBIGRID_ERROR_TEST_EQUAL( combigrid->eval(coords) , testfunction3D1( coords[0] , coords[1] , coords[2] ) , 1e-4 , "" );
        coords[0] = 0.999;
        coords[1] = 0.001;
        coords[2] = 0.999;
        COMBIGRID_ERROR_TEST_EQUAL( combigrid->eval(coords) , testfunction3D1( coords[0] , coords[1] , coords[2] ) , 1e-4 , "" );

        // delete the allocated variables
        delete alphas;
        delete minAlpha;
        delete maxAlpha;
        delete gridstorageSGpp;
        delete combigrid;
        delete combischeme;
      }


      static void test21() {

        std::vector<double> coords( 3 , 0.0);
        std::vector<int> levels( 3 , 9);
        std::vector<int> levels_trunc( 3 , 3);
        levels[0] = 9;
        levels[1] = 6;
        levels[2] = 8;
        levels_trunc[0] = 1;
        levels_trunc[1] = 3;
        levels_trunc[2] = 1;
        CombiSchemeBasis* combischeme = new S_CT( 3 , levels , levels_trunc );
        AbstractCombiGrid* combigrid = new SerialCombiGrid(combischeme);
        SGPP::base::GridStorage* gridstorageSGpp;
        SGPP::base::DataVector* alphas , *minAlpha , *maxAlpha ;

        // creates the full grids in the combination grid
        combigrid->createFullGrids();

        // for each full grid set the values
        for ( int fg = 0 ; fg < combigrid->getNrFullGrid() ; fg++) {
          FullGridD* fgp = combigrid->getFullGrid(fg);

          // for this full grid set the values
          for (int i = 0 ; i < fgp->getNrElements() ; i++) {
            fgp->getCoords( i , coords );
            fgp->getElementVector()[i] = testfunction3D1( coords[0] , coords[1] , coords[2]);
          }
        }

        // create the grid storage
        gridstorageSGpp = combigrid->createSGppGridStorage();

        // create the coefficient vectors for the SGpp grid
        alphas = new SGPP::base::DataVector(gridstorageSGpp->size());
        minAlpha = new SGPP::base::DataVector(gridstorageSGpp->size());
        maxAlpha = new SGPP::base::DataVector(gridstorageSGpp->size());

        // recompose full grids -> SGpp
        combigrid->reCompose( gridstorageSGpp , alphas , minAlpha , maxAlpha );
        combigrid->reCompose( gridstorageSGpp , alphas );

        double minA, maxA, minMin , maxMax;
        maxMax = (*maxAlpha)[0];
        minMin = (*minAlpha)[0];
        minA = maxA = (*alphas)[0];

        for (int ind = 1 ; ind < (int)gridstorageSGpp->size() ; ind++) {
          maxMax = ( (*maxAlpha)[ind] > maxMax) ? (*maxAlpha)[ind] : maxMax;
          minMin = ( (*minAlpha)[ind] < minMin) ? (*minAlpha)[ind] : minMin;
          minA = ( (*minAlpha)[ind] < minA) ? (*minAlpha)[ind] : minA;
          maxA = ( (*minAlpha)[ind] > maxA) ? (*minAlpha)[ind] : maxA;
        }

        coords[0] = 0.0;
        coords[1] = 0.0;
        coords[2] = 0.0;
        COMBIGRID_ERROR_TEST_EQUAL( minMin , testfunction3D1( coords[0] , coords[1] , coords[2]) , 1e-10 , "" );
        COMBIGRID_ERROR_TEST_EQUAL( minA , testfunction3D1( coords[0] , coords[1] , coords[2]) , 1e-10 , "" );
        coords[0] = 1.0;
        coords[1] = 1.0;
        coords[2] = 1.0;
        COMBIGRID_ERROR_TEST_EQUAL( maxMax , testfunction3D1( coords[0] , coords[1] , coords[2]) , 1e-10 , "" );
        COMBIGRID_ERROR_TEST_EQUAL( maxA , testfunction3D1( coords[0] , coords[1] , coords[2]) , 1e-10 , "" );

        // for each full grid set the values
        for ( int fg = 0 ; fg < combigrid->getNrFullGrid() ; fg++) {
          FullGridD* fgp = combigrid->getFullGrid(fg);

          for (int i = 0 ; i < fgp->getNrElements() ; i++) {
            // set all the full grid values with dummy values
            // so that we can test if the decomposition really works
            fgp->getElementVector()[i] = 912.34;
          }
        }

        // --- decompose SGpp -> full grids
        combigrid->deCompose( gridstorageSGpp , alphas );

        coords[0] = 1.0 / 3.0;
        coords[1] = 1.0 / 3.0;
        coords[2] = 1.0 / 3.0;
        COMBIGRID_ERROR_TEST_EQUAL( combigrid->eval(coords) , testfunction3D1( coords[0] , coords[1] , coords[2] ) , 1e-3 , "" );
        coords[0] = 0.0001;
        coords[1] = 1.0 / 3.0;
        coords[2] = 1.0 / 3.0;
        COMBIGRID_ERROR_TEST_EQUAL( combigrid->eval(coords) , testfunction3D1( coords[0] , coords[1] , coords[2] ) , 1e-3 , "" );
        coords[0] = 1.0 / 3.0;
        coords[1] = 0.0001;
        coords[2] = 0.0001;
        COMBIGRID_ERROR_TEST_EQUAL( combigrid->eval(coords) , testfunction3D1( coords[0] , coords[1] , coords[2] ) , 1e-3 , "" );
        coords[0] = 0.0;
        coords[1] = 0.0;
        coords[2] = 0.0;
        COMBIGRID_ERROR_TEST_EQUAL( combigrid->eval(coords) , testfunction3D1( coords[0] , coords[1] , coords[2] ) , 1e-3 , "" );
        coords[0] = 1.0;
        coords[1] = 1.0;
        coords[2] = 1.0;
        COMBIGRID_ERROR_TEST_EQUAL( combigrid->eval(coords) , testfunction3D1( coords[0] , coords[1] , coords[2] ) , 1e-3 , "" );
        coords[0] = 0.999;
        coords[1] = 0.95;
        coords[2] = 0.95;
        COMBIGRID_ERROR_TEST_EQUAL( combigrid->eval(coords) , testfunction3D1( coords[0] , coords[1] , coords[2] ) , 1e-3 , "" );
        coords[0] = 0.999;
        coords[1] = 0.001;
        coords[2] = 0.999;
        COMBIGRID_ERROR_TEST_EQUAL( combigrid->eval(coords) , testfunction3D1( coords[0] , coords[1] , coords[2] ) , 1e-3 , "" );

        // delete the allocated variables
        delete alphas;
        delete minAlpha;
        delete maxAlpha;
        delete gridstorageSGpp;
        delete combigrid;
        delete combischeme;
      }


      static void test22() {

        std::vector<double> coords( 3 , 0.0);
        std::vector<int> levels( 3 , 9);
        std::vector<int> levels_trunc( 3 , 3);
        levels[0] = 6;
        levels[1] = 5;
        levels[2] = 6;
        levels_trunc[0] = 4;
        levels_trunc[1] = 3;
        levels_trunc[2] = 3;
        CombiSchemeBasis* combischeme = new S_CT( 3 , levels , levels_trunc );
        AbstractCombiGrid* combigrid = new SerialCombiGrid( combischeme , false );
        SGPP::base::GridStorage* gridstorageSGpp;
        SGPP::base::DataVector* alphas , *minAlpha , *maxAlpha ;

        // creates the full grids in the combination grid
        combigrid->createFullGrids();

        // for each full grid set the values
        for ( int fg = 0 ; fg < combigrid->getNrFullGrid() ; fg++) {
          FullGridD* fgp = combigrid->getFullGrid(fg);

          // for this full grid set the values
          for (int i = 0 ; i < fgp->getNrElements() ; i++) {
            fgp->getCoords( i , coords );
            fgp->getElementVector()[i] = testfunction3D1( coords[0] , coords[1] , coords[2]);
          }
        }

        // create the grid storage
        gridstorageSGpp = combigrid->createSGppGridStorage();

        // create the coefficient vectors for the SGpp grid
        alphas = new SGPP::base::DataVector(gridstorageSGpp->size());
        minAlpha = new SGPP::base::DataVector(gridstorageSGpp->size());
        maxAlpha = new SGPP::base::DataVector(gridstorageSGpp->size());

        // recompose full grids -> SGpp
        combigrid->reCompose( gridstorageSGpp , alphas , minAlpha , maxAlpha );
        combigrid->reCompose( gridstorageSGpp , alphas );

        double minA, maxA, minMin , maxMax;
        maxMax = (*maxAlpha)[0];
        minMin = (*minAlpha)[0];
        minA = maxA = (*alphas)[0];

        for (int ind = 1 ; ind < (int)gridstorageSGpp->size() ; ind++) {
          maxMax = ( (*maxAlpha)[ind] > maxMax) ? (*maxAlpha)[ind] : maxMax;
          minMin = ( (*minAlpha)[ind] < minMin) ? (*minAlpha)[ind] : minMin;
          minA = ( (*minAlpha)[ind] < minA) ? (*minAlpha)[ind] : minA;
          maxA = ( (*minAlpha)[ind] > maxA) ? (*minAlpha)[ind] : maxA;
        }

        // since this grid has no boundary points the tolerance has to be higher
        coords[0] = 0.0;
        coords[1] = 0.0;
        coords[2] = 0.0;
        COMBIGRID_ERROR_TEST_EQUAL( minMin , testfunction3D1( coords[0] , coords[1] , coords[2]) , 1e-1 , "" );
        COMBIGRID_ERROR_TEST_EQUAL( minA , testfunction3D1( coords[0] , coords[1] , coords[2]) , 1e-1 , "" );
        coords[0] = 1.0;
        coords[1] = 1.0;
        coords[2] = 1.0;
        COMBIGRID_ERROR_TEST_EQUAL( maxMax , testfunction3D1( coords[0] , coords[1] , coords[2]) , 3e-0 , "" );
        COMBIGRID_ERROR_TEST_EQUAL( maxA , testfunction3D1( coords[0] , coords[1] , coords[2]) , 3e-0 , "" );

        // for each full grid set the values
        for ( int fg = 0 ; fg < combigrid->getNrFullGrid() ; fg++) {
          FullGridD* fgp = combigrid->getFullGrid(fg);

          for (int i = 0 ; i < fgp->getNrElements() ; i++) {
            // set all the full grid values with dummy values
            // so that we can test if the decomposition really works
            fgp->getElementVector()[i] = 912.34;
          }
        }

        // --- decompose SGpp -> full grids
        combigrid->deCompose( gridstorageSGpp , alphas );

        // since this grid has no boundary points the tolerance has to be higher
        coords[0] = 1.0 / 3.0;
        coords[1] = 1.0 / 3.0;
        coords[2] = 1.0 / 3.0;
        COMBIGRID_ERROR_TEST_EQUAL( combigrid->eval(coords) , testfunction3D1( coords[0] , coords[1] , coords[2] ) , 1e-3 , "" );
        coords[0] = 0.0001;
        coords[1] = 1.0 / 3.0;
        coords[2] = 1.0 / 3.0;
        COMBIGRID_ERROR_TEST_EQUAL( combigrid->eval(coords) , testfunction3D1( coords[0] , coords[1] , coords[2] ) , 1e-3 , "" );
        coords[0] = 1.0 / 3.0;
        coords[1] = 0.0001;
        coords[2] = 0.0001;
        COMBIGRID_ERROR_TEST_EQUAL( combigrid->eval(coords) , testfunction3D1( coords[0] , coords[1] , coords[2] ) , 1e-2 , "" );
        coords[0] = 0.0;
        coords[1] = 0.0;
        coords[2] = 0.0;
        COMBIGRID_ERROR_TEST_EQUAL( combigrid->eval(coords) , testfunction3D1( coords[0] , coords[1] , coords[2] ) , 1e-2 , "" );
        coords[0] = 1.0;
        coords[1] = 1.0;
        coords[2] = 1.0;
        COMBIGRID_ERROR_TEST_EQUAL( combigrid->eval(coords) , testfunction3D1( coords[0] , coords[1] , coords[2] ) , 3e-2 , "" );
        coords[0] = 0.999;
        coords[1] = 0.95;
        coords[2] = 0.95;
        COMBIGRID_ERROR_TEST_EQUAL( combigrid->eval(coords) , testfunction3D1( coords[0] , coords[1] , coords[2] ) , 3e-2 , "" );
        coords[0] = 0.999;
        coords[1] = 0.001;
        coords[2] = 0.999;
        COMBIGRID_ERROR_TEST_EQUAL( combigrid->eval(coords) , testfunction3D1( coords[0] , coords[1] , coords[2] ) , 3e-2 , "" );

        // delete the allocated variables
        delete alphas;
        delete minAlpha;
        delete maxAlpha;
        delete gridstorageSGpp;
        delete combigrid;
        delete combischeme;
      }
  };
}

#endif /* TESTSGPPCONVERTER_HPP_ */