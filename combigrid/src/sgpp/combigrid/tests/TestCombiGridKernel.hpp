// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#ifndef TESTCOMBIGRIDKERNEL_HPP_
#define TESTCOMBIGRIDKERNEL_HPP_

#include "combigrid.hpp"

namespace combigrid {

  /** Class to thes the all the full grids functionalities */
  class TestCombiGridKernel {
    public:

      TestCombiGridKernel();

      static void test_all_cases() {
        COMBIGRID_OUT_LEVEL3( 4 , "Testing CombiGrid Kernel ... ");
        // 3D test
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

        CombiGridKernel<double>* CG_kernel = new CombiGridKernel<double>(3);

        std::vector<int> levels(3, 5);
        levels[0] = 2;
        levels[1] = 6;
        levels[2] = 5;
        std::vector<bool> hasBdrPoints(3, true);

        std::vector<int> levels1(3, 6);
        levels1[0] = 4;
        levels1[1] = 6;
        levels1[2] = 5;
        std::vector<bool> hasBdrPoints1(3, false);

        CG_kernel->addFullGrid( levels , hasBdrPoints , 1.0 );
        CG_kernel->addFullGrid( levels , hasBdrPoints , 2.0 );
        CG_kernel->addFullGrid( levels , hasBdrPoints , 3.0 );
        CG_kernel->addFullGrid( levels1 , hasBdrPoints1 , 4.0 );

        COMBIGRID_ERROR_TEST_EQUAL( (double)CG_kernel->getNrFullGrids() , 4.0 , 1e-10 , "" );

        // delete the duplicated spaces
        CG_kernel->deleteDuplicate();

        // test the combi grid state after the duplicated spaces were deleted
        COMBIGRID_ERROR_TEST_EQUAL( (double)CG_kernel->getNrFullGrids() , 2.0 , 1e-10 , "" );
        COMBIGRID_ERROR_TEST_EQUAL( CG_kernel->getCoef(0) , 1.0 , 1e-10 , "" );
        COMBIGRID_ERROR_TEST_EQUAL( CG_kernel->getCoef(1) , 4.0 , 1e-10 , "" );

        // test the full grid
        FullGrid<double>* fg = CG_kernel->getFullGrid(1);
        COMBIGRID_ERROR_TEST_EQUAL( fg->getLevels()[0] , 4.0 , 1e-10 , "" );
        COMBIGRID_ERROR_TEST_EQUAL( fg->getLevels()[1] , 6.0 , 1e-10 , "" );
        COMBIGRID_ERROR_TEST_EQUAL( fg->getLevels()[2] , 5.0 , 1e-10 , "" );

        // delete the kernel which should delete all data
        delete CG_kernel;
      }


      static void test2() {
        CombiGridKernel<double>* CG_kernel = new CombiGridKernel<double>(3);

        std::vector<int> levels(3, 5);
        levels[0] = 2;
        levels[1] = 6;
        levels[2] = 5;
        std::vector<bool> hasBdrPoints(3, true);

        std::vector<int> levels1(3, 6);
        levels1[0] = 2;
        levels1[1] = 7;
        levels1[2] = 8;
        std::vector<bool> hasBdrPoints1(3, false);

        std::vector<int> levels2(3, 6);
        levels2[0] = 2;
        levels2[1] = 8;
        levels2[2] = 3;
        std::vector<bool> hasBdrPoints2(3, false);

        CG_kernel->addFullGrid( levels , hasBdrPoints , 1.0 );
        CG_kernel->addFullGrid( levels , hasBdrPoints , 1.0 );
        CG_kernel->addFullGrid( levels1 , hasBdrPoints1 , 3.0 );
        CG_kernel->addFullGrid( levels1 , hasBdrPoints1 , 3.0 );
        CG_kernel->addFullGrid( levels2 , hasBdrPoints2 , 6.0 );
        CG_kernel->addFullGrid( levels2 , hasBdrPoints2 , 6.0 );

        COMBIGRID_ERROR_TEST_EQUAL( (double)CG_kernel->getNrFullGrids() , 6.0 , 1e-10 , "" );

        // delete the duplicated spaces
        CG_kernel->deleteDuplicate();

        // test the combi grid state after the duplicated spaces were deleted
        COMBIGRID_ERROR_TEST_EQUAL( (double)CG_kernel->getNrFullGrids() , 3.0 , 1e-10 , "" );
        COMBIGRID_ERROR_TEST_EQUAL( CG_kernel->getCoef(0) , 1.0 , 1e-10 , "" );
        COMBIGRID_ERROR_TEST_EQUAL( CG_kernel->getCoef(1) , 6.0 , 1e-10 , "" );
        COMBIGRID_ERROR_TEST_EQUAL( CG_kernel->getCoef(2) , 3.0 , 1e-10 , "" );

        // test the full grid
        FullGrid<double>* fg = CG_kernel->getFullGrid(2);
        COMBIGRID_ERROR_TEST_EQUAL( fg->getLevels()[0] , 2.0 , 1e-10 , "" );
        COMBIGRID_ERROR_TEST_EQUAL( fg->getLevels()[1] , 7.0 , 1e-10 , "" );
        COMBIGRID_ERROR_TEST_EQUAL( fg->getLevels()[2] , 8.0 , 1e-10 , "" );

        fg = CG_kernel->getFullGrid(1);
        COMBIGRID_ERROR_TEST_EQUAL( fg->getLevels()[0] , 2.0 , 1e-10 , "" );
        COMBIGRID_ERROR_TEST_EQUAL( fg->getLevels()[1] , 8.0 , 1e-10 , "" );
        COMBIGRID_ERROR_TEST_EQUAL( fg->getLevels()[2] , 3.0 , 1e-10 , "" );

        // delete the kernel which should delete all data
        delete CG_kernel;
      }

  };
}

#endif /* TESTCOMBIGRIDKERNEL_HPP_ */