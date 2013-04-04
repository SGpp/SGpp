/* ****************************************************************************
* Copyright (C) 2011 Technische Universitaet Muenchen                         *
* This file is part of the SG++ project. For conditions of distribution and   *
* use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
**************************************************************************** */
// @author Janos Benk (benk@in.tum.de)

#ifndef TESTTIKHONOVPROBLEM_HPP_
#define TESTTIKHONOVPROBLEM_HPP_

#include "combigrid.hpp"
#include "combigrid/multigridFG/utils/RunTikhonov.hpp"
#include "combigrid/plotter/GridPlotter.hpp"

namespace combigrid {

  /** Class to thes the all the full grids functionalities */
  class TestTikhonov {
    public:

      TestTikhonov();

      static void test_all_cases() {
        COMBIGRID_OUT_LEVEL3( 4 , "Tikhonov testing ... ");
        // 1D tests
        test1();
        test12();
        test13();
        test14();

        // 2D tests
        test2();
        test201();
        test21();
        test22();
        test2_cross();
        test2_cross_CG();

        // 3D test
        test3();
        test31();
      }

      static double testfunction2D( double x1, double x2) {
        return (x1 + x2);
      }

      static void test1() {

        std::vector<double> min(1);
        min[0] = 0.8;
        std::vector<double> max(1);
        max[0] = 1.2;
        std::vector<int> levels(1);
        levels[0] = 2;

        AbstractStretchingMaker* stretchingMaker = new UniformStretching();
        GridDomain* domain = new GridDomain( 1 , levels , min , max , (*stretchingMaker) );

        // make the Tikhonov computations
        FullGridD* fg = RunTikhonov::computeFGTikhonov( *domain , levels , 0.01 ,
                        "../../../datasets/regression/1DLine/X.txt" , "../../../datasets/regression/1DLine/Y.txt");

        // test directly the unknowns
        COMBIGRID_ERROR_TEST_EQUAL( fg->getElementVector()[0] , 0.8 , 1e-7 , "" );
        COMBIGRID_ERROR_TEST_EQUAL( fg->getElementVector()[1] , 0.9 , 1e-7 , "" );
        COMBIGRID_ERROR_TEST_EQUAL( fg->getElementVector()[2] , 1.0 , 1e-7 , "" );
        COMBIGRID_ERROR_TEST_EQUAL( fg->getElementVector()[3] , 1.1 , 1e-7 , "" );
        COMBIGRID_ERROR_TEST_EQUAL( fg->getElementVector()[4] , 1.2 , 1e-7 , "" );

        std::vector<double> globalC(1, 0.0);
        GridPlotter::plotFullGrid( "line1D.m" , fg , globalC );

        delete fg;
        delete domain;
        delete stretchingMaker;
      }

      static void test12() {

        std::vector<double> min(1);
        min[0] = 0.8;
        std::vector<double> max(1);
        max[0] = 1.2;
        std::vector<int> levels(1);
        levels[0] = 5;

        AbstractStretchingMaker* stretchingMaker = new AtanSpecialStretching(); //new UniformStretching();
        GridDomain* domain = new GridDomain( 1 , levels , min , max , (*stretchingMaker) );

        // make the Tikhonov computations
        FullGridD* fg = RunTikhonov::computeFGTikhonov( *domain , levels , 3e-5 ,
                        "../../../datasets/regression/1DOption/X.txt" , "../../../datasets/regression/1DOption/Y.txt");

        // test directly the unknowns
        std::vector<double> coord(1);
        coord[0] = 0.81;
        COMBIGRID_ERROR_TEST_EQUAL( fg->eval(coord) , 0.05 , 1e-1 , "" );

        std::vector<double> globalC(1, 0.0);
        GridPlotter::plotFullGrid( "option1D.m" , fg , globalC );

        delete fg;
        delete domain;
        delete stretchingMaker;
      }

      static void test13() {

        std::vector<double> min(1);
        min[0] = 0.6;
        std::vector<double> max(1);
        max[0] = 1.4;
        std::vector<int> levels(1);
        levels[0] = 5;

        AbstractStretchingMaker* stretchingMaker = new AtanSpecialStretching(); //new UniformStretching();
        GridDomain* domain = new GridDomain( 1 , levels , min , max , (*stretchingMaker) );

        // make the Tikhonov computations
        FullGridD* fg = RunTikhonov::computeFGTikhonov( *domain , levels , 0.01 ,
                        "../../../datasets/regression/1DLine/X.txt" , "../../../datasets/regression/1DLine/Y.txt");

        // test directly the unknowns
        COMBIGRID_ERROR_TEST_EQUAL( fg->getElementVector()[0] , 0.6 , 1e-7 , "" );
        COMBIGRID_ERROR_TEST_EQUAL( fg->getElementVector()[4] , 1.4 , 1e-7 , "" );

        delete fg;
        delete domain;
        delete stretchingMaker;
      }

      static void test14() {

        std::vector<double> min(1);
        min[0] = 0.5;
        std::vector<double> max(1);
        max[0] = 1.5;
        std::vector<int> levels(1);
        levels[0] = 5;

        AbstractStretchingMaker* stretchingMaker = new TanStretching(0.5); //new AtanSpecialStretching(); //new UniformStretching(); //new AtanSpecialStretching();
        GridDomain* domain = new GridDomain( 1 , levels , min , max , (*stretchingMaker) );

        // make the Tikhonov computations
        FullGridD* fg = RunTikhonov::computeFGTikhonov( *domain , levels , 3e-5 ,
                        "../../../datasets/regression/1DOption_1Y/X.txt" , "../../../datasets/regression/1DOption_1Y/Y.txt");

        // test directly the unknowns
        std::vector<double> coord(1);
        coord[0] = 0.81;
        COMBIGRID_ERROR_TEST_EQUAL( fg->eval(coord) , 0.05 , 1e-1 , "" );

        std::vector<double> globalC(1, 0.0);
        GridPlotter::plotFullGrid( "option1D_1Y.m" , fg , globalC );

        delete fg;
        delete domain;
        delete stretchingMaker;
      }

      static void test2() {

        std::vector<double> min(2);
        min[0] = 0.9;
        min[1] = 0.9;
        std::vector<double> max(2);
        max[0] = 1.1;
        max[1] = 1.1;
        std::vector<int> levels(2);
        levels[0] = 3;
        levels[1] = 3;

        AbstractStretchingMaker* stretchingMaker = new UniformStretching();
        GridDomain* domain = new GridDomain( 2 , levels , min , max , (*stretchingMaker) );

        // make the Tikhonov computations
        FullGridD* fg = RunTikhonov::computeFGTikhonov( *domain , levels , 0.001 ,
                        "../../../datasets/regression/2DPlane/X.txt" , "../../../datasets/regression/2DPlane/Y.txt");

        // test the unkowns
        COMBIGRID_ERROR_TEST_EQUAL( fg->getElementVector()[0] , 1.8 , 1e-2 , "" );
        COMBIGRID_ERROR_TEST_EQUAL( fg->getElementVector()[fg->getElementVector().size() - 1] , 2.2 , 1e-5 , "" );

        std::vector<double> globalC(2, 0.0);
        GridPlotter::plotFullGrid( "plane2D.m" , fg , globalC );

        delete fg;
        delete domain;
        delete stretchingMaker;
      }

      static void test201() {

        std::vector<double> min(2);
        min[0] = 0.9;
        min[1] = 0.9;
        std::vector<double> max(2);
        max[0] = 1.1;
        max[1] = 1.1;
        std::vector<int> levels(2);
        levels[0] = 4;
        levels[1] = 3;

        AbstractStretchingMaker* stretchingMaker = new AtanSpecialStretching();
        GridDomain* domain = new GridDomain( 2 , levels , min , max , (*stretchingMaker) );

        // make the Tikhonov computations
        FullGridD* fg = RunTikhonov::computeFGTikhonov( *domain , levels , 0.001 ,
                        "../../../datasets/regression/2DPlane/X.txt" , "../../../datasets/regression/2DPlane/Y.txt");

        // test the unkowns
        COMBIGRID_ERROR_TEST_EQUAL( fg->getElementVector()[0] , 1.8 , 1e-4 , "" );
        COMBIGRID_ERROR_TEST_EQUAL( fg->getElementVector()[fg->getElementVector().size() - 1] , 2.2 , 1e-4 , "" );

        delete fg;
        delete domain;
        delete stretchingMaker;
      }

      static void test21() {

        std::vector<double> min(2);
        min[0] = 0.9;
        min[1] = 0.9;
        std::vector<double> max(2);
        max[0] = 1.1;
        max[1] = 1.1;
        std::vector<int> levels(2);
        levels[0] = 4;
        levels[1] = 3;

        AbstractStretchingMaker* stretchingMaker = new AtanSpecialStretching();
        GridDomain* domain = new GridDomain( 2 , levels , min , max , (*stretchingMaker) );

        // make the Tikhonov computations
        FullGridD* fg = RunTikhonov::computeFGTikhonov( *domain , levels , 2e-5 ,
                        "../../../datasets/regression/2DOption/X.txt" , "../../../datasets/regression/2DOption/Y.txt");

        std::vector<double> coord(2);
        coord[0] = 0.81;
        coord[1] = 0.81;
        COMBIGRID_ERROR_TEST_EQUAL( fg->eval(coord) , 0.05 , 1e-1 , "" );

        std::vector<double> globalC(2, 0.0);
        GridPlotter::plotFullGrid( "option2D.m" , fg , globalC );

        delete fg;
        delete domain;
        delete stretchingMaker;
      }

      static void test22() {

        std::vector<double> min(2);
        min[0] = 0.6;
        min[1] = 0.6;
        std::vector<double> max(2);
        max[0] = 1.4;
        max[1] = 1.4;
        std::vector<int> levels(2);
        levels[0] = 4;
        levels[1] = 4;

        AbstractStretchingMaker* stretchingMaker = new AtanSpecialStretching();
        GridDomain* domain = new GridDomain( 2 , levels , min , max , (*stretchingMaker) );

        // make the Tikhonov computations
        FullGridD* fg = RunTikhonov::computeFGTikhonov( *domain , levels , 2e-5 ,
                        "../../../datasets/regression/2DOption_1Y/X.txt" , "../../../datasets/regression/2DOption_1Y/Y.txt");

        std::vector<double> coord(2);
        coord[0] = 0.81;
        coord[1] = 0.81;
        COMBIGRID_ERROR_TEST_EQUAL( fg->eval(coord) , 0.05 , 1e-1 , "" );

        std::vector<double> globalC(2, 0.0);
        GridPlotter::plotFullGrid( "option2D_1Y.m" , fg , globalC );

        delete fg;
        delete domain;
        delete stretchingMaker;
      }

      static void test2_cross() {
        std::vector<double> min(2);
        min[0] = 0.6;
        min[1] = 0.6;
        std::vector<double> max(2);
        max[0] = 1.4;
        max[1] = 1.4;
        std::vector<int> levels(2);
        levels[0] = 3;
        levels[1] = 4;
        int dim , nrPoints;
        double lambda;

        AbstractStretchingMaker* stretchingMaker = new AtanSpecialStretching();
        GridDomain* domain = new GridDomain( 2 , levels , min , max , (*stretchingMaker) );
        FullGridD* fg = new FullGridD( 2 , levels );
        fg->createFullGrid();
        fg->setDomain( domain );

        std::vector<double> XCoords;
        std::vector<double> YPoint;
        RunTikhonov::readInInput( "../../../datasets/regression/2DOption_1Y/X.txt" , "../../../datasets/regression/2DOption_1Y/Y.txt" ,
                                  dim , nrPoints , XCoords, YPoint);

        // make the Tikhonov computations
        RunTikhonov::computeTikhonov_FG_crossvalidation( fg, lambda, XCoords, YPoint);

        std::vector<double> globalC(2, 0.0);
        GridPlotter::plotFullGrid( "option2D_cross.m" , fg , globalC , 100 );

        delete fg;
        delete domain;
        delete stretchingMaker;
      }

      static void test2_cross_CG() {
        std::vector<double> min(2);
        min[0] = 0.6;
        min[1] = 0.6;
        std::vector<double> max(2);
        max[0] = 1.4;
        max[1] = 1.4;
        std::vector<int> minlevels(2);
        minlevels[0] = 2;
        minlevels[1] = 2;
        std::vector<int> maxlevels(2);
        maxlevels[0] = 4;
        maxlevels[1] = 4;
        int dim , nrPoints;
        double lambda;

        AbstractStretchingMaker* stretchingMaker = new AtanSpecialStretching();
        GridDomain* domain = new GridDomain( 2 , maxlevels , min , max , (*stretchingMaker) );

        CombiSchemeBasis* combischeme = new TS_CT( minlevels, maxlevels);
        AbstractCombiGrid* combigrid = new SerialCombiGrid(combischeme);
        // creates the full grids in the combination grid
        combigrid->createFullGrids();
        combigrid->setDomainAllFG( domain );

        std::vector<double> XCoords;
        std::vector<double> YPoint;
        RunTikhonov::readInInput( "../../../datasets/regression/2DOption_1Y/X.txt" , "../../../datasets/regression/2DOption_1Y/Y.txt" ,
                                  dim , nrPoints , XCoords, YPoint);

        // make the Tikhonov computations
        RunTikhonov::computeTikhonov_CG_crossvalidation( combigrid , lambda, XCoords, YPoint);

        std::vector<double> globalC(2, 0.0);
        combigrid->setDomainAllFG( NULL );
        combigrid->setDomain( domain );
        GridPlotter::plotCombiGrid( "option2D_CG_cross.m" , combigrid , globalC , 100 );

        delete stretchingMaker;
        delete domain;
        delete combigrid;
        delete combischeme;
      }

      static void test3() {

        std::vector<double> min(3);
        min[0] = 0.9;
        min[1] = 0.9;
        min[2] = 0.9;
        std::vector<double> max(3);
        max[0] = 1.1;
        max[1] = 1.1;
        max[2] = 1.1;
        std::vector<int> levels(3);
        levels[0] = 2;
        levels[1] = 4;
        levels[2] = 3;

        //AbstractStretchingMaker* stretchingMaker = new UniformStretching();
        AbstractStretchingMaker* stretchingMaker = new AtanSpecialStretching();
        GridDomain* domain = new GridDomain( 3 , levels , min , max , (*stretchingMaker) );

        // make the Tikhonov computations
        FullGridD* fg = RunTikhonov::computeFGTikhonov( *domain , levels , 0.0001 ,
                        "../../../datasets/regression/3DPlane/X.txt" , "../../../datasets/regression/3DPlane/Y.txt");

        std::vector<double> globalC1(3, 1.0);
        std::vector<double> globalC2(3, 0.97);
        double tmp;
        // test the unkowns
        tmp = fg->eval( globalC1 );
        COMBIGRID_ERROR_TEST_EQUAL( tmp , 3.0 , 1e-6 , "" );
        tmp = fg->eval( globalC2 );
        COMBIGRID_ERROR_TEST_EQUAL( tmp , 2.91 , 1e-6 , "" );


        delete fg;
        delete domain;
        delete stretchingMaker;
      }

      static void test31() {

        std::vector<double> min(3);
        min[0] = 0.9;
        min[1] = 0.9;
        min[2] = 0.9;
        std::vector<double> max(3);
        max[0] = 1.1;
        max[1] = 1.1;
        max[2] = 1.1;
        std::vector<int> levels(3);
        levels[0] = 2;
        levels[1] = 4;
        levels[2] = 3;

        //AbstractStretchingMaker* stretchingMaker = new UniformStretching();
        AbstractStretchingMaker* stretchingMaker = new AtanSpecialStretching();
        GridDomain* domain = new GridDomain( 3 , levels , min , max , (*stretchingMaker) );

        // make the Tikhonov computations
        FullGridD* fg = RunTikhonov::computeFGTikhonov( *domain , levels , 1e-5 ,
                        "../../../datasets/regression/3DOption/X.txt" , "../../../datasets/regression/3DOption/Y.txt");

        std::vector<double> globalC1(3, 1.0);
        double tmp;
        // test the unkowns
        tmp = fg->eval( globalC1 );
        COMBIGRID_ERROR_TEST_EQUAL( tmp , 0.1 , 1e-1 , "" );

        delete fg;
        delete domain;
        delete stretchingMaker;
      }
  };
}

#endif /* TESTTIKHONOVPROBLEM_HPP_ */
