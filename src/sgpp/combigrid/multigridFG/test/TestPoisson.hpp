/* ****************************************************************************
* Copyright (C) 2011 Technische Universitaet Muenchen                         *
* This file is part of the SG++ project. For conditions of distribution and   *
* use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
**************************************************************************** */
// @author Janos Benk (benk@in.tum.de)

#ifndef TESTPOISSONPROBLEM_HPP_
#define TESTPOISSONPROBLEM_HPP_

#include "combigrid.hpp"
#include "combigrid/multigridFG/utils/RunPoisson.hpp"

namespace combigrid {

  /** Class to thes the all the full grids functionalities */
  class TestPoisson {
    public:

      TestPoisson();

      static void test_all_cases() {
        COMBIGRID_OUT_LEVEL3( 4 , "Poisson testing ... ");
        // 1D tests
        test1();

        // 2D tests
        test2();

        //3D test
        test3();

        //4D test
        test4();
      }

      static double testfunction2D( double x1, double x2) {
        return (x1 + x2);
      }

      static void test1() {

        std::vector<double> min(1);
        min[0] = 0.0;
        std::vector<double> max(1);
        max[0] = 1.0;
        std::vector<double> sigma(1);
        sigma[0] = 1.0;
        std::vector<int> levels(1);
        levels[0] = 7;

        AbstractStretchingMaker* stretchingMaker = new UniformStretching();
        GridDomain* domain = new GridDomain( 1 , levels , min , max , (*stretchingMaker) );
        ConstRHS* rhs = new ConstRHS(1.0);

        // make the Tikhonov computations
        FullGridD* fg = RunPoisson::computeFGPoisson( *domain , levels , sigma , 1.0 , rhs );

        // test directly the unknowns
        COMBIGRID_ERROR_TEST_EQUAL( fg->getElementVector()[0] , 0.0 , 1e-5 , "" );
        COMBIGRID_ERROR_TEST_EQUAL( fg->getElementVector()[fg->getElementVector().size() - 1] , 0.0 , 1e-5 , "" );

        std::vector<double> coord(1);
        coord[0] = 0.5;
        COMBIGRID_ERROR_TEST_EQUAL( fg->eval(coord) , 0.125 , 1e-5 , "" );

        delete fg;
        delete domain;
        delete stretchingMaker;
        delete rhs;
      }


      static void test2() {

        std::vector<double> min(2);
        min[0] = 0.4;
        min[1] = 0.4;
        std::vector<double> max(2);
        max[0] = 1.6;
        max[1] = 1.6;
        std::vector<double> sigma(2);
        sigma[0] = 1.0;
        sigma[1] = 1.0;
        std::vector<int> levels(2);
        levels[0] = 5;
        levels[1] = 6;

        AbstractStretchingMaker* stretchingMaker = new TanStretching(0.7); //new UniformStretching(); //new AtanSpecialStretching(); //new UniformStretching();
        GridDomain* domain = new GridDomain( 2 , levels , min , max , (*stretchingMaker) );
        ConstRHS* rhs = new ConstRHS(2.0);

        // make the Poisson computations
        FullGridD* fg = RunPoisson::computeFGPoisson( *domain , levels , sigma , 1.0 , rhs );

        // test the unkowns
        COMBIGRID_ERROR_TEST_EQUAL( fg->getElementVector()[0] , 0.0 , 1e-5 , "" );
        COMBIGRID_ERROR_TEST_EQUAL( fg->getElementVector()[fg->getElementVector().size() - 1] , 0.0 , 1e-5 , "" );

        std::vector<double> coord(2);
        coord[0] = 1.1001;
        coord[1] = 1.1001;
        double val = fg->eval(coord);
        COMBIGRID_ERROR_TEST_EQUAL( val , 0.202233 , 1e-3 , "" );

        delete fg;
        delete domain;
        delete stretchingMaker;
        delete rhs;
      }

      static void test3() {

        std::vector<double> min(3);
        min[0] = 0.0;
        min[1] = 0.0;
        min[2] = 0.0;
        std::vector<double> max(3);
        max[0] = 1.0;
        max[1] = 1.0;
        max[2] = 1.0;
        std::vector<double> sigma(3);
        sigma[0] = 1.0;
        sigma[1] = 1.0;
        sigma[2] = 1.0;
        std::vector<int> levels(3);
        levels[0] = 4;
        levels[1] = 5;
        levels[2] = 3;
        std::vector<int> levelsDom(3);
        levelsDom[0] = 7;
        levelsDom[1] = 7;
        levelsDom[2] = 7;

        AbstractStretchingMaker* stretchingMaker = new TanStretching(0.7); //new UniformStretching(); //new AtanSpecialStretching(); //new UniformStretching();
        GridDomain* domain = new GridDomain( 3 , levelsDom , min , max , (*stretchingMaker) );
        ConstRHS* rhs = new ConstRHS(1.0);

        // make the Poisson computations
        FullGridD* fg = RunPoisson::computeFGPoisson( *domain , levels , sigma , 1.0 , rhs );

        // test the unkowns
        COMBIGRID_ERROR_TEST_EQUAL( fg->getElementVector()[0] , 0.0 , 1e-5 , "" );
        COMBIGRID_ERROR_TEST_EQUAL( fg->getElementVector()[fg->getElementVector().size() - 1] , 0.0 , 1e-5 , "" );

        std::vector<double> coord(3);
        coord[0] = 0.5;
        coord[1] = 0.5;
        coord[2] = 0.5;
        double val = fg->eval(coord);
        COMBIGRID_ERROR_TEST_EQUAL( val , 0.0559004 , 1e-3 , "" );

        delete fg;
        delete domain;
        delete stretchingMaker;
        delete rhs;
      }

      static void test4() {

        std::vector<double> min(4);
        min[0] = 0.0;
        min[1] = 0.0;
        min[2] = 0.0;
        min[3] = 0.0;
        std::vector<double> max(4);
        max[0] = 1.0;
        max[1] = 1.0;
        max[2] = 1.0;
        max[3] = 1.0;
        std::vector<double> sigma(4);
        sigma[0] = 1.0;
        sigma[1] = 1.0;
        sigma[2] = 1.0;
        sigma[3] = 1.0;
        std::vector<int> levels(4);
        levels[0] = 3;
        levels[1] = 3;
        levels[2] = 4;
        levels[3] = 3;
        std::vector<int> levelsDom(4);
        levelsDom[0] = 7;
        levelsDom[1] = 7;
        levelsDom[2] = 7;
        levelsDom[3] = 7;

        AbstractStretchingMaker* stretchingMaker = new TanStretching(0.7); //new UniformStretching(); //new AtanSpecialStretching(); //new UniformStretching();
        GridDomain* domain = new GridDomain( 4 , levelsDom , min , max , (*stretchingMaker) );
        ConstRHS* rhs = new ConstRHS(1.0);

        // make the Poisson computations
        FullGridD* fg = RunPoisson::computeFGPoisson( *domain , levels , sigma , 1.0 , rhs );

        // test the unkowns
        COMBIGRID_ERROR_TEST_EQUAL( fg->getElementVector()[0] , 0.0 , 1e-5 , "" );
        COMBIGRID_ERROR_TEST_EQUAL( fg->getElementVector()[fg->getElementVector().size() - 1] , 0.0 , 1e-5 , "" );

        std::vector<double> coord(4);
        coord[0] = 0.5;
        coord[1] = 0.5;
        coord[2] = 0.5;
        coord[3] = 0.5;
        double val = fg->eval(coord);
        COMBIGRID_ERROR_TEST_EQUAL( val , 0.0454 , 1e-3 , "" );

        delete fg;
        delete domain;
        delete stretchingMaker;
        delete rhs;
      }
  };
}

#endif /* TESTPOISSONPROBLEM_HPP_ */
