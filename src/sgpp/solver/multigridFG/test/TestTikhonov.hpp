/*
 * TestTikhonov.hpp
 *
 *  Created on: Feb 23, 2011
 *      Author: benk
 */

#ifndef TESTTIKHONOVPROBLEM_HPP_
#define TESTTIKHONOVPROBLEM_HPP_

#include "combigrid.hpp"
#include "solver/multigridFG/utils/RunTikhonov.hpp"

namespace combigrid{

/** Class to thes the all the full grids functionalities */
class TestTikhonov{
public:

	TestTikhonov();

    static void test_all_cases() {
    	COMBIGRID_OUT_LEVEL3( 4 , "Tikhonov testing ... ");
    	// 1D tests
    	test1();
    	test12();
    	test13();

    	// 2D tests
    	test2();
    	test201();
    	test21();
    }

    static double testfunction2D( double x1, double x2){
    	return (x1+x2);
    }

    static void test1(){

    	std::vector<double> min(1); min[0] = 0.8;
    	std::vector<double> max(1); max[0] = 1.2;
    	std::vector<int> levels(1); levels[0] = 2;

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

    	delete fg;
    	delete domain;
    	delete stretchingMaker;
    }

    static void test12(){

    	std::vector<double> min(1); min[0] = 0.8;
    	std::vector<double> max(1); max[0] = 1.2;
    	std::vector<int> levels(1); levels[0] = 5;

    	AbstractStretchingMaker* stretchingMaker = new UniformStretching();
    	GridDomain* domain = new GridDomain( 1 , levels , min , max , (*stretchingMaker) );

    	// make the Tikhonov computations
    	FullGridD* fg = RunTikhonov::computeFGTikhonov( *domain , levels , 3e-6 ,
    			"../../../datasets/regression/1DOption/X.txt" , "../../../datasets/regression/1DOption/Y.txt");

    	// test directly the unknowns
    	std::vector<double> coord(1); coord[0] = 0.81;
    	COMBIGRID_ERROR_TEST_EQUAL( fg->eval(coord) , 0.05 , 1e-1 , "" );

    	delete fg;
    	delete domain;
    	delete stretchingMaker;
    }

    static void test13(){

    	std::vector<double> min(1); min[0] = 0.6;
    	std::vector<double> max(1); max[0] = 1.4;
    	std::vector<int> levels(1); levels[0] = 5;

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

    static void test2(){

    	std::vector<double> min(2); min[0] = 0.9; min[1] = 0.9;
    	std::vector<double> max(2); max[0] = 1.1; max[1] = 1.1;
    	std::vector<int> levels(2); levels[0] = 3; levels[1] = 3;

    	AbstractStretchingMaker* stretchingMaker = new UniformStretching();
    	GridDomain* domain = new GridDomain( 2 , levels , min , max , (*stretchingMaker) );

    	// make the Tikhonov computations
    	FullGridD* fg = RunTikhonov::computeFGTikhonov( *domain , levels , 0.001 ,
    			"../../../datasets/regression/2DPlane/X.txt" , "../../../datasets/regression/2DPlane/Y.txt");

    	// test the unkowns
    	COMBIGRID_ERROR_TEST_EQUAL( fg->getElementVector()[0] , 1.8 , 1e-2 , "" );
    	COMBIGRID_ERROR_TEST_EQUAL( fg->getElementVector()[fg->getElementVector().size()-1] , 2.2 , 1e-5 , "" );

    	delete fg;
    	delete domain;
    	delete stretchingMaker;
    }

    static void test201(){

    	std::vector<double> min(2); min[0] = 0.9; min[1] = 0.9;
    	std::vector<double> max(2); max[0] = 1.1; max[1] = 1.1;
    	std::vector<int> levels(2); levels[0] = 4; levels[1] = 3;

    	AbstractStretchingMaker* stretchingMaker = new AtanSpecialStretching();
    	GridDomain* domain = new GridDomain( 2 , levels , min , max , (*stretchingMaker) );

    	// make the Tikhonov computations
    	FullGridD* fg = RunTikhonov::computeFGTikhonov( *domain , levels , 0.001 ,
    			"../../../datasets/regression/2DPlane/X.txt" , "../../../datasets/regression/2DPlane/Y.txt");

    	// test the unkowns
    	COMBIGRID_ERROR_TEST_EQUAL( fg->getElementVector()[0] , 1.8 , 1e-4 , "" );
    	COMBIGRID_ERROR_TEST_EQUAL( fg->getElementVector()[fg->getElementVector().size()-1] , 2.2 , 1e-4 , "" );

    	delete fg;
    	delete domain;
    	delete stretchingMaker;
    }

    static void test21(){

    	std::vector<double> min(2); min[0] = 0.9; min[1] = 0.9;
    	std::vector<double> max(2); max[0] = 1.1; max[1] = 1.1;
    	std::vector<int> levels(2); levels[0] = 4; levels[1] = 3;

    	AbstractStretchingMaker* stretchingMaker = new AtanSpecialStretching();
    	GridDomain* domain = new GridDomain( 2 , levels , min , max , (*stretchingMaker) );

    	// make the Tikhonov computations
    	FullGridD* fg = RunTikhonov::computeFGTikhonov( *domain , levels , 2e-5 ,
    			"../../../datasets/regression/2DOption/X.txt" , "../../../datasets/regression/2DOption/Y.txt");

    	std::vector<double> coord(2); coord[0] = 0.81; coord[1] = 0.81;
    	COMBIGRID_ERROR_TEST_EQUAL( fg->eval(coord) , 0.05 , 1e-1 , "" );

    	delete fg;
    	delete domain;
    	delete stretchingMaker;
    }

};
}

#endif /* TESTTIKHONOVPROBLEM_HPP_ */
