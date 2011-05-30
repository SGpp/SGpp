/*
 * TestCombiGridKernel.hpp
 *
 *  Created on: Feb 23, 2011
 *      Author: benk
 */

#ifndef TESTCOMBIGRIDKERNEL_HPP_
#define TESTCOMBIGRIDKERNEL_HPP_

#include "combigrid.hpp"
#include "solver/multigridFG/utils/RunTikhonov.hpp"

namespace combigrid{

/** Class to thes the all the full grids functionalities */
class TestTikhonov{
public:

	TestTikhonov();

    static void test_all_cases() {
    	COMBIGRID_OUT_LEVEL3( 4 , "Tikhonov testing ... ");
    	// 1D test
    	test1();

    	// 1D test
    	test2();
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

    	// test the combi grid state after the duplicated spaces were deleted
    	COMBIGRID_ERROR_TEST_EQUAL( fg->getElementVector()[0] , 0.8 , 1e-7 , "" );
    	COMBIGRID_ERROR_TEST_EQUAL( fg->getElementVector()[1] , 0.9 , 1e-7 , "" );
    	COMBIGRID_ERROR_TEST_EQUAL( fg->getElementVector()[2] , 1.0 , 1e-7 , "" );
    	COMBIGRID_ERROR_TEST_EQUAL( fg->getElementVector()[3] , 1.1 , 1e-7 , "" );
    	COMBIGRID_ERROR_TEST_EQUAL( fg->getElementVector()[4] , 1.2 , 1e-7 , "" );

    	delete fg;
    }

    static void test2(){
    	// test the combi grid state after the duplicated spaces were deleted
    	COMBIGRID_ERROR_TEST_EQUAL( 1.0 , 1.0 , 1e-10 , "" );
    }

};
}

#endif /* TESTCOMBIGRIDKERNEL_HPP_ */
