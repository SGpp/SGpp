/*
 * test_main.cpp
 *
 *  Created on: Feb 22, 2011
 *      Author: benk
 */


#include "solver/multigridFG/test/TestTikhonov.hpp"
#include "solver/multigridFG/test/TestPoisson.hpp"


int main(int argc, char** argv)
{

	// test the Poisson problem solver
	combigrid::TestPoisson::test_all_cases();

	// test the Tikhonov solver
	combigrid::TestTikhonov::test_all_cases();

	return 0;
}
