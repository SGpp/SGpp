/*
 * test_main.cpp
 *
 *  Created on: Feb 22, 2011
 *      Author: benk
 */


#include "solver/multigridFG/test/TestTikhonov.hpp"

using namespace combigrid;

int main(int argc, char** argv)
{
	// test the Tikhonov solver
	TestTikhonov::test_all_cases();

	return 0;
}
