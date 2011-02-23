/*
 * test_main.cpp
 *
 *  Created on: Feb 22, 2011
 *      Author: benk
 */

#include "combigrid/utils/combigrid_ultils.hpp"
#include "combigrid/tests/TestFullGrid.hpp"
#include "combigrid/tests/TestCombiGridKernel.hpp"

using namespace combigrid;

int main(int argc, char** argv)
{
	// test the full grid
	TestFullGrid::test_all_cases();

	// test the combi grid kernel
	TestCombiGridKernel::test_all_cases();

	return 0;
}
