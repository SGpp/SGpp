/*
 * test_main.cpp
 *
 *  Created on: Feb 22, 2011
 *      Author: benk
 */


#include "combigrid/tests/TestFullGrid.hpp"
#include "combigrid/tests/TestCombiGridKernel.hpp"
#include "combigrid/tests/TestSerialCombiGrid.hpp"
#include "combigrid/tests/TestSerialS_CTGrid.hpp"
#include "combigrid/tests/TestSGppConverter.hpp"
#include "combigrid/tests/TestStretching.hpp"
#include "combigrid/tests/TestStretchingCombi.hpp"


int main(int argc, char** argv)
{
	// test the full grid
	combigrid::TestFullGrid::test_all_cases();

	// test the combi grid kernel
	combigrid::TestCombiGridKernel::test_all_cases();

	// test serial combi grid
	combigrid::TestSerialCombiGrid::test_all_cases();

	// test serial combi grid
	combigrid::TestSerialS_CTGrid::test_all_cases();

	// test the converter
	combigrid::TestSGppConverter::test_all_cases();

	// test the stretched full grids
	combigrid::TestStretching::test_all_cases();

	// test the stretched combi grids
	combigrid::TestStretchingCombi::test_all_cases();

	return 0;
}
