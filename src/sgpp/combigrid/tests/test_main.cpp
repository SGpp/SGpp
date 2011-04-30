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

using namespace combigrid;

int main(int argc, char** argv)
{
	// test the full grid
	TestFullGrid::test_all_cases();

	// test the combi grid kernel
	TestCombiGridKernel::test_all_cases();

	// test serial combi grid
	TestSerialCombiGrid::test_all_cases();

	// test serial combi grid
	TestSerialS_CTGrid::test_all_cases();

	// test the converter
	TestSGppConverter::test_all_cases();

	// test the stretched full grids
	TestStretching::test_all_cases();

	// test the stretched combi grids
	TestStretchingCombi::test_all_cases();

	return 0;
}
