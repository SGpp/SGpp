/* ****************************************************************************
* Copyright (C) 2011 Technische Universitaet Muenchen                         *
* This file is part of the SG++ project. For conditions of distribution and   *
* use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
**************************************************************************** */
// @author Janos Benk (benk@in.tum.de)


#include <sgpp/combigrid/tests/TestFullGrid.hpp>
#include <sgpp/combigrid/tests/TestCombiGridKernel.hpp>
#include <sgpp/combigrid/tests/TestSerialCombiGrid.hpp>
#include <sgpp/combigrid/tests/TestSerialS_CTGrid.hpp>
#include <sgpp/combigrid/tests/TestSGppConverter.hpp>
#include <sgpp/combigrid/tests/TestStretching.hpp>
#include <sgpp/combigrid/tests/TestStretchingCombi.hpp>


int main(int argc, char** argv) {
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
