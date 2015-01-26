/* ****************************************************************************
* Copyright (C) 2011 Technische Universitaet Muenchen                         *
* This file is part of the SG++ project. For conditions of distribution and   *
* use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
**************************************************************************** */
// @author Janos Benk (benk@in.tum.de)

#include <sgpp/combigrid/multigridFG/test/TestTikhonov.hpp>
#include <sgpp/combigrid/multigridFG/test/TestPoisson.hpp>


int main(int argc, char** argv) {

  // test the Poisson problem solver
  combigrid::TestPoisson::test_all_cases();

  // test the Tikhonov solver
  combigrid::TestTikhonov::test_all_cases();

  return 0;
}
