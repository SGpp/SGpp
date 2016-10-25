// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <sgpp_base.hpp>

#include <iostream>

using sgpp::base::Grid;
using sgpp::base::DataVector;
using sgpp::base::SurplusRefinementFunctor;
using sgpp::base::HashRefinementBoundaries;
using sgpp::base::SubspaceRefinement;

int main() {
  const size_t d = 2;
  const size_t l = 3;

  std::unique_ptr<Grid> grid = Grid::createLinearBoundaryGrid(d);
  grid->getGenerator().regular(l);
  const size_t N = grid->getSize();

  std::cout << "grid = " << grid->getStorage().toString() << "\n";
  std::cout << "grid size = " << N << "\n";

  DataVector tmp(N, 0.0);
  const size_t indexToRefine = 5;
  // const size_t indexToRefine = 35;
  tmp[indexToRefine] = 1.0;
  const size_t numberOfPointsToRefine = 1;

  // grid->refine(tmp, refine_points);

  SurplusRefinementFunctor functor(tmp, numberOfPointsToRefine);
  HashRefinementBoundaries hashRefinement;
  SubspaceRefinement subspaceRefinement(&hashRefinement);

  /*AbstractRefinement::refinement_container_type refinablePoints;
   subspaceRefinement.collectRefinablePoints(grid->getStorage(), functor, refinablePoints);
   std::cout << "refinablePoints.size() = " << refinablePoints.size() << "\n";
   std::cout << "refinablePoints[0] seq = " << refinablePoints[0].first->getSeq() << "\n";
   std::cout << "refinablePoints[0] coord = [" <<
   refinablePoints[0].first->getIndex().getCoordsString() << "]\n";
   std::cout << "refinablePoints[0] value = " << refinablePoints[0].second << "\n";*/

  subspaceRefinement.free_refine(grid->getStorage(), functor);
}
