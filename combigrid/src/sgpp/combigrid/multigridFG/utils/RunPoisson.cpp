// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <sgpp/combigrid/multigridFG/utils/RunPoisson.hpp>
#include <sgpp/combigrid/multigridFG/multigrid/Multigrid.hpp>
#include <sgpp/combigrid/multigridFG/multigrid/MultigridFAS.hpp>

using namespace combigrid;

FullGridD* combigrid::RunPoisson::computeFGPoisson(
  combigrid::GridDomain& domain,
  const std::vector<int>& levels,
  const std::vector<double>& sigma ,
  double startValue ,
  const combigrid::CallBackRHS* callbackRHS) {

  int dimensions = (int)sigma.size();
  int verb = 6;
  std::vector<double> unknowns;
  FullGridD* fg;
  // first read in the results

  // create the full grid
  COMBIGRID_OUT_LEVEL2( verb , "RunPoisson::computeFGPoisson ... create Full Grid");
  fg = new FullGridD( dimensions , levels );
  fg->createFullGrid();
  fg->setDomain( &domain );
  //fg->getElementVector()

  // create the Tikhonov operator
  COMBIGRID_OUT_LEVEL2( verb , "RunPoisson::computeFGPoisson ... create combigrid::PoissonOperator");
  combigrid::PoissonOperator op( fg , sigma, callbackRHS );

  // create the multigrid method
  COMBIGRID_OUT_LEVEL2( verb , "RunPoisson::computeFGPoisson ... create combigrid::Multigrid");
  combigrid::Multigrid multigrid( &(op) , fg );

  // solve using only the smoother
  COMBIGRID_OUT_LEVEL2( verb , "RunPoisson::computeFGPoisson ... solve multigird");
  unknowns.resize( fg->getNrElements() , startValue );
  multigrid.solveCS( unknowns , 1e-8 , true );
  //multigrid.solveSmoothing( unknowns , 1e-8);
  //multigrid.solveCG( unknowns , 1e-6 );
  //multigridFAS.solveFAS( unknowns , 1e-8 );

  // copy the solution back
  COMBIGRID_OUT_LEVEL2( verb , "RunPoisson::computeFGPoisson ... write solution back and return");

  for (int i = 0 ; i < fg->getNrElements() ; i++) {
    fg->getElementVector()[i] = unknowns[i];
    //COMBIGRID_OUT_LEVEL2( verb , "unknowns["<<i<<"]="<<unknowns[i]);
  }

  // return full grid
  return fg;
}

void combigrid::RunPoisson::computeFGPoisson_FG(
  FullGridD* fg ,
  const std::vector<double>& sigma ,
  std::vector<double>& unknowns ,
  const CallBackRHS* callbackRHS ) {

  int verb = 6;

  // create the Tikhonov operator
  COMBIGRID_OUT_LEVEL2( verb , "RunPoisson::computeFGPoisson ... create combigrid::PoissonOperator");
  combigrid::PoissonOperator op( fg , sigma, callbackRHS );

  // create the multigrid method
  COMBIGRID_OUT_LEVEL2( verb , "RunPoisson::computeFGPoisson ... create combigrid::Multigrid");
  combigrid::Multigrid multigrid( &(op) , fg );

  // solve using only the smoother
  COMBIGRID_OUT_LEVEL2( verb , "RunPoisson::computeFGPoisson ... solve multigird");
  multigrid.solveCS( unknowns , 1e-8 , true );
  //multigrid.solveSmoothing( unknowns , 1e-8);
  //multigrid.solveCG( unknowns , 1e-6 );
  //multigridFAS.solveFAS( unknowns , 1e-8 );

  // copy the solution back
  COMBIGRID_OUT_LEVEL2( verb , "RunPoisson::computeFGPoisson ... write solution back and return");

  for (int i = 0 ; i < fg->getNrElements() ; i++) {
    fg->getElementVector()[i] = unknowns[i];
    //COMBIGRID_OUT_LEVEL2( verb , "unknowns["<<i<<"]="<<unknowns[i]);
  }
}