/*
 * RunPoisson.cpp
 *
 *  Created on: May 30, 2011
 *      Author: benk
 */

#include "RunPoisson.hpp"
#include "solver/multigridFG/multigrid/Multigrid.hpp"

using namespace combigrid;


FullGridD* RunPoisson::computeFGPoisson(
		GridDomain& domain,
		const std::vector<int>& levels,
		const std::vector<double>& sigma ,
		double startValue ,
		const CallBackRHS* callbackRHS) {

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
	COMBIGRID_OUT_LEVEL2( verb , "RunPoisson::computeFGPoisson ... create PoissonOperator");
	PoissonOperator op( fg , sigma, callbackRHS );

	// create the multigrid method
	COMBIGRID_OUT_LEVEL2( verb , "RunPoisson::computeFGPoisson ... create Multigrid");
	Multigrid multigrid( &(op) , fg );

	// solve using only the smoother
	COMBIGRID_OUT_LEVEL2( verb , "RunPoisson::computeFGPoisson ... solve multigird");
	unknowns.resize( fg->getNrElements() , startValue );
	multigrid.solveCS( unknowns , 1e-8 );

	// copy the solution back
	COMBIGRID_OUT_LEVEL2( verb , "RunPoisson::computeFGPoisson ... write solution back and return");
	for (int i = 0 ; i < fg->getNrElements() ; i++){
		fg->getElementVector()[i] = unknowns[i];
		//COMBIGRID_OUT_LEVEL2( verb , "unknowns["<<i<<"]="<<unknowns[i]);
	}

	// return full grid
	return fg;
}
