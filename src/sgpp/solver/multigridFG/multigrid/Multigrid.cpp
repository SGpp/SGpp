/*
 * Multigrid.cpp
 *
 *  Created on: May 16, 2011
 *      Author: benk
 */

#include "Multigrid.hpp"
#include "solver/multigridFG/multigrid/ProlongationRestriction.hpp"

using namespace combigrid;

Multigrid::Multigrid(OperatorFG* op, FullGridD* fg) {

	// - create the hierarchy of grids
	//     - in each direction there should be at least 3 points
	//     - if in each direction the number of points is <= 3 then stop with the hierarchy

	FullGridD* fg_tmp = fg;
	OperatorFG* op_tmp = op;
	GridDomain* domain = fg->getDomain();
	depth_ = 0;
	nrGSPre_ = 1;
	nrGSPost_ = 1;
	bool refine = false;

	const std::vector<bool>& boundaryFlag = fg_tmp->returnBoundaryFlags();
	std::vector<int> levels = fg_tmp->getLevels();
	int dim = fg_tmp->getDimension() , i , maxNrPoint , maxLevel;

	// create the grids for the
	do {
		fullgrids_.push_back(fg_tmp);
		operators_.push_back(op_tmp);
		unknowns_.push_back( new std::vector<double>( fg_tmp->getNrElements() * op_tmp->getNrSpace() , 0.0) );
		correction_.push_back( new std::vector<double>( fg_tmp->getNrElements() * op_tmp->getNrSpace() , 0.0) );
		rhs_.push_back( new std::vector<double>( fg_tmp->getNrElements() * op_tmp->getNrSpace() , 0.0) );
		op_tmp->getRHS( *(rhs_[depth_]) , i );

		depth_++;
		// look for the levels
		refine = false;
		maxNrPoint = fg_tmp->length( 0 );
		maxLevel = (boundaryFlag[0]) ? levels[0]-1 : levels[0] ;
		for ( i = 1 ; i < dim ; i++){
			maxNrPoint = ( fg_tmp->length(i) > maxNrPoint) ? fg_tmp->length(i) : maxNrPoint ;
			maxLevel = (boundaryFlag[i]) ? levels[i]-1 : levels[i] ;
		}
		// only if we have more than
		if (maxNrPoint > 3) { refine = true; }
		// create a new full grid
        if (refine){
        	for ( i = 0 ; i < dim ; i++){
        		// decrease the levels for cases
        		if ( (levels[i] >= maxLevel) || (fg_tmp->length(i) > maxNrPoint) ) {
        			levels[i] = levels[i] - 1;
        		}
        	}
        	// create a full grid with smaller levels
        	fg_tmp = new FullGridD( dim , levels , boundaryFlag );
        	fg_tmp->setDomain( domain );
        	op_tmp = op_tmp->factory(fg_tmp);
        }
	} while (refine);
}

Multigrid::~Multigrid() {
	for (unsigned int i = 1 ; i < fullgrids_.size() ; i++ ){
		// delete the full grids which have been created
		delete fullgrids_[i];
	}
}


void Multigrid::solveCS( std::vector<double>& unknowns , double errorTol){

	int  i , vCycle = 0 , verb = 4;
	double  error , errorAct;

	COMBIGRID_OUT_LEVEL2(verb,"Multigrid::solveCS ... START");

	// copy the unknown vector
	for (i = 0 ; i < (int)unknowns.size() ; i++) { (*unknowns_[0])[i] = unknowns[i]; }

	std::vector<double> errVect( unknowns_[0]->size() );
	// measure the error at the beginning
	operators_[0]->multiplyVector( *(unknowns_[0]) , errVect );
	combigrid::vect_diff( &errVect , rhs_[0] );
	errorAct = error = combigrid::l2_norm( &errVect );
	COMBIGRID_OUT_LEVEL2(verb," error : " << error)

	// do the VCycles as long the residuum is below threshold
	while ( error > errorTol ){
	    // ---- the V-cycle is done with a for loop ----
		// declining branch
		for (i = 0; i < depth_-1 ; i++){
			// pre smooth
			operators_[i]->doSmoothing( nrGSPre_ , *(unknowns_[i]) , *(rhs_[i]) );
			COMBIGRID_OUT_LEVEL2(verb," GS level" << i)
			// restriction
			operators_[i]->multiplyVector( *(unknowns_[i]) , *(correction_[i]) );
			combigrid::vect_diff( correction_[i] , rhs_[i] );
			ProlongationRestriction::restriction( fullgrids_[i] , *(correction_[i]) , 1.0, fullgrids_[i+1] , *(rhs_[i+1]) , 0.0, operators_[i]->getNrSpace() );
			COMBIGRID_OUT_LEVEL2(verb," Restriction level" << i)
		}
		// presmooth & postsmooth on the coarsest level
		operators_[depth_]->doSmoothing( nrGSPre_ + nrGSPost_ , *(unknowns_[depth_]) , *(rhs_[depth_]) );
		COMBIGRID_OUT_LEVEL2(verb," GS lowest level")
		// rising branch
		for (i = depth_-2 ; i >= 0 ; i-- ){
			// prolongation
			ProlongationRestriction::prolongation( fullgrids_[i] , *(unknowns_[i]) , 1.0, fullgrids_[i+1] , *(unknowns_[i+1]) , 1.0, operators_[i]->getNrSpace() );
			COMBIGRID_OUT_LEVEL2(verb," Prolongation level" << i)
			// post smoothing
			operators_[depth_]->doSmoothing( nrGSPost_ , *(unknowns_[i]) , *(rhs_[i]) );
			COMBIGRID_OUT_LEVEL2(verb," GS level" << i)
		}
		// calculate error
		operators_[0]->multiplyVector( (*unknowns_[0]) , errVect);
		combigrid::vect_diff( &errVect , rhs_[0] );
		errorAct = combigrid::l2_norm( &errVect );
		// if the smoother is not working properly than increase the
		if ( errorAct > error ){
			nrGSPost_ = nrGSPost_ + 5;
			nrGSPre_++;
		}
		if (vCycle > 10) {
			nrGSPost_ = nrGSPost_ + 2;
			nrGSPre_ ++;
		}
		// actual error update, and increase the number of iterations
		error = errorAct;
		vCycle++;
	}
}

void Multigrid::solveFAS( std::vector<double>& unknowns , double errorTol){
  // todo: implement this

	// the V-cycle should be done with a for loop
}
