/* ****************************************************************************
* Copyright (C) 2011 Technische Universitaet Muenchen                         *
* This file is part of the SG++ project. For conditions of distribution and   *
* use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
**************************************************************************** */
// @author Janos Benk (benk@in.tum.de)

#include "MultigridFAS.hpp"
#include <sgpp/combigrid/multigridFG/multigrid/ProlongationRestriction.hpp>


combigrid::MultigridFAS::MultigridFAS(combigrid::OperatorFG* op , FullGridD* fg , bool createHierarchy) {

  // - create the hierarchy of grids
  //     - in each direction there should be at least 3 points
  //     - if in each direction the number of points is <= 3 then stop with the hierarchy

  FullGridD* fg_tmp = fg;
  combigrid::OperatorFG* op_tmp = op;
  combigrid::GridDomain* domain = fg->getDomain();
  depth_ = 0;
  nrGSPre_ = 2;
  nrGSPost_ = 2;
  bool refine = false , doneRef = true;
  int verb = 6;

  const std::vector<bool>& boundaryFlag = fg_tmp->returnBoundaryFlags();
  std::vector<int> levels = fg_tmp->getLevels();
  int dim = fg_tmp->getDimension() , i , maxNrPoint , maxLevel;

  COMBIGRID_OUT_LEVEL2(verb, "MultigridFAS::MultigridFAS ... START");

  // create the grids for the
  do {
    fullgrids_.push_back(fg_tmp);
    operators_.push_back(op_tmp);
    unknowns_.push_back( new std::vector<double>( fg_tmp->getNrElements() * op_tmp->getNrSpace() , 0.0) );
    correction_.push_back( new std::vector<double>( fg_tmp->getNrElements() * op_tmp->getNrSpace() , 0.0) );
    rhs_.push_back( new std::vector<double>( fg_tmp->getNrElements() * op_tmp->getNrSpace() , 0.0) );
    rhs_tmp_.push_back( new std::vector<double>( fg_tmp->getNrElements() * op_tmp->getNrSpace() , 0.0) );
    u_hH_.push_back( new std::vector<double>( fg_tmp->getNrElements() * op_tmp->getNrSpace() , 0.0) );
    lh_.push_back( new std::vector<double>( fg_tmp->getNrElements() * op_tmp->getNrSpace() , 0.0) );
    op_tmp->getRHS( *(rhs_[depth_]) , i );
    combigrid::vect_add_mul( 0.0 , rhs_tmp_[depth_] , 1.0 , rhs_[depth_]);
    COMBIGRID_OUT_LEVEL2(verb, "Created depth_:" << depth_);

    depth_++;
    // look for the levels
    refine = false;
    maxNrPoint = fg_tmp->length( 0 );
    maxLevel = (boundaryFlag[0]) ? levels[0] - 1 : levels[0] ;

    for ( i = 1 ; i < dim ; i++) {
      maxNrPoint = ( fg_tmp->length(i) > maxNrPoint) ? fg_tmp->length(i) : maxNrPoint ;
      maxLevel = (boundaryFlag[i]) ? levels[i] - 1 : levels[i] ;
    }

    // only if we have more than
    if (maxNrPoint > 3) {
      refine = true;
    }

    COMBIGRID_OUT_LEVEL2(verb, "refine =" << refine);
    // create a new full grid
    doneRef = false;

    if (refine && createHierarchy) {
      COMBIGRID_OUT_LEVEL2( verb , "Refine Grid ");

      for ( i = 0 ; i < dim ; i++) {
        // decrease the levels for cases
        if ( (levels[i] >= maxLevel) && (fg_tmp->length(i) >= maxNrPoint) ) {
          levels[i] = levels[i] - 1;
          doneRef = true;
          COMBIGRID_OUT_LEVEL2(verb, "REFINE levels[" << i << "]=" << levels[i] << " , maxLevel=" << maxLevel);
        }
      }

      // create a full grid with smaller levels
      fg_tmp = new FullGridD( dim , levels , boundaryFlag );
      fg_tmp->setDomain( domain );
      op_tmp = op_tmp->factory(fg_tmp);
    }
  } while (refine && doneRef && createHierarchy);

  COMBIGRID_OUT_LEVEL2(verb, "MultigridFAS::MultigridFAS ... END depth_=" << depth_);
}

combigrid::MultigridFAS::~MultigridFAS() {
  for (unsigned int i = 1 ; i < fullgrids_.size() ; i++ ) {
    // delete the full grids which have been created
    delete unknowns_[i];
    delete correction_[i];
    delete rhs_[i];
    delete rhs_tmp_[i];
    delete u_hH_[i];
    delete lh_[i];
    delete fullgrids_[i];
    delete operators_[i];
  }

  delete unknowns_[0];
  delete correction_[0];
  delete rhs_[0];
  delete rhs_tmp_[0];
  delete u_hH_[0];
  delete lh_[0];
}


void combigrid::MultigridFAS::solveFAS( std::vector<double>& unknowns , double errorTol) {

  // todo: get good starting value

  int  i , vCycle = 0 , verb = 4;
  double  error , errorAct;

  COMBIGRID_OUT_LEVEL2(verb, "MultigridFAS::solveFAS ... START");

  // copy the unknown vector
  for (i = 0 ; i < (int)unknowns.size() ; i++) {
    (*unknowns_[0])[i] = unknowns[i];
  }

  std::vector<double> errVect( unknowns_[0]->size() );
  // measure the error at the beginning
  operators_[0]->multiplyVector( *(unknowns_[0]) , errVect );
  combigrid::vect_diff( &errVect , rhs_tmp_[0] );
  errorAct = error = combigrid::l2_norm( &errVect );
  COMBIGRID_OUT_LEVEL2(verb, " error : " << error)

  // do the VCycles as long the residuum is below threshold
  while ( error > errorTol ) {
    // ---- the V-cycle is done with a for loop ----
    // declining branch
    for (i = 0; i < depth_ - 1 ; i++) {
      // pre smooth
      operators_[i]->doSmoothing( nrGSPre_ , *(unknowns_[i]) , *(rhs_tmp_[i]) );
      COMBIGRID_OUT_LEVEL2(verb, " GS level" << i)
      // restriction to the right hand side
      operators_[i]->multiplyVector( *(unknowns_[i]) , *(correction_[i]) );
      combigrid::ProlongationRestriction::restriction( fullgrids_[i], *(correction_[i]), -1.0, fullgrids_[i + 1], *(correction_[i + 1]), 0.0, operators_[i]->getNrSpace() );
      combigrid::vect_add_mul( 0.0 , rhs_tmp_[i + 1] , 1.0 , rhs_[i + 1]);
      combigrid::vect_add_mul( 1.0 , rhs_tmp_[i + 1] , -1.0 , correction_[i + 1]);
      combigrid::ProlongationRestriction::restriction( fullgrids_[i], *(unknowns_[i]), 1.0, fullgrids_[i + 1], *(u_hH_[i + 1]), 0.0, operators_[i]->getNrSpace() );
      combigrid::vect_add_mul( 0.0 , correction_[i + 1] , 1.0 , u_hH_[i + 1]);
      operators_[i + 1]->multiplyVector( *(correction_[i + 1]) , *(lh_[i + 1]) );
      combigrid::vect_add_mul( 1.0 , rhs_tmp_[i + 1] , 1.0 , lh_[i + 1] );
      COMBIGRID_OUT_LEVEL2(verb, " Restriction level" << i)
    }

    // presmooth & postsmooth on the coarsest level
    operators_[depth_ - 1]->doSmoothing( nrGSPre_ + nrGSPost_ , *(unknowns_[depth_ - 1]) , *(rhs_tmp_[depth_ - 1]) );
    COMBIGRID_OUT_LEVEL2(verb, " GS lowest level")

    // rising branch
    for (i = depth_ - 2 ; i >= 0 ; i-- ) {
      // prolongation
      combigrid::vect_add_mul( -1.0 , u_hH_[i + 1] , 1.0 , unknowns_[i + 1] );
      combigrid::ProlongationRestriction::prolongation( fullgrids_[i], *(u_hH_[i]), 0.0, fullgrids_[i + 1], *(u_hH_[i + 1]), 1.0, operators_[i]->getNrSpace() );
      combigrid::vect_add_mul( 1.0 , unknowns_[i] , 1.0 , u_hH_[i]);
      COMBIGRID_OUT_LEVEL2(verb, " Prolongation level" << i)
      // post smoothing
      operators_[i]->doSmoothing( nrGSPost_ , *(unknowns_[i]) , *(rhs_tmp_[i]) );
      COMBIGRID_OUT_LEVEL2(verb, " GS level" << i)
    }

    // calculate error
    operators_[0]->multiplyVector( (*unknowns_[0]) , errVect);
    combigrid::vect_diff( &errVect , rhs_tmp_[0] );
    errorAct = combigrid::l2_norm( &errVect );
    COMBIGRID_OUT_LEVEL2(verb, " ===== END V-cycle ========= error:" << errorAct << " , vCycle:" << vCycle )

    // if the smoother is not working properly than increase the
    if ( errorAct > error ) {
      nrGSPost_ = 2 * nrGSPost_;
      nrGSPre_ = nrGSPre_ + 1;
    }

    if (vCycle > 10) {
      nrGSPost_ = nrGSPost_ + 2;
      nrGSPre_ = nrGSPre_ + 1;
    }

    // actual error update, and increase the number of iterations
    error = errorAct;
    vCycle++;
  }

  // copy the solution vector back
  for (i = 0 ; i < (int)unknowns.size() ; i++) {
    unknowns[i] = (*unknowns_[0])[i];
  }

  COMBIGRID_OUT_LEVEL2(verb, "MultigridFAS::solveFAS ... END");
}
