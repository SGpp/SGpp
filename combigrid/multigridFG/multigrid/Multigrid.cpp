/* ****************************************************************************
* Copyright (C) 2011 Technische Universitaet Muenchen                         *
* This file is part of the SG++ project. For conditions of distribution and   *
* use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
**************************************************************************** */
// @author Janos Benk (benk@in.tum.de)

#include "Multigrid.hpp"
#include "combigrid/multigridFG/multigrid/ProlongationRestriction.hpp"


combigrid::Multigrid::Multigrid(combigrid::OperatorFG* op , FullGridD* fg , bool createHierarchy) {

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
  int verb = 3;

  const std::vector<bool>& boundaryFlag = fg_tmp->returnBoundaryFlags();
  std::vector<int> levels = fg_tmp->getLevels();
  int dim = fg_tmp->getDimension() , i , maxNrPoint , maxLevel;

  COMBIGRID_OUT_LEVEL2(verb, "Multigrid::Multigrid ... START");

  // create the grids for the
  do {
    fullgrids_.push_back(fg_tmp);
    operators_.push_back(op_tmp);
    unknowns_.push_back( new std::vector<double>( fg_tmp->getNrElements() * op_tmp->getNrSpace() , 0.0) );
    correction_.push_back( new std::vector<double>( fg_tmp->getNrElements() * op_tmp->getNrSpace() , 0.0) );
    rhs_.push_back( new std::vector<double>( fg_tmp->getNrElements() * op_tmp->getNrSpace() , 0.0) );
    op_tmp->getRHS( *(rhs_[depth_]) , i );
    COMBIGRID_OUT_LEVEL4(verb, "Created depth_:" << depth_);

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

    COMBIGRID_OUT_LEVEL4(verb, "refine =" << refine);
    // create a new full grid
    doneRef = false;

    if (refine && createHierarchy) {
      COMBIGRID_OUT_LEVEL4( verb , "Refine Grid ");

      for ( i = 0 ; i < dim ; i++) {
        // decrease the levels for cases
        if ( (levels[i] >= maxLevel) && (fg_tmp->length(i) >= maxNrPoint) ) {
          levels[i] = levels[i] - 1;
          doneRef = true;
          COMBIGRID_OUT_LEVEL4(verb, "REFINE levels[" << i << "]=" << levels[i] << " , maxLevel=" << maxLevel);
        }
      }

      // create a full grid with smaller levels
      fg_tmp = new FullGridD( dim , levels , boundaryFlag );
      fg_tmp->setDomain( domain );
      op_tmp = op_tmp->factory(fg_tmp);
    }
  } while (refine && doneRef && createHierarchy);

  COMBIGRID_OUT_LEVEL2(verb, "Multigrid::Multigrid ... END depth_=" << depth_);
}

combigrid::Multigrid::~Multigrid() {
  for (unsigned int i = 1 ; i < fullgrids_.size() ; i++ ) {
    // delete the full grids which have been created
    delete unknowns_[i];
    delete correction_[i];
    delete rhs_[i];
    delete fullgrids_[i];
    delete operators_[i];
  }

  delete unknowns_[0];
  delete correction_[0];
  delete rhs_[0];
}


void combigrid::Multigrid::solveCS( std::vector<double>& unknowns ,
                                    double errorTol , bool makeFullMG) {

  int  i , vCycle = 0 , verb = 3;
  double  error , errorAct;

  COMBIGRID_OUT_LEVEL2(verb, "Multigrid::solveCS ... START");

  // if it is needed then make a full multigrid cycle
  if (makeFullMG) {
    makeFullMultigrid( unknowns );
  }

  // copy the unknown vector
  for (i = 0 ; i < (int)unknowns.size() ; i++) {
    (*unknowns_[0])[i] = unknowns[i];
  }

  std::vector<double> errVect( unknowns_[0]->size() );
  // measure the error at the beginning
  operators_[0]->multiplyVector( *(unknowns_[0]) , errVect );
  combigrid::vect_diff( &errVect , rhs_[0] );
  errorAct = error = combigrid::l2_norm( &errVect );
  COMBIGRID_OUT_LEVEL2(verb, " Initial error : " << error)

  // do the VCycles as long the residuum is below threshold
  while ( error > errorTol ) {
    // ---- the V-cycle is done with a for loop ----
    // declining branch
    for (i = 0; i < depth_ - 1 ; i++) {
      // pre smooth
      //combigrid::plot_vect(3, verb , unknowns_[i] , "u[i]");
      operators_[i]->doSmoothing( nrGSPre_ , *(unknowns_[i]) , *(rhs_[i]) );
      //combigrid::plot_vect(3, verb , unknowns_[i] , "u[i]");
      COMBIGRID_OUT_LEVEL4(verb, " GS level" << i)
      // restriction
      operators_[i]->multiplyVector( *(unknowns_[i]) , *(correction_[i]) );
      combigrid::vect_diff( correction_[i] , rhs_[i] );
      //combigrid::plot_vect(3, verb , correction_[i] , "Correction[i]");
      combigrid::ProlongationRestriction::restriction( fullgrids_[i], *(correction_[i]), -1.0, fullgrids_[i + 1], *(rhs_[i + 1]), 0.0, operators_[i]->getNrSpace() );
      //combigrid::plot_vect(3, verb , rhs_[i+1] , "rhs[i+1]");
      // it is very important to set the lower level unknowns to "0" this is esential to the correction scheme
      vect_setvalue(unknowns_[i + 1] , 0.0 );
      COMBIGRID_OUT_LEVEL4(verb, " Restriction level" << i)
    }

    // presmooth & postsmooth on the coarsest level
    operators_[depth_ - 1]->doSmoothing( nrGSPre_ + nrGSPost_ , *(unknowns_[depth_ - 1]) , *(rhs_[depth_ - 1]) );
    //combigrid::plot_vect(3, verb , unknowns_[i] , "u[i]");
    COMBIGRID_OUT_LEVEL4(verb, " GS lowest level")

    // rising branch
    for (i = depth_ - 2 ; i >= 0 ; i-- ) {
      // prolongation
      //combigrid::plot_vect(3, verb , unknowns_[i+1] , "u[i+1]");
      //ProlongationRestriction::prolongation( fullgrids_[i], *(unknowns_[i]), 1.0, fullgrids_[i+1], *(unknowns_[i+1]), -1.0, operators_[i]->getNrSpace() );
      combigrid::ProlongationRestriction::prolongation( fullgrids_[i], *(correction_[i]), 0.0, fullgrids_[i + 1], *(unknowns_[i + 1]), 1.0, operators_[i]->getNrSpace() );
      //combigrid::plot_vect(3, verb , correction_[i] , "correction[i]");
      combigrid::vect_add_mul( 1.0 , unknowns_[i] , 1.0 , correction_[i]);
      COMBIGRID_OUT_LEVEL4(verb, " Prolongation level" << i)
      // post smoothing
      operators_[i]->doSmoothing( nrGSPost_ , *(unknowns_[i]) , *(rhs_[i]) );
      //combigrid::plot_vect(3, verb , unknowns_[i] , "u[i]");
      COMBIGRID_OUT_LEVEL4(verb, " GS level" << i)
    }

    // calculate error
    operators_[0]->multiplyVector( (*unknowns_[0]) , errVect);
    combigrid::vect_diff( &errVect , rhs_[0] );
    errorAct = combigrid::l2_norm( &errVect );
    COMBIGRID_OUT_LEVEL2(verb, " ===== END V-cycle ========= error:" << errorAct << " , vCycle:" << vCycle )

    // if the smoother is not working properly than increase the
    if ( errorAct > error ) {
      nrGSPost_ = 2 * nrGSPost_;
      nrGSPre_ = nrGSPre_ + 1;
    }

    if ( (vCycle >= 10) && (vCycle % 10 == 0)) {
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
}

void combigrid::Multigrid::makeFullMultigrid( std::vector<double>& unknowns ) {

  int  i , verb = 3;
  COMBIGRID_OUT_LEVEL2(verb, " combigrid::Multigrid::makeFullMultigrid ... START ");

  for ( int depth = depth_ - 1; depth >= 1 ; depth--) {
    // here make a V cycle up to level "depth"
    for (i = depth; i < depth_ - 1 ; i++) {
      // pre smooth
      operators_[i]->doSmoothing( nrGSPre_ , *(unknowns_[i]) , *(rhs_[i]) );
      COMBIGRID_OUT_LEVEL4(verb, " GS level" << i)
      // restriction
      operators_[i]->multiplyVector( *(unknowns_[i]) , *(correction_[i]) );
      combigrid::vect_diff( correction_[i] , rhs_[i] );
      combigrid::ProlongationRestriction::restriction( fullgrids_[i], *(correction_[i]), -1.0, fullgrids_[i + 1], *(rhs_[i + 1]), 0.0, operators_[i]->getNrSpace() );
      // it is very important to set the lower level unknowns to "0" this is essential to the correction scheme
      vect_setvalue(unknowns_[i + 1] , 0.0 );
      COMBIGRID_OUT_LEVEL4(verb, " Restriction level" << i)
    }

    // presmooth & postsmooth on the coarsest level
    operators_[depth_ - 1]->doSmoothing( nrGSPre_ + nrGSPost_ , *(unknowns_[depth_ - 1]) , *(rhs_[depth_ - 1]) );
    //combigrid::plot_vect(3, verb , unknowns_[i] , "u[i]");
    COMBIGRID_OUT_LEVEL4(verb, " GS lowest level")

    // rising branch
    for (i = depth_ - 2 ; i >= depth ; i-- ) {
      // prolongation
      combigrid::ProlongationRestriction::prolongation( fullgrids_[i], *(correction_[i]), 0.0, fullgrids_[i + 1], *(unknowns_[i + 1]), 1.0, operators_[i]->getNrSpace() );
      combigrid::vect_add_mul( 1.0 , unknowns_[i] , 1.0 , correction_[i]);
      COMBIGRID_OUT_LEVEL4(verb, " Prolongation level" << i)
      // post smoothing
      operators_[i]->doSmoothing( nrGSPost_ , *(unknowns_[i]) , *(rhs_[i]) );
      //combigrid::plot_vect(3, verb , unknowns_[i] , "u[i]");
      COMBIGRID_OUT_LEVEL4(verb, " GS level" << i)
    }

    // project the solution from level depth to depth-1, and make some additional smoothing
    combigrid::ProlongationRestriction::prolongation( fullgrids_[depth - 1], *(unknowns_[depth - 1]), 0.0,
        fullgrids_[depth], *(unknowns_[depth]), 1.0, operators_[depth]->getNrSpace() );
    // additionally make some smoothing
    operators_[depth - 1]->doSmoothing( nrGSPre_ + nrGSPost_ , *(unknowns_[depth - 1]) , *(rhs_[depth - 1]) );
  }

  combigrid::vect_add_mul( 0.0 , &unknowns , 1.0 , unknowns_[0] );

  COMBIGRID_OUT_LEVEL2(verb, " combigrid::Multigrid::makeFullMultigrid ... END ");
}

void combigrid::Multigrid::solveSmoothing( std::vector<double>& unknowns , double errorTol) {
  int  i , vCycle = 0 , verb = 4;
  double  error , errorAct;

  COMBIGRID_OUT_LEVEL2(verb, "Multigrid::solveSmoothing ... START");

  // copy the unknown vector
  for (i = 0 ; i < (int)unknowns.size() ; i++) {
    (*unknowns_[0])[i] = unknowns[i];
  }

  std::vector<double> errVect( unknowns_[0]->size() );
  // measure the error at the beginning
  operators_[0]->multiplyVector( *(unknowns_[0]) , errVect );
  combigrid::vect_diff( &errVect , rhs_[0] );
  errorAct = error = combigrid::l2_norm( &errVect );
  COMBIGRID_OUT_LEVEL2(verb, "Initial error : " << error)

  // do the VCycles as long the residuum is below threshold
  while ( error > errorTol ) {
    // do smoothing
    operators_[0]->doSmoothing( 1 , *(unknowns_[0]) , *(rhs_[0]) );
    //combigrid::plot_vect(3, verb , unknowns_[0] , "u[0]");
    //combigrid::plot_vect(3, verb , rhs_[0] , "rhs[0]");
    // calculate error
    operators_[0]->multiplyVector( (*unknowns_[0]) , errVect);
    combigrid::vect_diff( &errVect , rhs_[0] );
    errorAct = combigrid::l2_norm( &errVect );
    COMBIGRID_OUT_LEVEL2(verb, " error:" << errorAct << " vCycle:" << vCycle );
    error = errorAct;
    vCycle++;
  }

  // copy the solution vector back
  for (i = 0 ; i < (int)unknowns.size() ; i++) {
    unknowns[i] = (*unknowns_[0])[i];
  }

  COMBIGRID_OUT_LEVEL2(verb, "Multigrid::solveSmoothing ... END error:" << errorAct);
}


void combigrid::Multigrid::solveCG( std::vector<double>& unknowns , double errorTol) {
  int  ii , vCycle = 0 , verb = 4;
  double  error , errorAct;
  COMBIGRID_OUT_LEVEL2(verb, "Multigrid::solveCG ... START");

  // copy the unknown vector
  for (ii = 0 ; ii < (int)unknowns.size() ; ii++) {
    (*unknowns_[0])[ii] = unknowns[ii];
  }

  std::vector<double> errVect( unknowns_[0]->size() , 0.0);
  std::vector<double> p( unknowns_[0]->size() , 0.0 );
  std::vector<double> Ap( unknowns_[0]->size() , 0.0 );
  double rsold , rsnew , alpha , tmp;
  // measure the error at the beginning
  //r=b-A*x;
  operators_[0]->multiplyVector( *(unknowns_[0]) , errVect );
  combigrid::vect_diff( &errVect , rhs_[0] );
  combigrid::vect_add_mul( -1.0 , &errVect , 0.0 , &errVect);
  errorAct = error = combigrid::l2_norm( &errVect );
  //p=r;
  p = errVect;
  //rsold=r'*r;
  combigrid::scalar_product( &errVect , &errVect , rsold);
  COMBIGRID_OUT_LEVEL2(verb, "Initial error : " << error)

  // do the VCycles as long the residuum is below threshold
  while ( error > errorTol ) {
    //Ap=A*p;
    operators_[0]->multiplyVector( p , Ap );
    //alpha=rsold/(p'*Ap);
    combigrid::scalar_product( &p , &Ap , tmp );
    alpha = rsold / tmp;
    //x=x+alpha*p;
    combigrid::vect_add_mul( 1.0 , unknowns_[0] , alpha , &p );
    //r=r-alpha*Ap;
    combigrid::vect_add_mul( 1.0 , &errVect , (-1)*alpha , &Ap );
    //rsnew=r'*r;
    combigrid::scalar_product( &errVect , &errVect , rsnew );

    //if sqrt(rsnew)<1e-10 break; end
    if (::sqrt(rsnew) < errorTol) {
      break;
    }

    //p=r+rsnew/rsold*p;
    combigrid::vect_add_mul( rsnew / rsold , &p , 0.0 , &Ap );
    combigrid::vect_add_mul( 1.0 , &p , 1.0 , &errVect );
    //rsold=rsnew;
    rsold = rsnew;
    error = ::sqrt(rsnew);
    COMBIGRID_OUT_LEVEL2(verb, "====== error : " << error << "  it:" << vCycle );
    vCycle++;
  }

  // copy the solution vector back
  for (ii = 0 ; ii < (int)unknowns.size() ; ii++) {
    unknowns[ii] = (*unknowns_[0])[ii];
  }

  COMBIGRID_OUT_LEVEL2(verb, "Multigrid::solveCG ... END error:" << errorAct);

}
