// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include "ProlongationRestriction.hpp"


void combigrid::ProlongationRestriction::prolongation(
  const FullGridD* fgFine ,
  std::vector<double>& vectFine ,
  double coefFine ,
  const FullGridD* fgCoarse ,
  const std::vector<double>& vectCoarse ,
  double coefCoarse ,
  int nrSpace) {

  // consistency check
  COMBIGRID_ERROR_TEST((int)vectFine.size() == (int)nrSpace * fgFine->getNrElements(), "ProlongationRestriction::prolongation "
                       << " SIZE MUST MATCH vectFine.size():" << vectFine.size() <<
                       " nrSpace:" << nrSpace << " fgFine->getNrElements():" << fgFine->getNrElements() );
  COMBIGRID_ERROR_TEST((int)vectCoarse.size() == (int)nrSpace * fgCoarse->getNrElements(), "ProlongationRestriction::prolongation "
                       << " SIZE MUST MATCH vectCoarse.size():" << vectCoarse.size() <<
                       " nrSpace:" << nrSpace << " fgCoarse->getNrElements():" << fgCoarse->getNrElements() );

  const int dim = fgFine->getDimension();
  const int nrFinePoints = fgFine->getNrElements();
  const int nrCellPoints = combigrid::powerOfTwo[dim];
  bool hasDomain = true;
  int tmp_I , linearIndexC , tmp_I2 , verb = 0 , i , axis , d;
  double tmp_D;
  int tOffs = 0;
  std::vector<int> axisIndexC(dim);
  std::vector<int> axisIndexF(dim);
  std::vector<int> mulFactors(dim, 1);
  std::vector<int> mulScalingC(dim, 1);
  std::vector<int> mulScalingF(dim, 1);
  std::vector<double> axisDiv(dim);
  std::vector<double> val(nrSpace, 0.0);
  const std::vector<bool>& boundaryFlag = fgFine->returnBoundaryFlags();
  // both grid have the same domain
  const combigrid::GridDomain* domain = fgFine->getDomain();

  COMBIGRID_OUT_LEVEL2( verb , "ProlongationRestriction::prolongation ... START");

  // first see which axis index (on the fine mesh) have to be multiplied by two
  for (int i = 0 ; i < dim ; i++) {
    mulFactors[i] = (fgFine->length(i) + 1) / fgCoarse->length(i);
  }

  if (domain == 0) {
    // we are on the unit square
    hasDomain = false;
  } else {
    // test if the axis are scaled, if not then do just as on unit square
    // we take in consideration only the first axis (they should be homogeneous )
    hasDomain = domain->get1DDomain(0).isAxisScaled();
  }

  // here we iterate on all fine grid points
  for (int pInd = 0 ; pInd < nrFinePoints ; pInd++) {
    tmp_I = pInd;
    //       - locate the fine point on the coarse mesh (fine the N-dim cell)
    //       - make interpolation (eval of the basis function) on the given coord (caz stetching this is not in the middle)
    //            - take the min per axis
    //       - set the value of the fine mesh point

    // calc fine grid point index
    tmp_I = pInd;
    linearIndexC = 0;

    for (i = dim - 1 ; i >= 0 ; i--) {
      axisIndexF[i] = tmp_I / (fgFine->getOffset(i));
      tmp_I = tmp_I % fgFine->getOffset(i);

      // locate the lower corner point of the coarse grid
      if (fgCoarse->getLevels()[i] == fgFine->getLevels()[i]) {
        // NO LEVEL reduction
        axisIndexC[i] = (axisIndexF[i] >= (fgCoarse->length(i) - 1)) ? (fgCoarse->length(i) - 2) : (axisIndexF[i]) ;
      } else {
        // LEVEL reduction , we subtract 1 because of the boundary
        if (boundaryFlag[i]) {
          axisIndexC[i] = (axisIndexF[i] / 2 >= (fgCoarse->length(i) - 1)) ? (fgCoarse->length(i) - 2) : (axisIndexF[i] / 2) ;
        } else {
          axisIndexC[i] = (axisIndexF[i] / 2 >= (fgCoarse->length(i) - 1)) ? (fgCoarse->length(i) - 2) : ( COMBIGRID_IMAX(axisIndexF[i] - 1, 0) / 2);
        }
      }

      linearIndexC = linearIndexC + axisIndexC[i] * fgCoarse->getOffset(i);
      //COMBIGRID_OUT_LEVEL2( verb , "makeLinearProlongation LOC i:" << i << ", axisIndexF[i]:" << axisIndexF[i]
      //                << ", axisIndexC[i]:" << axisIndexC[i] << " Lc:" << fgCoarse->level(i) << " Lf:" << gFine->level(i));
      //COMBIGRID_OUT_LEVEL2( verb , "makeLinearProlongation LOC i:" << i << " sizeC:" << fgCoarse->axisScaling(i).size() <<
      //                       " sizeF:" << fgFine->axisScaling(i).size() );
    }

    if (hasDomain) {
      // get the intersection per axis
      for (i = 0 ; i < dim ; i++) {
        // get the scaling vector
        const std::vector<double>& scalingAxis = domain->get1DDomain(i).axisScaling();
        // multiplicator for the scaling, the domain has potentially higher resolution
        mulScalingF[i] = combigrid::powerOfTwo[domain->get1DDomain(i).getLevel() - fgFine->getLevels()[i]];
        mulScalingC[i] = combigrid::powerOfTwo[domain->get1DDomain(i).getLevel() - fgCoarse->getLevels()[i]];

        if (boundaryFlag[i]) {
          axisDiv[i] = scalingAxis[mulScalingF[i] * axisIndexF[i]] - scalingAxis[mulScalingC[i] * axisIndexC[i]];
          //COMBIGRID_OUT_LEVEL2( verb , "makeLinearProlongation INT i:" << i << " F:" << gFine->axisScaling(i)[axisIndexF[i]]
          //          << " C1:" << gCoarse->axisScaling(i)[axisIndexC[i]] << " C2:" << gCoarse->axisScaling(i)[axisIndexC[i]+1]);
          axisDiv[i] = axisDiv[i] / (scalingAxis[mulScalingC[i] * (axisIndexC[i] + 1)] - scalingAxis[mulScalingC[i] * axisIndexC[i]]);
        } else {
          // the scaling starts at position 1 not at position 0
          axisDiv[i] = scalingAxis[mulScalingF[i] * (1 + axisIndexF[i])] - scalingAxis[mulScalingC[i] * (1 + axisIndexC[i])];
          //COMBIGRID_OUT_LEVEL2( verb , "makeLinearProlongation INT i:" << i << " F:" << gFine->axisScaling(i)[1+axisIndexF[i]]
          //            << " C1:" << gCoarse->axisScaling(i)[1+axisIndexC[i]] << " C2:" << gCoarse->axisScaling(i)[1+axisIndexC[i]+1]);
          axisDiv[i] = axisDiv[i] / (scalingAxis[mulScalingC[i] * (2 + axisIndexC[i])] - scalingAxis[mulScalingC[i] * (1 + axisIndexC[i])]);
          // since the finer grid boundary points can be outside the cell
          axisDiv[i] = COMBIGRID_DMAX( COMBIGRID_DMIN(axisDiv[i], 1.0) , 0.0);
        }
      }
    } else {
      // no stretching
      for (i = 0 ; i < dim ; i++) {
        // no stretching, so just 0.0, or 0.5
        if ( mulFactors[i] * axisIndexC[i] == axisIndexF[i] ) {
          axisDiv[i] = 0.0;
        } else {
          axisDiv[i] = 0.5;
        }

        // in special cases we need different values than 0.0 or 0.5 , e.g 1.0 or -0.5, 1.0 or 1.5
        if ( boundaryFlag[i] && ( axisIndexF[i] == fgFine->length(i) - 1) ) {
          axisDiv[i] = 1.0;
        }

        if ( (!boundaryFlag[i]) && ( axisIndexF[i] == fgFine->length(i) - 2) ) {
          axisDiv[i] = 1.0;
        }

        if ( (!boundaryFlag[i]) && ( axisIndexF[i] == fgFine->length(i) - 1) ) {
          axisDiv[i] = 1.5;
        }

        if ( (!boundaryFlag[i]) && ( axisIndexF[i] <= 0 )) {
          axisDiv[i] = -0.5;
        }

        // limit the value
        axisDiv[i] = COMBIGRID_DMAX( COMBIGRID_DMIN(axisDiv[i], 1.0) , 0.0);
      }
    }

    // loop over 2^D points and calculate the value of the interpolation
    for (i = 0 ; i < nrSpace ; i++ ) val[i] = 0.0 ;

    tOffs = 0;

    for (i = 0 ; i < nrCellPoints ; i++) {
      tmp_D = 1.0;
      tOffs = 0;

      // calculate the basis function value in N-dim
      for (axis = 0 ; axis < dim ; axis++) {
        // calculate the index on this axis
        tmp_I2 = (i / combigrid::powerOfTwo[axis]) % 2;
        // N-linear basis function
        tmp_D =  tmp_D * ((double)(tmp_I2) * axisDiv[axis] + (double)(1 - tmp_I2) * (1.0 - axisDiv[axis]));
        tOffs = tOffs + tmp_I2 * fgCoarse->getOffset(axis);
        //COMBIGRID_OUT_LEVEL2( verb , " Calc basis function, axis:"<< axis << " , tmp_I2:" << tmp_I2 << ",tmp_D:"<<tmp_D);
      }

      // add the contribution from this point
      //COMBIGRID_OUT_LEVEL2( verb , "pInd:"<<pInd<<", basis value:" << tmp_D << ", linearIndexC:" << linearIndexC << ", tOffs:" << tOffs);
      for (d = 0 ; d < nrSpace ; d++ ) {
        val[d] = val[d] + tmp_D * vectCoarse[ nrSpace * (linearIndexC + tOffs) + d ];
      }
    }

    // set the values
    for (i = 0 ; i < nrSpace ; i++ ) {
      vectFine[ nrSpace * pInd + i ] = coefFine * vectFine[ nrSpace * pInd + i ] + coefCoarse * val[i];
    }

    //COMBIGRID_OUT_LEVEL2( verb , " combigrid::ProlongationRestriction::makeLinearProlongation , indexFine:"
    //     << pInd << " , val:"<<val);
  }// end loop over fine grid points

  COMBIGRID_OUT_LEVEL2( verb , "ProlongationRestriction::prolongation ... END");
}


void combigrid::ProlongationRestriction::restriction( const FullGridD* fgFine ,
    const std::vector<double>& vectFine ,
    double coefFine ,
    const FullGridD* fgCoarse ,
    std::vector<double>& vectCoarse ,
    double coefCoarse ,
    int nrSpace) {

  // consistency check
  COMBIGRID_ERROR_TEST((int)vectFine.size() == (int)nrSpace * fgFine->getNrElements(), "ProlongationRestriction::restriction "
                       << " SIZE MUST MATCH vectFine.size():" << vectFine.size() <<
                       " nrSpace:" << nrSpace << " fgFine->getNrElements():" << fgFine->getNrElements() );
  COMBIGRID_ERROR_TEST((int)vectCoarse.size() == (int)nrSpace * fgCoarse->getNrElements(), "ProlongationRestriction::restriction "
                       << " SIZE MUST MATCH vectCoarse.size():" << vectCoarse.size() <<
                       " nrSpace:" << nrSpace << " fgCoarse->getNrElements():" << fgCoarse->getNrElements() );

  // variable declaration
  const int verb = 0;
  const int dim = fgFine->getDimension();
  const int nrCoarsePoints = fgCoarse->getNrElements();
  int tmp_I , linearIndexC , linearIndexF , i , s ;
  std::vector<int> axisIndex(dim);
  std::vector<int> mulFactors(dim, 1);
  std::vector<int> start_offs(dim, 0);
  const std::vector<bool>& boundaryFlag = fgFine->returnBoundaryFlags();

  COMBIGRID_OUT_LEVEL2( verb , "ProlongationRestriction::restriction ... START");

  // first see which axis index (on the fine mesh) have to be multiplied by two
  for (i = 0 ; i < dim ; i++) {
    mulFactors[i] = (fgFine->length(i) + 1) / fgCoarse->length(i); // this is either 1 or 2
    start_offs[i] = 0;

    if (!boundaryFlag[i]) {
      start_offs[i] = (mulFactors[i] > 1) ? 1 : 0;
    }

    COMBIGRID_OUT_LEVEL2( verb , " combigrid::ProlongationRestriction::makeDirectRestriction , mulFactors[i]:" << mulFactors[i] );
  }

  //COMBIGRID_OUT_LEVEL2( verb , " nrCoarsePoints:" << nrCoarsePoints );

  // here we iterate on all coarse grid points
  for (int pInd = 0 ; pInd < nrCoarsePoints ; pInd++) {
    tmp_I = pInd;
    // calc coarse grid point index
    linearIndexF = 0;

    for (i = dim - 1 ; i >= 0 ; i--) {
      axisIndex[i] = (tmp_I / (fgCoarse->getOffset(i)));
      tmp_I = tmp_I % fgCoarse->getOffset(i);
      // calculate the corresponding fine grid point index , restriction starts with point 1
      linearIndexF = linearIndexF + (start_offs[i] + mulFactors[i] * axisIndex[i]) * fgFine->getOffset(i);
      //COMBIGRID_OUT_LEVEL2( verb , " i:" << i << " axisIndex[i]:" << axisIndex[i] << " , mulFactors[i]:" << mulFactors[i] );
      //COMBIGRID_OUT_LEVEL2( verb , " i:" << i << " linearIndexF:" << linearIndexF << " , gFine->offset()[i]:" << gFine->offset()[i] );
    }

    // make the direct restriction by adding the value with a given coefficient
    //COMBIGRID_OUT_LEVEL2( verb , " pInd:" << pInd );
    linearIndexC = pInd;

    for (s = 0 ; s < nrSpace ; s++) {
      vectCoarse[nrSpace * linearIndexC + s] = coefCoarse * vectCoarse[nrSpace * linearIndexC + s] + coefFine * vectFine[nrSpace * linearIndexF + s];
    }

    //COMBIGRID_OUT_LEVEL2( verb , " combigrid::ProlongationRestriction::makeDirectRestriction , linearIndexC:"
    //     << linearIndexC << " , linearIndexF:"<<linearIndexF);
  }// end loop over coarse grid points

  COMBIGRID_OUT_LEVEL2( verb , "ProlongationRestriction::restriction ... END");
}
