// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <sgpp/combigrid/multigridFG/operators/PoissonOperator.hpp>
#include <math.h>


combigrid::PoissonOperator::PoissonOperator(const FullGridD* fg ,
    const std::vector<double>& sigma ,
    const combigrid::CallBackRHS* callbackRHS) : combigrid::OperatorFG(fg , 1) , sigma_(sigma) , callbackRHS_(callbackRHS) {

  const int verb = 0;
  COMBIGRID_OUT_LEVEL2( verb , "PoissonOperator::PoissonOperator START");
  // the full grid must have boundary points
  // consistency check for the grid (boundary points, at least 3 points per axis)
  // - the grid must have a domain, even if it is a unit square
  COMBIGRID_ERROR_TEST( getFG()->getDomain() != NULL , " FULL grid must have a valid domain , even it is unit square");
  dim_ = fg->getDimension();

  for (int d = 0 ; d < dim_ ; d++ ) {
    COMBIGRID_OUT_LEVEL2( verb , "PoissonOperator::PoissonOperator test level of dimension d = " << d );
    COMBIGRID_ERROR_TEST( getFG()->getLevels()[d] >= 1 ,
                          " Each dimension must have at least level 2 getFG()->getLevels()[" << d << "]=" << getFG()->getLevels()[d]);
    COMBIGRID_OUT_LEVEL2( verb , "PoissonOperator::PoissonOperator test flag = " << d );
    COMBIGRID_ERROR_TEST( getFG()->returnBoundaryFlags()[d] == true , " The full grid must have in each dimension boundary points d=" << d);
  }

  // the number of unknowns is the number of nodes in the full grids
  nrElem_ = fg->getNrElements();
  row_ptr_.resize(nrElem_ + 1, 0);
  rhs_.resize(nrElem_, 0.0);
  // first count the number of possible non zero elements of the
  nrMatrixElem_ = 0;
  COMBIGRID_OUT_LEVEL2( verb , "PoissonOperator::PoissonOperator dim_=" << dim_ << " , nrElem_=" << nrElem_ );
  COMBIGRID_OUT_LEVEL2( verb , "PoissonOperator::PoissonOperator count the number of possible non-zero matrix entries");
  int e = -2 , i = -2 , tmp = -2;
  std::vector<int> axisIndex(dim_, 0);
  std::vector<bool> isBoundaryPoint( nrElem_ , false );
  std::vector<double> coords(dim_, 0.0);

  for (e = 0 ; e < nrElem_ ; e++) {
    tmp = e;

    // get the index for each axis
    for (i = dim_ - 1 ; i >= 0 ; i--) {
      axisIndex[i] = tmp / (getFG()->getOffset(i));
      tmp = tmp % getFG()->getOffset(i);
    }

    // after we decomposed the index, we see how many overlaping
    tmp = combigrid::powerOfTwo[dim_] + 1;

    for (i = 0 ; i < dim_ ; i++) {
      if ( (axisIndex[i] <= 0) || (axisIndex[i] >= getFG()->length(i) - 1)) {
        tmp = 1;
        isBoundaryPoint[e] = true;
      }
    }

    nrMatrixElem_ = nrMatrixElem_ + tmp;
  }

  // now we know the number of possible non-zero elements in the matrix and so we can create the CRS structure
  matrixVal_.resize( nrMatrixElem_ , 0.0);
  col_ind_.resize( nrMatrixElem_ , -2);
  COMBIGRID_OUT_LEVEL2( verb , "PoissonOperator::PoissonOperator set up the CRS structure nrMatrixElem_=" << nrMatrixElem_);
  int actualP = 0 , nrElPerRow = -2;

  // loop over each node, to create first the CRS structure
  for (e = 0 ; e < nrElem_ ; e++) {
    tmp = e;

    // get the index for each axis
    for (i = dim_ - 1 ; i >= 0 ; i--) {
      axisIndex[i] = tmp / (getFG()->getOffset(i));
      tmp = tmp % getFG()->getOffset(i);
      COMBIGRID_OUT_LEVEL5( verb , "axisIndex[" << i << "]=" << axisIndex[i] << "  tmp = " << tmp);
    }

    row_ptr_[e] = actualP;
    COMBIGRID_OUT_LEVEL5( verb , " row_ptr_[" << e << "] = " << actualP );
    nrElPerRow = 0;

    // add elements to the matrix
    if (!isBoundaryPoint[e]) {
      for (i = 0 ; i < dim_ ; i++) {
        col_ind_[actualP] = e - getFG()->getOffset(i);
        COMBIGRID_OUT_LEVEL5( verb , " i = " << i << " col_ind_[" << actualP << "]=" << col_ind_[actualP] );
        actualP = actualP + 1;
        nrElPerRow = nrElPerRow + 1;
      }
    }

    col_ind_[actualP] = e;
    COMBIGRID_OUT_LEVEL5( verb , " i = " << i << " col_ind_[" << actualP << "]=" << col_ind_[actualP] );
    actualP = actualP + 1;
    nrElPerRow = nrElPerRow + 1;

    if (!isBoundaryPoint[e]) {
      for (i = dim_ - 1 ; i >= 0 ; i--) {
        col_ind_[actualP] = e + getFG()->getOffset(i);
        COMBIGRID_OUT_LEVEL5( verb , " i = " << i << " col_ind_[" << actualP << "]=" << col_ind_[actualP] );
        actualP = actualP + 1;
        nrElPerRow = nrElPerRow + 1;
      }
    }
  }

  // the last element is the end
  row_ptr_[nrElem_] = row_ptr_[nrElem_ - 1] + nrElPerRow;
  COMBIGRID_OUT_LEVEL2( verb , "LAST row_ptr_[" << nrElem_ << "] = " << row_ptr_[nrElem_] );

  //
  int p_e , p_e_m , p_e_p ; // the positions of the nodes
  double h0 = 1.0 , h1 = 1.0;
  COMBIGRID_OUT_LEVEL2( verb , "PoissonOperator::PoissonOperator  add the laplace ");

  for ( e = 0 ; e < nrElem_ ; e++) {
    // at the boundary we have constant Dirichlet values
    if (isBoundaryPoint[e]) {
      COMBIGRID_OUT_LEVEL2( verb , "Row " << e << " is constant");
      matrixVal_[row_ptr_[e]] = 1.0;

      // constant values at the boundary
      if ( fg->getElementVector().size() > 0) {
        rhs_[e] = fg->getElementVector()[e];
      } else {
        rhs_[e] = 0.0;
      }
    } else {
      // here we look for the BBT(i,i) element,
      p_e = row_ptr_[ e ] + (row_ptr_[ e + 1] - row_ptr_[ e ] - 1) / 2;
      COMBIGRID_OUT_LEVEL5( verb , "element position e=" << e << " , p_e=" << p_e );

      // set right hand side
      getFG()->getCoords( e , coords );
      rhs_[e] = callbackRHS_->eval( coords );
      tmp = e;

      for (i = dim_ - 1 ; i >= 0 ; i--) {
        axisIndex[i] = tmp / (getFG()->getOffset(i));
        tmp = tmp % getFG()->getOffset(i);
        // add to the direction only if it is not boundary direction
        COMBIGRID_OUT_LEVEL5( verb , "axisIndex[" << i << "]=" << axisIndex[i] << " , getFG()->length(i)=" << getFG()->length(i));

        // look for e - getFG()->getOffset(i) and for e+ getFG()->getOffset(i) elements
        p_e_m = p_e - (dim_ - i);
        p_e_p = p_e + (dim_ - i);

        // get h0 and h1
        // h0 = S_scaling(1,i) - S_scaling(1,i-1);
        // h1 = S_scaling(1,i+1) - S_scaling(1,i);
        getFG()->getDomain()->get1DDomain(i).getMeshWidth( axisIndex[i] , getFG()->getLevels()[i] , h0 , h1 );
        COMBIGRID_OUT_LEVEL5( verb , " add diffusion , sigma_[" << i << "]=" << sigma_[i] << " , h0=" << h0
                              << " , h1=" << h1 << " , p_e_m=" << p_e_m << " , p_e=" << p_e << " , p_e_p=" << p_e_p);

        //  here we add the laplace operator
        //  C(ind,ind-nr_points(2)) = C(ind,ind-nr_points(2)) - h1/(0.5*h0*h1*(h0+h1));
        //  C(ind,ind) = C(ind,ind) + (h0+h1)/(0.5*h0*h1*(h0+h1));
        //  C(ind,ind+nr_points(2)) = C(ind,ind+nr_points(2)) - h0/(0.5*h0*h1*(h0+h1));
        matrixVal_[p_e_m] = matrixVal_[p_e_m] - sigma_[i] * h1 / (0.5 * h0 * h1 * (h0 + h1));
        matrixVal_[p_e]   = matrixVal_[p_e] + sigma_[i] * (h0 + h1) / (0.5 * h0 * h1 * (h0 + h1));
        matrixVal_[p_e_p] = matrixVal_[p_e_p] - sigma_[i] * h0 / (0.5 * h0 * h1 * (h0 + h1));
      }
    }
  }

  if (verb > 3) {
    // visialize the matrix and the right hand side
    int startSearch = 0 , endSearch = 0;
    tmp = 0;
    COMBIGRID_OUT_LEVEL5( verb , " ============ START MATRIX PLOTT ============ ");

    for ( e = 0 ; e < nrElem_ ; e++) {
      startSearch = row_ptr_[ e ];
      endSearch = row_ptr_[ e + 1];

      // for each element in the row
      for ( ; startSearch < endSearch ; startSearch++) {
        COMBIGRID_OUT_LEVEL5( verb , " A(" << e + 1 << "," << col_ind_[tmp] + 1 << ") = " << matrixVal_[tmp] << ";");
        tmp++;
      }
    }

    for ( e = 0 ; e < nrElem_ ; e++) {
      COMBIGRID_OUT_LEVEL5( verb , "b(" << e + 1 << ")= " << rhs_[e] << ";");
    }

    COMBIGRID_OUT_LEVEL5( verb , " ============ END MATRIX PLOTT ============ ");
  }
}

combigrid::PoissonOperator::~PoissonOperator() {
  // nothing to do, no object created on the heap
}

combigrid::OperatorFG* combigrid::PoissonOperator::factory(const FullGridD* fg) const {
  return (new combigrid::PoissonOperator( fg , sigma_ , callbackRHS_ ));
}

void combigrid::PoissonOperator::getRHS(std::vector<double>& rhs , int& nrSpace) const {
  nrSpace = this->getNrSpace();
  rhs.resize(nrElem_);
  COMBIGRID_OUT_LEVEL2( 0 , "PoissonOperator::getRHS ... START");

  for (int i = 0 ; i < nrElem_ ; i++)  {
    rhs[i] = rhs_[i];
  }

  COMBIGRID_OUT_LEVEL2( 0 , "PoissonOperator::getRHS ... END");
}

void combigrid::PoissonOperator::multiplyVector(std::vector<double>& inVect , std::vector<double>& outVect) const {
  // implement this
  int i , of , verb = 0;
  double tmp = 0.0;
  COMBIGRID_OUT_LEVEL2( verb , "PoissonOperator::multiplyVector ... START");

  for (i = 0 ; i < nrElem_ ; i++) {
    tmp = 0.0;

    for ( of = row_ptr_[ i ] ; of < row_ptr_[ i + 1] ; of++ ) {
      tmp = tmp + matrixVal_[ of ] * inVect[ col_ind_[of] ];
    }

    outVect[i] = tmp;
  }

  COMBIGRID_OUT_LEVEL2( verb , "PoissonOperator::multiplyVector ... END");
}

void combigrid::PoissonOperator::doSmoothing(int nrIt ,
    std::vector<double>& u, std::vector<double>& rhs) const {
  // make a simple Gauss-Seidel
  int i , of , aii , verb = 0;
  double tmp = 0.0;
  COMBIGRID_OUT_LEVEL2( verb , "PoissonOperator::doSmoothing ... START it:" << nrIt);

  for (int it = 0 ; it < nrIt ; it++) {
    // one forward
    for (i = 0 ; i < nrElem_ ; i++) {
      tmp = 0.0;
      aii = -2;

      for ( of = row_ptr_[i] ; of < row_ptr_[i + 1] ; of++ ) {
        tmp = tmp + matrixVal_[of] * u[col_ind_[of]];

        if ( col_ind_[of] == i ) {
          aii = of;
        }
      }

      // because this matrix is symmetric
      u[i] = ( 1.0 / matrixVal_[aii] ) * ( rhs[i] - tmp + u[i] * matrixVal_[aii] );
    }

    /*// one backward
    for (i = nrElem_-1 ; i >= 0 ; i--)
    {
      tmp = 0.0;
      aii = -2;
      for ( of = row_ptr_[i] ; of < row_ptr_[i + 1] ; of++ )
      {
        tmp = tmp + matrixVal_[of] * u[col_ind_[of]];
        if ( col_ind_[of] == i ) { aii = of; }
      }
      // because this matrix is symmetric
      u[i] = ( 1.0/matrixVal_[aii] ) * ( rhs[i] - tmp + u[i]*matrixVal_[aii] );
    }*/
  }

  COMBIGRID_OUT_LEVEL2( verb , "PoissonOperator::doSmoothing ... END");
}