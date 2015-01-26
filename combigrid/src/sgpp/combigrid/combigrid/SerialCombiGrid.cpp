/* ****************************************************************************
* Copyright (C) 2011 Technische Universitaet Muenchen                         *
* This file is part of the SG++ project. For conditions of distribution and   *
* use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
**************************************************************************** */
// @author Janos Benk (benk@in.tum.de)
// @author Christoph Kowitz (kowitz@in.tum.de)


#include <sgpp/combigrid/combigrid/SerialCombiGrid.hpp>

using namespace std;



void combigrid::SerialCombiGrid::createFullGrids() {
  // iterate over each grid and create the all the full grid
  for ( int i = 0 ; i < combikernel_->getNrFullGrids() ; i++) {
    combikernel_->getFullGrid(i)->createFullGrid();
  }
}


void combigrid::SerialCombiGrid::setDomainAllFG( combigrid::GridDomain* gridDomain ) const {
  // sets the domain for all fullgrids
  for ( int i = 0 ; i < combikernel_->getNrFullGrids() ; i++) {
    combikernel_->getFullGrid(i)->setDomain(gridDomain);
  }
}


double combigrid::SerialCombiGrid::eval( std::vector<double>& coords ) const {
  double result = 0.0;
  double coef = 0.0;
  std::vector<double> coords_tmp = coords;

  //COMBIGRID_OUT_LEVEL3( 4 , "SerialCombiGrid::eval");
  // if there is a transformation then transform to the unit coordinates
  if (gridDomain_ != NULL ) {
    gridDomain_->transformRealToUnit( coords , combischeme_->getMaxLevel() , this->getBoundaryFlags() );
  }

  // we evaluate each full grid and multiply with the coefficient, and sum the result up
  for ( int i = 0 ; i < combikernel_->getNrFullGrids() ; i++) {
    coords_tmp = coords;
    coef = combikernel_->getCoef(i);

    if (coef != 0.0) {
      result = result + combikernel_->getCoef(i) * combikernel_->getFullGrid(i)->eval(coords_tmp);
      //      std::cout<<coords[2]<<'\t'<<combikernel_->getCoef(i)<<'\t'<<combikernel_->getFullGrid(i)->eval(coords_tmp)<<std::endl;
    }
  }

  //COMBIGRID_OUT_LEVEL3( 4 , "SerialCombiGrid::eval result=" << result);
  return result;
}

double combigrid::SerialCombiGrid::evalSingleGrid(int index, std::vector<double>& coords) const {
  std::vector<double> coords_temp = coords;

  if (gridDomain_ != NULL ) {
    gridDomain_->transformRealToUnit( coords , combischeme_->getMaxLevel() , this->getBoundaryFlags() );
  }

  return combikernel_->getFullGrid(index)->eval(coords_temp);
}

void combigrid::SerialCombiGrid::eval( std::vector< std::vector<double> >& coords , std::vector<double>& results ) const {
  // just iterate over each point and call the serial evaluation function
  for ( int i = 0 ; i < (int)results.size() ; i++) {
    results[i] = eval(coords[i]);
  }
}

SGPP::base::GridStorage* combigrid::SerialCombiGrid::createSGppGridStorage() const {
  SGPP::base::GridStorage* gridStoreSGpp = new SGPP::base::GridStorage( combikernel_->getDim() );
  combigrid::CombiSGppConverter::createSGpp( gridStoreSGpp , combikernel_ );
  return gridStoreSGpp;
}

void combigrid::SerialCombiGrid::reCompose(SGPP::base::GridStorage* gridstorageSGpp , SGPP::base::DataVector* alpha,
    SGPP::base::DataVector* minAlpha , SGPP::base::DataVector* maxAlpha ) const {
  // for each full grid call the converter
  if (minAlpha != NULL && maxAlpha != NULL) {
    // min and max will be computed
    for ( int i = 0 ; i < (int)alpha->getSize() ; i++) {
      // we set the vector to zero and the min and max values
      (*alpha)[i] = 0.0;
      (*minAlpha)[i] = 1e+100;
      (*maxAlpha)[i] = -1e+100;
    }

    // for each full grid just call the converter function which makes the conversion automatically
    for ( int i = 0 ; i < combikernel_->getNrFullGrids() ; i++) {
      const FullGridD* fg = combikernel_->getFullGrid(i);
      combigrid::CombiSGppConverter::FullGridToSGpp( fg , combikernel_->getCoef(i) , gridstorageSGpp , alpha , minAlpha , maxAlpha );
    }
  } else {
    // this is the case when no min or max should be calculated
    for ( int i = 0 ; i < (int)alpha->getSize() ; i++) {
      // we set the vector to zero
      (*alpha)[i] = 0.0;
    }

    // for each full grid just call the converter function which makes the conversion automatically
    for ( int i = 0 ; i < combikernel_->getNrFullGrids() ; i++) {
      const FullGridD* fg = combikernel_->getFullGrid(i);
      combigrid::CombiSGppConverter::FullGridToSGpp( fg , combikernel_->getCoef(i) , gridstorageSGpp , alpha );
    }
  }
}

void combigrid::SerialCombiGrid::deCompose(SGPP::base::GridStorage* gridstorageSGpp , SGPP::base::DataVector* alpha) {
  for ( int i = 0 ; i < combikernel_->getNrFullGrids() ; i++) {
    FullGridD* fg = combikernel_->getFullGrid(i);
    combigrid::CombiSGppConverter::SGppToFullGrid( gridstorageSGpp , alpha , fg );
  }
}
