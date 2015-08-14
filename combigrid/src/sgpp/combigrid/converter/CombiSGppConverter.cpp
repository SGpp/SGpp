// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <sgpp/combigrid/converter/CombiSGppConverter.hpp>

using namespace std;

void combigrid::CombiSGppConverter::createSGpp( SGPP::base::GridStorage* storage , const CombiGridKernelD* combikernel ) {

  if (storage->size() > 0) {
    storage->emptyStorage();
    //COMBIGRID_ERROR_EXIT("CombiSgppConverter::createSGpp storage not empty , but size = " << storage->size() );
  }

  // create the SGPP::base::GridStorage, and insert all points from the each full grid into the hashmap (but only once)
  int dim = combikernel->getDim();
  std::vector<int> levelsLI( dim );
  std::vector<int> indexsLI( dim );
  SGPP::base::GridIndex* hgi = new SGPP::base::GridIndex( dim );
  int k , sgppIndex ;

  // use the stored SGpp index if they exist
  for (int nrfg = 0 ; nrfg < combikernel->getNrFullGrids() ; nrfg++) {
    const FullGridD* fg = combikernel->getFullGrid(nrfg);
    // for loop over each full grid points
    fg->getSGppIndex().resize(fg->getNrElements());

    for (int nrp = 0 ; nrp < fg->getNrElements() ; nrp++) {
      // ... get the index and level
      fg->getLI( nrp , levelsLI , indexsLI);

      //COMBIGRID_OUT_LEVEL3(4 , "createSGpp nrp = " << nrp );
      for (k = 0 ; k < dim ; k++ ) {
        hgi->push( k , levelsLI[k] , indexsLI[k] );
        //std::cout << ",l:" <<levelsLI[k] << ",i:" << indexsLI[k];
      }

      //std::cout << std::endl;
      // rehash for this index
      hgi->rehash();

      // test if it already in, if no then add the point
      if ( !(storage->has_key( hgi )) ) {
        storage->insert( (*hgi) );
        //COMBIGRID_OUT_LEVEL3(4 , "createSGpp  add point nrp:" << nrp << " , sgppIndex:" << sgppIndex);
      }

      sgppIndex = static_cast<int>(storage->seq(hgi));
      fg->getSGppIndex()[nrp] = sgppIndex;
    }
  }

  // delete the dynamically created variable
  delete hgi;
}

// --------- the old and wrong FullGridToSGpp method --------------
/*
 void combigrid::CombiSGppConverter::FullGridToSGpp(const FullGridD* fg , double coef , SGPP::base::GridStorage* storage , SGPP::base::DataVector *alpha){
        int dim = fg->getDimension() , sgppIndex , k;
        std::vector<int> levelsLI( dim );
        std::vector<int> indexsLI( dim );
        SGPP::base::GridIndex *hgi = new SGPP::base::GridIndex( dim );

        // use the stored SGpp index if they exist
        if ( fg->getSGppIndex().size() < 1 )
        {
                // each FG value will be added to the SGPP::base::DataVector, with the specified coefficient
                fg->getSGppIndex().resize(fg->getNrElements());
                for (int nrp = 0 ; nrp < fg->getNrElements() ; nrp++){
                        // ... get the index and level
                        fg->getLI( nrp , levelsLI , indexsLI);
                        for (k = 0 ; k < dim ; k++ ){
                                hgi->push( k , levelsLI[k] , indexsLI[k] );
                        }
                        // rehash for this index
                        hgi->rehash();
                        //todo: we do not test if this is present in the hashmap
                        sgppIndex = (*storage)[hgi];
                        (*alpha)[sgppIndex] = (*alpha)[sgppIndex] + coef * fg->getElementVector()[nrp];
                        fg->getSGppIndex()[nrp] = sgppIndex;
                        //COMBIGRID_OUT_LEVEL3(4 , "FullGridToSGpp nrp:" << nrp << " , sgppIndex:" << sgppIndex << " , val:" << fg->getElementVector()[nrp] <<
                        //              " , (*alpha)[sgppIndex]:" << (*alpha)[sgppIndex]);
                }
        }
        else
        {
                // just get the index
                for (int nrp = 0 ; nrp < fg->getNrElements() ; nrp++){
                        sgppIndex = fg->getSGppIndex()[nrp];
                        (*alpha)[sgppIndex] = (*alpha)[sgppIndex] + coef * fg->getElementVector()[nrp];
                }
        }
        delete hgi;
}

 */

void combigrid::CombiSGppConverter::FullGridToSGpp(const FullGridD* fg , double coef , SGPP::base::GridStorage* storage , SGPP::base::DataVector* alpha) {
  int dim = fg->getDimension();
  SGPP::base::GridIndex* hgi = new SGPP::base::GridIndex( dim );

  // ---- evaluate the full grid at each point of the SGpp and add to the SGpp vector ----

  const GridDomain* domain = fg->getDomain();
  std::vector<double> coords( dim, 0.0), minV(dim, 0.0), scaleV(dim, 1.0);
  double evalVal = 0.0;
  int d;

  if (domain == 0) {
    for (int nrp = 0 ; nrp < (int)storage->size() ; nrp++) {
      // get the coordinates in the unit square
      for (d = 0 ; d < dim ; d++) {
        coords[d] = (*storage)[(size_t)nrp]->getCoordBB( d , 1.0 , 0.0 );
      }

      // evaluate the FG at this position
      evalVal = fg->eval(coords);
      (*alpha)[nrp] = (*alpha)[nrp] + coef * evalVal;
    }
  } else {
    // for
    for (int nrp = 0 ; nrp < (int)storage->size() ; nrp++) {
      for (d = 0 ; d < dim ; d++) {
        coords[d] = (*storage)[(size_t)nrp]->getCoordBB( d , 1.0 , 0.0 );
        //fg->getDomain()->get1DDomain(d).getMaxDomain() - fg->getDomain()->get1DDomain(d).getMinDomain() ,
        //fg->getDomain()->get1DDomain(d).getMinDomain() );
        // transform the coordinates
        domain->get1DDomain(d).transformRealToUnit( coords[d] , coords[d] );
      }

      // evaluate the FG at this position
      evalVal = fg->eval(coords);
      (*alpha)[nrp] = (*alpha)[nrp] + coef * evalVal;
    }
  }

  delete hgi;
}


void combigrid::CombiSGppConverter::FullGridToSGpp(const FullGridD* fg , double coef , SGPP::base::GridStorage* storage ,
    SGPP::base::DataVector* alpha , SGPP::base::DataVector* minAlpha , SGPP::base::DataVector* maxAlpha ) {

  int dim = fg->getDimension();
  SGPP::base::GridIndex* hgi = new SGPP::base::GridIndex( dim );

  // ---- evaluate the full grid at each point of the SGpp and add to the SGpp vector ----

  const GridDomain* domain = fg->getDomain();
  std::vector<double> coords( dim, 0.0), minV(dim, 0.0), scaleV(dim, 1.0);
  double evalVal = 0.0;
  int d;

  if (domain == 0) {
    for (int nrp = 0 ; nrp < (int)storage->size() ; nrp++) {
      // get the coordinates in the unit square
      for (d = 0 ; d < dim ; d++) {
        coords[d] = (*storage)[(size_t)nrp]->getCoordBB( d , 1.0 , 0.0 );
      }

      // evaluate the FG at this position
      evalVal = fg->eval(coords);
      (*alpha)[nrp] = (*alpha)[nrp] + coef * evalVal;
      (*minAlpha)[nrp] = (evalVal < (*minAlpha)[nrp]) ?
                         (evalVal) : ((*minAlpha)[nrp]);
      (*maxAlpha)[nrp] = (evalVal > (*maxAlpha)[nrp]) ?
                         (evalVal) : ((*maxAlpha)[nrp]);
      //COMBIGRID_OUT_LEVEL3(4 , "FullGridToSGpp1(minmax) nrp:" << nrp << " , val:" << evalVal);
      //COMBIGRID_OUT_LEVEL3(4 , "FullGridToSGpp (*minAlpha)[nrp]:" << (*minAlpha)[nrp]);
      //COMBIGRID_OUT_LEVEL3(4 , "FullGridToSGpp (*maxAlpha)[nrp]:" << (*maxAlpha)[nrp]);
    }
  } else {
    // for
    for (int nrp = 0 ; nrp < (int)storage->size() ; nrp++) {
      for (d = 0 ; d < dim ; d++) {
        coords[d] = (*storage)[(size_t)nrp]->getCoordBB( d , 1.0 , 0.0 );
        //fg->getDomain()->get1DDomain(d).getMaxDomain() - fg->getDomain()->get1DDomain(d).getMinDomain() ,
        //fg->getDomain()->get1DDomain(d).getMinDomain() );
        // transform the coordinates
        domain->get1DDomain(d).transformRealToUnit( coords[d] , coords[d] );
      }

      // evaluate the FG at this position
      evalVal = fg->eval(coords);
      (*alpha)[nrp] = (*alpha)[nrp] + coef * evalVal;
      (*minAlpha)[nrp] = (evalVal < (*minAlpha)[nrp]) ?
                         (evalVal) : ((*minAlpha)[nrp]);
      (*maxAlpha)[nrp] = (evalVal > (*maxAlpha)[nrp]) ?
                         (evalVal) : ((*maxAlpha)[nrp]);
      //COMBIGRID_OUT_LEVEL3(4 , "FullGridToSGpp2(minmax) nrp:" << nrp << " , val:" << evalVal);
      //COMBIGRID_OUT_LEVEL3(4 , "FullGridToSGpp (*minAlpha)[nrp]:" << (*minAlpha)[nrp]);
      //COMBIGRID_OUT_LEVEL3(4 , "FullGridToSGpp (*maxAlpha)[nrp]:" << (*maxAlpha)[nrp]);
    }
  }

  delete hgi;
}


void combigrid::CombiSGppConverter::SGppToFullGrid( SGPP::base::GridStorage* storage , SGPP::base::DataVector* alpha , FullGridD* fg ) {
  // for each full grid point get the corresponding value from SGpp and just set the FG value with that value
  int dim = fg->getDimension() , sgppIndex , k;
  std::vector<int> levelsLI( dim , 0.0);
  std::vector<int> indexsLI( dim , 0.0);
  SGPP::base::GridIndex* hgi = new SGPP::base::GridIndex( dim );

  // use the stored SGpp index if they exist
  if ( fg->getSGppIndex().size() < 1 ) {
    fg->getSGppIndex().resize(fg->getNrElements());

    // each FG value will be added to the SGPP::base::DataVector, with the specified coefficient
    for (int nrp = 0 ; nrp < fg->getNrElements() ; nrp++) {
      // ... get the index and level
      fg->getLI( nrp , levelsLI , indexsLI);

      for (k = 0 ; k < dim ; k++ ) {
        hgi->push( k , levelsLI[k] , indexsLI[k] );
      }

      // rehash for this index
      hgi->rehash();
      //todo: we do not test if this is present in the hashmap
      sgppIndex = static_cast<int>(storage->seq(hgi));
      fg->getElementVector()[nrp] = (*alpha)[sgppIndex];
      fg->getSGppIndex()[nrp] = sgppIndex;
    }
  } else {
    // just get the index from the full grid
    for (int nrp = 0 ; nrp < fg->getNrElements() ; nrp++) {
      // get the SGpp index
      sgppIndex = fg->getSGppIndex()[nrp];
      fg->getElementVector()[nrp] = (*alpha)[sgppIndex];
    }
  }

  delete hgi;
}
