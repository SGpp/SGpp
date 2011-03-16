/*
 * CombiSGppConverter.cpp
 *
 *  Created on: Feb 24, 2011
 *      Author: benk
 */


#include "combigrid/converter/CombiSGppConverter.hpp"
using namespace sg::base;

using namespace std;
using namespace combigrid;

void CombiSGppConverter::createSGpp( GridStorage* storage , const CombiGridKernelD* combikernel ){

	if(storage->size() > 0)
	{
		storage->emptyStorage();
		//COMBIGRID_ERROR_EXIT("CombiSgppConverter::createSGpp storage not empty , but size = " << storage->size() );
	}

	// create the GridStorage, and insert all points from the each full grid into the hashmap (but only once)
	int dim = combikernel->getDim();
	std::vector<int> levelsLI( dim );
	std::vector<int> indexsLI( dim );
	GridIndex *hgi = new GridIndex( dim );
	int k , sgppIndex ;

	for (int nrfg = 0 ; nrfg < combikernel->getNrFullGrids() ; nrfg++){
		const FullGridD* fg = combikernel->getFullGrid(nrfg);
		// for loop over each full grid points
		for (int nrp = 0 ; nrp < fg->getNrElements() ; nrp++){
			// ... get the index and level
			fg->getLI( nrp , levelsLI , indexsLI);
			//COMBIGRID_OUT_LEVEL3(4 , "createSGpp nrp = " << nrp );
			for (k = 0 ; k < dim ; k++ ){
				hgi->push( k , levelsLI[k] , indexsLI[k] );
				//std::cout << ",l:" <<levelsLI[k] << ",i:" << indexsLI[k];
			}
			//std::cout << std::endl;
			// rehash for this index
			hgi->rehash();
			// test if it already in, if no then add the point
			if ( !(storage->has_key( hgi )) ){
				storage->insert( (*hgi) );
				sgppIndex = (*storage)[hgi];
				//COMBIGRID_OUT_LEVEL3(4 , "createSGpp  add point nrp:" << nrp << " , sgppIndex:" << sgppIndex);
			}
		}
	}
	// delete the dynamically created variable
	delete hgi;
}


void CombiSGppConverter::FullGridToSGpp(const FullGridD* fg , double coef , GridStorage* storage , DataVector *alpha){
	int dim = fg->getDimension() , sgppIndex , k;
	std::vector<int> levelsLI( dim );
	std::vector<int> indexsLI( dim );
	GridIndex *hgi = new GridIndex( dim );

	// each FG value will be added to the DataVector, with the specified coefficient
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
		//COMBIGRID_OUT_LEVEL3(4 , "FullGridToSGpp nrp:" << nrp << " , sgppIndex:" << sgppIndex << " , val:" << fg->getElementVector()[nrp] <<
		//		" , (*alpha)[sgppIndex]:" << (*alpha)[sgppIndex]);
	}
	delete hgi;
}


void CombiSGppConverter::FullGridToSGpp(const FullGridD* fg , double coef , GridStorage* storage ,
				DataVector *alpha , DataVector *minAlpha , DataVector *maxAlpha ){

	int dim = fg->getDimension() , sgppIndex , k;
	std::vector<int> levelsLI( dim );
	std::vector<int> indexsLI( dim );
	GridIndex *hgi = new GridIndex( dim );

	// each FG value will be added to the DataVector, with the specified coefficient
	// as a byproduct the min and maximum Vector will be updated
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
		(*minAlpha)[sgppIndex] = (fg->getElementVector()[nrp] < (*minAlpha)[sgppIndex]) ?
						(fg->getElementVector()[nrp]) : ((*minAlpha)[sgppIndex]);
		(*maxAlpha)[sgppIndex] = (fg->getElementVector()[nrp] > (*maxAlpha)[sgppIndex]) ?
						(fg->getElementVector()[nrp]) : ((*maxAlpha)[sgppIndex]);
		//COMBIGRID_OUT_LEVEL3(4 , "FullGridToSGpp(minmax) nrp:" << nrp << " , sgppIndex:" << sgppIndex << " , val:" << fg->getElementVector()[nrp]);
		//COMBIGRID_OUT_LEVEL3(4 , "FullGridToSGpp (*minAlpha)[sgppIndex]:" << (*minAlpha)[sgppIndex]);
		//COMBIGRID_OUT_LEVEL3(4 , "FullGridToSGpp (*maxAlpha)[sgppIndex]:" << (*maxAlpha)[sgppIndex]);
	}

	delete hgi;
}


void CombiSGppConverter::SGppToFullGrid( GridStorage* storage , DataVector *alpha , FullGridD* fg ){
	// for each full grid point get the corresponding value from SGpp and just set the FG value with that value
	int dim = fg->getDimension() , sgppIndex , k;
	std::vector<int> levelsLI( dim , 0.0);
	std::vector<int> indexsLI( dim , 0.0);
	GridIndex *hgi = new GridIndex( dim );

	// each FG value will be added to the DataVector, with the specified coefficient
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
		fg->getElementVector()[nrp] = (*alpha)[sgppIndex];
	}
	delete hgi;
}
