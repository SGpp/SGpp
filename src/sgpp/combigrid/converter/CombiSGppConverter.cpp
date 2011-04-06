/*
 * CombiSGppConverter.cpp
 *
 *  Created on: Feb 24, 2011
 *      Author: benk
 */


#include "combigrid/converter/CombiSGppConverter.hpp"

using namespace std;
using namespace combigrid;

void CombiSGppConverter::createSGpp( sg::GridStorage* storage , const CombiGridKernelD* combikernel ){

	if(storage->size() > 0)
	{
		storage->emptyStorage();
		//COMBIGRID_ERROR_EXIT("CombiSgppConverter::createSGpp storage not empty , but size = " << storage->size() );
	}

	// create the GridStorage, and insert all points from the each full grid into the hashmap (but only once)
	int dim = combikernel->getDim();
	std::vector<int> levelsLI( dim );
	std::vector<int> indexsLI( dim );
	sg::GridIndex *hgi = new sg::GridIndex( dim );
	int k , sgppIndex ;

	// use the stored SGpp index if they exist
	for (int nrfg = 0 ; nrfg < combikernel->getNrFullGrids() ; nrfg++){
		const FullGridD* fg = combikernel->getFullGrid(nrfg);
		// for loop over each full grid points
		fg->getSGppIndex().resize(fg->getNrElements());
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
				//COMBIGRID_OUT_LEVEL3(4 , "createSGpp  add point nrp:" << nrp << " , sgppIndex:" << sgppIndex);
			}
			sgppIndex = (*storage)[hgi];
			fg->getSGppIndex()[nrp] = sgppIndex;
		}
	}
	// delete the dynamically created variable
	delete hgi;
}


void CombiSGppConverter::FullGridToSGpp(const FullGridD* fg , double coef , sg::GridStorage* storage , DataVector *alpha){
	int dim = fg->getDimension() , sgppIndex , k;
	std::vector<int> levelsLI( dim );
	std::vector<int> indexsLI( dim );
	sg::GridIndex *hgi = new sg::GridIndex( dim );

	// use the stored SGpp index if they exist
	if ( fg->getSGppIndex().size() < 1 )
	{
		// each FG value will be added to the DataVector, with the specified coefficient
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
			//		" , (*alpha)[sgppIndex]:" << (*alpha)[sgppIndex]);
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


void CombiSGppConverter::FullGridToSGpp(const FullGridD* fg , double coef , sg::GridStorage* storage ,
				DataVector *alpha , DataVector *minAlpha , DataVector *maxAlpha ){

	int dim = fg->getDimension() , sgppIndex , k;
	std::vector<int> levelsLI( dim );
	std::vector<int> indexsLI( dim );
	sg::GridIndex *hgi = new sg::GridIndex( dim );

	// use the stored SGpp index if they exist
	if ( fg->getSGppIndex().size() < 1 ){
		// each FG value will be added to the DataVector, with the specified coefficient
		// as a byproduct the min and maximum Vector will be updated
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
			(*minAlpha)[sgppIndex] = (fg->getElementVector()[nrp] < (*minAlpha)[sgppIndex]) ?
						(fg->getElementVector()[nrp]) : ((*minAlpha)[sgppIndex]);
			(*maxAlpha)[sgppIndex] = (fg->getElementVector()[nrp] > (*maxAlpha)[sgppIndex]) ?
						(fg->getElementVector()[nrp]) : ((*maxAlpha)[sgppIndex]);
			//COMBIGRID_OUT_LEVEL3(4 , "FullGridToSGpp(minmax) nrp:" << nrp << " , sgppIndex:" << sgppIndex << " , val:" << fg->getElementVector()[nrp]);
			//COMBIGRID_OUT_LEVEL3(4 , "FullGridToSGpp (*minAlpha)[sgppIndex]:" << (*minAlpha)[sgppIndex]);
			//COMBIGRID_OUT_LEVEL3(4 , "FullGridToSGpp (*maxAlpha)[sgppIndex]:" << (*maxAlpha)[sgppIndex]);
			fg->getSGppIndex()[nrp] = sgppIndex;
		}
	}
	else
	{
		for (int nrp = 0 ; nrp < fg->getNrElements() ; nrp++){
			// get the SGpp index
			sgppIndex = fg->getSGppIndex()[nrp];
			// use the index
			(*alpha)[sgppIndex] = (*alpha)[sgppIndex] + coef * fg->getElementVector()[nrp];
			(*minAlpha)[sgppIndex] = (fg->getElementVector()[nrp] < (*minAlpha)[sgppIndex]) ?
						(fg->getElementVector()[nrp]) : ((*minAlpha)[sgppIndex]);
			(*maxAlpha)[sgppIndex] = (fg->getElementVector()[nrp] > (*maxAlpha)[sgppIndex]) ?
						(fg->getElementVector()[nrp]) : ((*maxAlpha)[sgppIndex]);
		}
	}

	delete hgi;
}


void CombiSGppConverter::SGppToFullGrid( sg::GridStorage* storage , DataVector *alpha , FullGridD* fg ){
	// for each full grid point get the corresponding value from SGpp and just set the FG value with that value
	int dim = fg->getDimension() , sgppIndex , k;
	std::vector<int> levelsLI( dim , 0.0);
	std::vector<int> indexsLI( dim , 0.0);
	sg::GridIndex *hgi = new sg::GridIndex( dim );

	// use the stored SGpp index if they exist
	if ( fg->getSGppIndex().size() < 1 )
	{
		fg->getSGppIndex().resize(fg->getNrElements());
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
			fg->getSGppIndex()[nrp] = sgppIndex;
		}
	}
	else
	{
		// just get the index from the full grid
		for (int nrp = 0 ; nrp < fg->getNrElements() ; nrp++){
			// get the SGpp index
			sgppIndex = fg->getSGppIndex()[nrp];
			fg->getElementVector()[nrp] = (*alpha)[sgppIndex];
		}
	}
	delete hgi;
}
