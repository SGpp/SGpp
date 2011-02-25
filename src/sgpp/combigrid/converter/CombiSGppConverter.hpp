/*
 * CombiSGppConverter.hpp
 *
 *  Created on: Feb 24, 2011
 *      Author: benk
 */

#ifndef COMBISGPPCONVERTER_HPP_
#define COMBISGPPCONVERTER_HPP_

// ---- the combi includes --------
#include "combigrid/utils/combigrid_ultils.hpp"
#include "combigrid/fullgrid/CombiFullGrid.hpp"
#include "combigrid/combigridkernel/CombiCombiGridKernel.hpp"
#include "combigrid/combigrid/AbstractCombiGrid.hpp"

// ---- the SGpp includes --------
#include "data/DataVector.hpp"
#include "grid/Grid.hpp"
#include "grid/GridStorage.hpp"

using namespace std;
using namespace combigrid;

namespace combigrid {

/** calss which converts a combi grid to a SGpp sparse grid and back <br>
 * */
	class CombiSGppConverter {

public:

		/** empty Ctor*/
		CombiSGppConverter() {;}

		/** creates a SGpp grid out of a combi grid
		 * @param storage the grid storage
		 * @param combikernel contains all the information to create a SGpp grid*/
		static void createSGpp( sg::GridStorage* storage , const CombiGridKernelD* combikernel ){

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

		/** add the value of the FG to the SGpp grid. We assume that the array has been already resized
		 * @param fg the full grid
		 * @param coef the multiplicator coefficient (the combination coefficient)*
		 * @param storage SGpp grid storage
		 * @param alpha coefficient vector */
		static void FullGridToSGpp(const FullGridD* fg , double coef , sg::GridStorage* storage , DataVector *alpha){
			int dim = fg->getDimension() , sgppIndex , k;
			std::vector<int> levelsLI( dim );
			std::vector<int> indexsLI( dim );
			sg::GridIndex *hgi = new sg::GridIndex( dim );

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

		/** similar as the previous method, but it updates the min and max vector <br>
		 * these vectors can be used to estimate the goodnes of the combination of the FGs */
		static void FullGridToSGpp(const FullGridD* fg , double coef , sg::GridStorage* storage ,
				DataVector *alpha , DataVector *minAlpha , DataVector *maxAlpha ){

			int dim = fg->getDimension() , sgppIndex , k;
			std::vector<int> levelsLI( dim );
			std::vector<int> indexsLI( dim );
			sg::GridIndex *hgi = new sg::GridIndex( dim );

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

		/** Sets the values of the FG from the SGpp grid values, this method is used by the decomposition
		 * of the SGpp grid into the combi grid
		 * @param storage SGpp grid storage
		 * @param alpha coefficient vector of the SGpp
		 * @param fg the full grid whos values will be set */
		static void SGppToFullGrid( sg::GridStorage* storage , DataVector *alpha , FullGridD* fg ){
			// for each full grid point get the corresponding value from SGpp and just set the FG value with that value
			int dim = fg->getDimension() , sgppIndex , k;
			std::vector<int> levelsLI( dim , 0.0);
			std::vector<int> indexsLI( dim , 0.0);
			sg::GridIndex *hgi = new sg::GridIndex( dim );

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

	};
}


#endif /* COMBISGPPCONVERTER_HPP_ */
