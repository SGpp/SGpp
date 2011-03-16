/*
 * SerialCombiGrid.hpp
 *
 *  Created on: Feb 24, 2011
 *      Author: benk
 */

#include "combigrid/combigrid/SerialCombiGrid.hpp"
using namespace sg::base;

using namespace std;

using namespace combigrid;


void SerialCombiGrid::createFullGrids() {
	// iterate over each grid and create the all the full grid
	for ( int i = 0 ; i < combikernel_->getNrFullGrids() ; i++){
		combikernel_->getFullGrid(i)->createFullGrid();
	}
}

double SerialCombiGrid::eval( const std::vector<double>& coords ) const {
	double result = 0.0;
	// we evaluate each full grid and multiply with the coefficient, and sum the result up
	for ( int i = 0 ; i < combikernel_->getNrFullGrids() ; i++){
		result = result + combikernel_->getCoef(i) * combikernel_->getFullGrid(i)->eval(coords);
	}
	return result;
}

void SerialCombiGrid::eval( const std::vector< std::vector<double> >& coords , std::vector<double>& results ) const {
	// just iterate over each point and call the serial evaluation function
	for ( int i = 0 ; i < (int)results.size() ; i++){
		results[i] = eval(coords[i]);
	}
}

GridStorage* SerialCombiGrid::createSGppGridStorage() const {
	GridStorage* gridStoreSGpp = new GridStorage( combikernel_->getDim() );
	CombiSGppConverter::createSGpp( gridStoreSGpp , combikernel_ );
	return gridStoreSGpp;
}

void SerialCombiGrid::reCompose(GridStorage* gridstorageSGpp , DataVector* alpha,
				DataVector* minAlpha , DataVector* maxAlpha ) const {
	// for each full grid call the converter
	if (minAlpha != NULL && maxAlpha != NULL){
		// min and max will be computed
		for ( int i = 0 ; i < (int)alpha->getSize() ; i++){
			// we set the vector to zero and the min and max values
			(*alpha)[i] = 0.0;
			(*minAlpha)[i] = 1e+100;
			(*maxAlpha)[i] = -1e+100;
		}
		// for each full grid just call the converter function which makes the conversion automatically
		for ( int i = 0 ; i < combikernel_->getNrFullGrids() ; i++){
			const FullGridD* fg = combikernel_->getFullGrid(i);
			CombiSGppConverter::FullGridToSGpp( fg , combikernel_->getCoef(i) , gridstorageSGpp , alpha , minAlpha , maxAlpha );
		}
	}
	else{ // this is the case when no min or max should be calculated
		for ( int i = 0 ; i < (int)alpha->getSize() ; i++){
			// we set the vector to zero
			(*alpha)[i] = 0.0;
		}
		// for each full grid just call the converter function which makes the conversion automatically
		for ( int i = 0 ; i < combikernel_->getNrFullGrids() ; i++){
			const FullGridD* fg = combikernel_->getFullGrid(i);
			CombiSGppConverter::FullGridToSGpp( fg , combikernel_->getCoef(i) , gridstorageSGpp , alpha );
		}
	}
}

void SerialCombiGrid::deCompose(GridStorage* gridstorageSGpp , DataVector* alpha) {
	for ( int i = 0 ; i < combikernel_->getNrFullGrids() ; i++){
		FullGridD* fg = combikernel_->getFullGrid(i);
		CombiSGppConverter::SGppToFullGrid( gridstorageSGpp , alpha , fg );
	}
}
