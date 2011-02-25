/*
 * SerialCombiGrid.hpp
 *
 *  Created on: Feb 24, 2011
 *      Author: benk
 */

#ifndef SERIALCOMBIGRID_HPP_
#define SERIALCOMBIGRID_HPP_

#include "combigrid/combigrid/AbstractCombiGrid.hpp"
#include "combigrid/converter/CombiSGppConverter.hpp"

using namespace std;

namespace combigrid{

	/** Simple serial implementation of the combi grid. */
	class SerialCombiGrid : public AbstractCombiGrid{

	public:

		/** Ctor as input argument a CombiScheme is needed */
		SerialCombiGrid(const CombiSchemeBasis* combischeme ,
				const std::vector<bool>& hasBoundaryPts ) : AbstractCombiGrid( combischeme , hasBoundaryPts) { ; }

		/** Ctor as input argument a CombiScheme is needed */
		SerialCombiGrid(const CombiSchemeBasis* combischeme , bool hasBoundaryPts = true ) :
			AbstractCombiGrid( combischeme , hasBoundaryPts ) { ; }

		/** see supercalss for docu */
		virtual void createFullGrids() {
			// iterate over each grid and create the all the full grid
			for ( int i = 0 ; i < combikernel_->getNrFullGrids() ; i++){
				combikernel_->getFullGrid(i)->createFullGrid();
			}
		}

		/** see supercalss for docu */
		virtual FullGridD* getFullGrid( int i ) { return combikernel_->getFullGrid(i); }

		/** see supercalss for docu */
		virtual const FullGridD* getFullGrid( int i ) const { return combikernel_->getFullGrid(i); }

		/** see supercalss for docu */
		virtual int getNrFullGrid() const { return combikernel_->getNrFullGrids(); }

		/** see supercalss for docu */
		virtual double eval( const std::vector<double>& coords ) const {
			double result = 0.0;
			// we evaluate each full grid and multiply with the coefficient, and sum the result up
			for ( int i = 0 ; i < combikernel_->getNrFullGrids() ; i++){
				result = result + combikernel_->getCoef(i) * combikernel_->getFullGrid(i)->eval(coords);
			}
			return result;
		}

		/** see supercalss for docu */
		virtual void eval( const std::vector< std::vector<double> >& coords , std::vector<double>& results ) const {
			// just iterate over each point and call the serial evaluation function
			for ( int i = 0 ; i < (int)results.size() ; i++){
				results[i] = eval(coords[i]);
			}
		}

		/** see supercalss for docu */
		virtual sg::GridStorage* createSGppGridStorage() const {
			sg::GridStorage* gridStoreSGpp = new sg::GridStorage( combikernel_->getDim() );
			CombiSGppConverter::createSGpp( gridStoreSGpp , combikernel_ );
			return gridStoreSGpp;
		}

		/** see supercalss for docu */
		virtual void reCompose(sg::GridStorage* gridstorageSGpp , DataVector* alpha,
				DataVector* minAlpha = NULL , DataVector* maxAlpha = NULL) const {
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

		/** see supercalss for docu */
		virtual void deCompose(sg::GridStorage* gridstorageSGpp , DataVector* alpha) {
			for ( int i = 0 ; i < combikernel_->getNrFullGrids() ; i++){
				FullGridD* fg = combikernel_->getFullGrid(i);
				CombiSGppConverter::SGppToFullGrid( gridstorageSGpp , alpha , fg );
			}
		}

    };
}

#endif /* SERIALCOMBIGRID_HPP_ */
