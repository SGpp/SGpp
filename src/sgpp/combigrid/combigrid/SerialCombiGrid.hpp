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
using namespace sg::base;

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
		virtual void createFullGrids() ;

		/** see supercalss for docu */
		virtual FullGridD* getFullGrid( int i ) { return combikernel_->getFullGrid(i); }

		/** see supercalss for docu */
		virtual const FullGridD* getFullGrid( int i ) const { return combikernel_->getFullGrid(i); }

		/** see supercalss for docu */
		virtual int getNrFullGrid() const { return combikernel_->getNrFullGrids(); }

		/** see supercalss for docu */
		virtual double eval( const std::vector<double>& coords ) const ;

		/** see supercalss for docu */
		virtual void eval( const std::vector< std::vector<double> >& coords , std::vector<double>& results ) const ;

		/** see supercalss for docu */
		virtual GridStorage* createSGppGridStorage() const ;

		/** see supercalss for docu */
		virtual void reCompose(GridStorage* gridstorageSGpp , DataVector* alpha,
				DataVector* minAlpha = NULL , DataVector* maxAlpha = NULL) const ;

		/** see supercalss for docu */
		virtual void deCompose(GridStorage* gridstorageSGpp , DataVector* alpha) ;

    };
}

#endif /* SERIALCOMBIGRID_HPP_ */
