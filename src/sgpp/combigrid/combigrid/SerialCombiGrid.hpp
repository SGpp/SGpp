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

		/** see superclass for docu */
		virtual void createFullGrids() ;

		/** see superclass for docu */
		virtual FullGridD* getFullGrid( int i ) { return combikernel_->getFullGrid(i); }

		/** see superclass for docu */
		virtual const FullGridD* getFullGrid( int i ) const { return combikernel_->getFullGrid(i); }

		/** see superclass for docu */
		virtual int getNrFullGrid() const { return combikernel_->getNrFullGrids(); }

		/** see superclass for docu */
		virtual double eval( const std::vector<double>& coords ) const ;

		/** see superclass for docu */
		virtual void eval( const std::vector< std::vector<double> >& coords , std::vector<double>& results ) const ;

		/** see superclass for docu */
		virtual sg::GridStorage* createSGppGridStorage() const ;

		/** see superclass for docu */
		virtual void reCompose(sg::GridStorage* gridstorageSGpp , DataVector* alpha,
				DataVector* minAlpha = NULL , DataVector* maxAlpha = NULL) const ;

		/** see superclass for docu */
		virtual void deCompose(sg::GridStorage* gridstorageSGpp , DataVector* alpha) ;

    };
}

#endif /* SERIALCOMBIGRID_HPP_ */
