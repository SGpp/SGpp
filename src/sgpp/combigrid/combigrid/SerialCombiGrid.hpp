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

		/** sets the domain for all fullgrids */
		virtual void setDomainAllFG( GridDomain* gridDomain ) const;

		/** see superclass for docu */
		virtual FullGridD* getFullGrid( int i ) { return combikernel_->getFullGrid(i); }

		/** see superclass for docu */
		virtual const FullGridD* getFullGrid( int i ) const { return combikernel_->getFullGrid(i); }

		/** see superclass for docu */
		virtual int getNrFullGrid() const { return combikernel_->getNrFullGrids(); }

		double evalSingleGrid(int index,std::vector<double> & coords) const;

		/** see superclass for docu */
		virtual double eval( std::vector<double>& coords ) const ;

		/** see superclass for docu */
		virtual void eval( std::vector< std::vector<double> >& coords , std::vector<double>& results ) const ;

		/** see superclass for docu */
		virtual sg::base::GridStorage* createSGppGridStorage() const ;

		/** see superclass for docu */
		virtual void reCompose(sg::base::GridStorage* gridstorageSGpp , sg::base::DataVector* alpha,
				sg::base::DataVector* minAlpha = NULL , sg::base::DataVector* maxAlpha = NULL) const ;

		/** see superclass for docu */
		virtual void deCompose(sg::base::GridStorage* gridstorageSGpp , sg::base::DataVector* alpha) ;
    };
}

#endif /* SERIALCOMBIGRID_HPP_ */
