/*
 * AbstractCombiGrid.hpp
 *
 *  Created on: Feb 23, 2011
 *      Author: benk
 */

#ifndef ABSTRACTCOMBIGRID_HPP_
#define ABSTRACTCOMBIGRID_HPP_

#include "combigrid/utils/combigrid_ultils.hpp"
#include "combigrid/fullgrid/CombiFullGrid.hpp"
#include "combigrid/combigridkernel/CombiGridKernel.hpp"
#include "combigrid/combischeme/CombiSchemeBasis.hpp"

// ------ SGpp includes -------------
#include "data/DataVector.hpp"
#include "grid/GridStorage.hpp"

using namespace std;

namespace combigrid{

	/** full grid type definition */
    typedef FullGrid< double > FullGridD;

    /** combi kernel type definition */
    typedef CombiGridKernel< double > CombiGridKernelD;

	/** The virtual class which defines the interface of the combi grid. <br>
	 * This interface can be implemented by several implementations, serial, OpenMP or even MPI */
	class AbstractCombiGrid {

	public:

    	/** Ctor
    	 * @param combischeme combi schme of the combi grid
    	 * @param hasBoundaryPts array of flag to indicate weather we have boundary points in the dimension*/
		AbstractCombiGrid(const CombiSchemeBasis* combischeme ,
				const std::vector<bool>& hasBoundaryPts ) :combischeme_(combischeme){
			// create the combi kernel which has the non initialized
			combikernel_ = new CombiGridKernelD( combischeme , hasBoundaryPts );
		}

    	/** Ctor
    	 * @param combischeme combi schme of the combi grid
    	 * @param hasBoundaryPts array of flag to indicate weather we have boundary points in the dimension*/
		AbstractCombiGrid(const CombiSchemeBasis* combischeme ,
				bool hasBoundaryPts = true ) :combischeme_(combischeme){
			// create the combi kernel which has the non initialized
			std::vector<bool> boundaryTmp( combischeme->getDim() , hasBoundaryPts);
			combikernel_ = new CombiGridKernelD( combischeme , boundaryTmp );
		}

		/** Dtor which deletes the kernel */
		~AbstractCombiGrid(){
			// delete the kernel, and with it all the full grids
			delete combikernel_;
		}

		/** create the actual vector for the full grids. <br>
		 * This is be different for serial and parallel implementations */
		virtual void createFullGrids() = 0;

		/** returns one full grid . <br>
		 * In case of the serial implementation this will be the whole
		 * combination scheme, but in the
		 * @param i is the local index of the full grid */
		virtual FullGridD* getFullGrid( int i ) = 0;

		/** same as the previous method but with const environment */
		virtual const FullGridD* getFullGrid( int i ) const = 0;

		/** get the number of full grids. <br>
		 * In serial case this is simple but in parallel (MPI) case this might be only a part
		 * of the combi scheme */
		virtual int getNrFullGrid() const = 0;

		/** evaluate the combi grid at one specified point
		 * @param coords , the coordinates on the unit form [0,1]^D */
		virtual double eval( const std::vector<double>& coords ) const = 0;

		/** evaluate the combi grid at one specified point. Buffered evaluation
		 * which might be faster than one evaluation point.
		 * @param coords , the coordinates on the unit form [0,1]^D
		 * @param results , the result vector */
		virtual void eval( const std::vector< std::vector<double> >& coords , std::vector<double>& results ) const = 0;

		/** create SGpp grid storage out of the combi grid <br>
		 * @return the created grid storage for the SGpp grid */
		virtual sg::GridStorage* createSGppGridStorage() const = 0;

		/** in this step the full grid values will be projected to the sparse grid space
		 * @param gridstorageSGpp  , SGpp grid storage
		 * @param alpha , the coefficient vector which will be overwritten
		 * @param minAlpha , the minimum of all points, if this argument is set
		 * @param maxAlpha , the maximum of all points, if this argument is set*/
		virtual void reCompose(sg::GridStorage* gridstorageSGpp , DataVector* alpha,
				DataVector* minAlpha = NULL , DataVector* maxAlpha = NULL) const = 0;

		/** takes the SGpp sparse grid space and projects to the combi grid
		 * @param gridstorageSGpp , SGpp grid storage
		 * @param alpha , the coefficient vector with which the combi grid values will be set */
		virtual void deCompose(sg::GridStorage* gridstorageSGpp , DataVector* alpha) = 0;

		/** return the combi kernel, needed for the combination of the full grids */
		inline const CombiGridKernelD* getCombiKernel() const { return combikernel_; }

		/** return the combi scheme */
		inline const CombiSchemeBasis* getCombiScheme() const { return combischeme_; }

protected:

		/** pointer to the combi scheme which is set from the constructor argument */
		const CombiSchemeBasis* combischeme_;

		/** the combi kernel which stores the full grids */
		CombiGridKernelD* combikernel_;

	};
}


#endif /* ABSTRACTCOMBIGRID_HPP_ */
