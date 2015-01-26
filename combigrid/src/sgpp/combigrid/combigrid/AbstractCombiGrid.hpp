/* ****************************************************************************
* Copyright (C) 2011 Technische Universitaet Muenchen                         *
* This file is part of the SG++ project. For conditions of distribution and   *
* use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
**************************************************************************** */
// @author Janos Benk (benk@in.tum.de)

#ifndef ABSTRACTCOMBIGRID_HPP_
#define ABSTRACTCOMBIGRID_HPP_

#include <sgpp/combigrid/utils/combigrid_ultils.hpp>
#include <sgpp/combigrid/fullgrid/CombiFullGrid.hpp>
#include <sgpp/combigrid/combigridkernel/CombiGridKernel.hpp>
#include <sgpp/combigrid/combischeme/CombiSchemeBasis.hpp>
#include <sgpp/combigrid/domain/CombiGridDomain.hpp>

// ------ SGpp includes -------------
#include <sgpp/base/datatypes/DataVector.hpp>
#include <sgpp/base/grid/GridStorage.hpp>

namespace combigrid {

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
                        const std::vector<bool>& hasBoundaryPts ) : gridDomain_(0) {
        // create the combi kernel which has the non initialized
        combischeme_ = new CombiSchemeBasis(*combischeme);
        combikernel_ = new CombiGridKernelD( combischeme_ , hasBoundaryPts );
      }

      /** Ctor
       * @param combischeme combi schme of the combi grid
       * @param hasBoundaryPts array of flag to indicate weather we have boundary points in the dimension*/
      AbstractCombiGrid(const  CombiSchemeBasis* combischeme ,
                        bool hasBoundaryPts = true ) : gridDomain_(0) {
        // create the combi kernel which has the non initialized
        std::vector<bool> boundaryTmp( combischeme->getDim() , hasBoundaryPts);
        combischeme_ = new CombiSchemeBasis(*combischeme);
        combikernel_ = new CombiGridKernelD( combischeme_ , boundaryTmp );
      }

      /** Dtor which deletes the kernel */
      ~AbstractCombiGrid() {
        // delete the kernel, and with it all the full grids
        delete combikernel_;
      }

      /** sets the domain of the combi grid, this sets globaly for the combi grid,
       *  and not for each full grid */
      void setDomain( GridDomain* gridDomain ) const {
        gridDomain_ = gridDomain;
      }

      /** sets the domain of all the full grids, this is the correct way for extrapolation */
      virtual void setDomainAllFG( GridDomain* gridDomain ) const = 0;

      /** returns the domain of the combi grid */
      const GridDomain* getDomain() const {
        return gridDomain_;
      }

      /** returns the array, size of the dimension, for each dimensions indicates weather there are
       * boundary points */
      const std::vector<bool>& getBoundaryFlags() const {
        return combikernel_->getBoundaryFlags();
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
       * @param coords the coordinates on the unit cube [0,1]^D
       * */
      virtual double eval( std::vector<double>& coords ) const = 0;

      /** evaluate the combi grid at one specified point. Buffered evaluation
       * which might be faster than one evaluation point.
       * @param coords the coordinates on the unit form [0,1]^D
       * @param results the result vector
       */
      virtual void eval( std::vector< std::vector<double> >& coords , std::vector<double>& results ) const = 0;

      /** create SGpp grid storage out of the combi grid <br>
       * @return the created grid storage for the SGpp grid
       *
       * */
      virtual SGPP::base::GridStorage* createSGppGridStorage() const = 0;

      /** in this step the full grid values will be projected to the sparse grid space
       * @param gridstorageSGpp SGpp grid storage
       * @param alpha the coefficient vector which will be overwritten
       * @param minAlpha the minimum of all points, if this argument is set
       * @param maxAlpha the maximum of all points, if this argument is set
       *
       * */
      virtual void reCompose(SGPP::base::GridStorage* gridstorageSGpp , SGPP::base::DataVector* alpha,
                             SGPP::base::DataVector* minAlpha = NULL , SGPP::base::DataVector* maxAlpha = NULL) const = 0;

      /** takes the SGpp sparse grid space and projects to the combi grid
       * @param gridstorageSGpp SGpp grid storage
       * @param alpha the coefficient vector with which the combi grid values will be set
       */
      virtual void deCompose(SGPP::base::GridStorage* gridstorageSGpp , SGPP::base::DataVector* alpha) = 0;

      /** return the combi scheme */
      inline const CombiSchemeBasis* getCombiScheme() const {
        return combischeme_;
      }

      /** return the combi kernel, needed for the combination of the full grids */
      inline const CombiGridKernelD* getCombiKernel() const {
        return combikernel_;
      }
    protected:



      /** pointer to the combi scheme which is set from the constructor argument */
      CombiSchemeBasis* combischeme_;

      /** the combi kernel which stores the full grids */
      CombiGridKernelD* combikernel_;

      /** grid domain for the combi grid */
      mutable GridDomain* gridDomain_;

  };
}


#endif /* ABSTRACTCOMBIGRID_HPP_ */
