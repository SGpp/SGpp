/* ****************************************************************************
* Copyright (C) 2011 Technische Universitaet Muenchen                         *
* This file is part of the SG++ project. For conditions of distribution and   *
* use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
**************************************************************************** */
// @author Janos Benk (benk@in.tum.de)


#ifndef COMBIFULLGRID_HPP_
#define COMBIFULLGRID_HPP_

#include "combigrid/utils/combigrid_ultils.hpp"
#include "combigrid/basisfunction/CombiBasisFunctionBasis.hpp"
#include "combigrid/basisfunction/CombiLinearBasisFunction.hpp"
#include "combigrid/domain/CombiGridDomain.hpp"

namespace combigrid {

  /** The full grid class which is the main building block of the combi grid <br>
   *  It is important that the grid will actually occupy memory when the createFullGrid()
   *  method is called. <br>
   *  The index of a gridpoint in the full grid is given by the formula : <br>
   *  ind = i0 + i1*N0 + i2*N0*N1 + ... + id*N0*N1*N2*...*Nd, where i0,i1,i2,... are the indexes in every dimension,
   *  and Nk is the number of gridpoints in direction k. <br>
   *  For more infos you can look at the "getVectorIndex" and "getLinearIndex" functions and their implementation. <br>
   *  Nk=2^level[k]+1 for every k for the directions with boundary points and Nk=2^level[k]-1 for the directions without boundary. <br>
   *  <br>
   *  The full grid can also be scaled which is done with a separate "combigrid::GridDomain" object. <br>
   * */
  template<typename FG_ELEMENT>
  class FullGrid {

    public:

      /** simplest Ctor with homogeneous  levels */
      FullGrid(int dim, int level, bool hasBdrPoints = true, const BasisFunctionBasis* basis = NULL) {
        // set the basis function for the full grid
        if (basis == NULL) basis_ = LinearBasisFunction::getDefaultBasis();
        else basis_ = basis;

        sgppIndex_ = new std::vector<int>(0);
        gridDomain_ = NULL;
        isFGcreated_ = false;
        dim_ = dim;
        levels_.resize( dim , level );
        hasBoundaryPoints_.resize( dim , hasBdrPoints);
        nrElements_ = 1;
        offsets_.resize(dim_);
        nrPoints_.resize(dim_);

        for (int j = 0 ; j < dim_ ; j++ ) {
          nrPoints_[j] = ( (hasBoundaryPoints_[j] == true) ? (powerOfTwo[levels_[j]] + 1) : (powerOfTwo[levels_[j]] - 1) );
          offsets_[j] = nrElements_;
          nrElements_ = nrElements_ * nrPoints_[j];
        }
      }

      /** dimension adaptive Ctor */
      FullGrid( int dim , const std::vector<int>& levels, bool hasBdrPoints = true, const BasisFunctionBasis* basis = NULL) {
        // set the basis function for the full grid
        if (basis == NULL) basis_ = LinearBasisFunction::getDefaultBasis();
        else basis_ = basis;

        sgppIndex_ = new std::vector<int>(0);
        gridDomain_ = NULL;
        dim_ = dim;
        isFGcreated_ = false;
        levels_ = levels;
        hasBoundaryPoints_.resize( dim , hasBdrPoints);
        nrElements_ = 1;
        offsets_.resize(dim_);
        nrPoints_.resize(dim_);

        for (int j = 0 ; j < dim_ ; j++ ) {
          nrPoints_[j] = ( (hasBoundaryPoints_[j] == true) ? (powerOfTwo[levels_[j]] + 1) : (powerOfTwo[levels_[j]] - 1) );
          offsets_[j] = nrElements_;
          nrElements_ = nrElements_ * nrPoints_[j];
        }
      }

      /** dimension adaptive Ctor */
      FullGrid( int dim , const std::vector<int>& levels, const std::vector<bool>& hasBdrPoints, const BasisFunctionBasis* basis = NULL) {
        // set the basis function for the full grid
        if (basis == NULL) basis_ = LinearBasisFunction::getDefaultBasis();
        else basis_ = basis;

        sgppIndex_ = new std::vector<int>(0);
        gridDomain_ = NULL;
        dim_ = dim;
        isFGcreated_ = false;
        levels_ = levels;
        hasBoundaryPoints_ = hasBdrPoints;
        nrElements_ = 1;
        offsets_.resize(dim_);
        nrPoints_.resize(dim_);

        for (int j = 0 ; j < dim_ ; j++ ) {
          nrPoints_[j] = ( (hasBoundaryPoints_[j] == true) ? (powerOfTwo[levels_[j]] + 1) : (powerOfTwo[levels_[j]] - 1) );
          offsets_[j] = nrElements_;
          nrElements_ = nrElements_ * nrPoints_[j];
        }
      }

      virtual ~FullGrid() {
        delete sgppIndex_;
      }

      /** allocates the memory for the element vector of the full grid <br>
       * Only after this needs the full grid considerable amount of memory */
      void createFullGrid() {
        isFGcreated_ = true;
        fullgridVector_.resize(nrElements_);
      }

      /** Deallocates the element vector, so the memory requirement should not be significant afer this*/
      void deleteFullGrid() {
        isFGcreated_ = false;
        //    fullgridVector_.flush();
        fullgridVector_.clear();
        delete sgppIndex_;
        sgppIndex_ = new std::vector<int>(0);
      }

      /** evaluates the full grid on the specified coordinates
       * @param coords ND coordinates on the unit square [0,1]^D*/
      FG_ELEMENT eval(std::vector<double>& coords) const {
        int ii, i, tmp_val, vv, nr;
        int jj;
        double baseVal;
        double normcoord;
        // this value will be reseted, but just to avoid compiler warnings
        FG_ELEMENT ret_val = fullgridVector_[0];
        double intersect[127]; // alternatievly [127]
        int aindex[63]; // alternatievly [63]
        //int verb = 6;

        // if there is a transformation then transform to the unit coordinates
        if (gridDomain_ != NULL ) {
          for (ii = dim_ - 1 ; ii >= 0; ii--) {
            (gridDomain_->get1DDomain(ii)).transformRealToUnit( coords[ii] , coords[ii] , levels_[ii] , (hasBoundaryPoints_[ii] == false));
          }
        }

        // the coordinates are on the unit square
        for ( ii = dim_ - 1 ; ii >= 0; ii--) {
          // the value has to be between 0.0 and 1.0
          //coords[ii] = (coords[ii] > 1.0) ? 1.0 : coords[ii];
          //coords[ii] = (coords[ii] < 0.0) ? 0.0 : coords[ii];
          //if ((coords[ii]<0.0)||(coords[ii]>1.0)) return 0.0;
          // scale to the reference [0,1] interval the intersection point
          normcoord = coords[ii] * powerOfTwo[levels_[ii]];
          aindex[ii] = static_cast<int>( std::floor(normcoord) );

          if (hasBoundaryPoints_[ii] == true ) {
            aindex[ii] = (aindex[ii] < 0) ? 0 : aindex[ii];
            aindex[ii] = (aindex[ii] >= (int)(nrPoints_[ii] - 1)) ? (int)(nrPoints_[ii] - 2) : aindex[ii];
            //calculate the coordinates
            normcoord =  normcoord - (double)aindex[ii];
          } else {
            aindex[ii] = (aindex[ii] < 0) ? 0 : aindex[ii];
            aindex[ii] = (aindex[ii] >= (int)(nrPoints_[ii])) ? (int)(nrPoints_[ii]) : aindex[ii];
            normcoord =  normcoord - (double)aindex[ii];

            if (aindex[ii] <= 0) {
              aindex[ii] = 1;
              // extrapolation to the left, must be a negative number
              normcoord = coords[ii] * powerOfTwo[levels_[ii]] - 1;
              //COMBIGRID_OUT_LEVEL3( verb , "FullG L , aindex[ii]:" << aindex[ii]);
            } else {
              // make extrapolation at the last cell (where there is no boundary point)
              if (aindex[ii] >= (int)nrPoints_[ii] - 1) {
                aindex[ii] = (int)nrPoints_[ii] - 1;
                // extrapolation to the right this will be a >=1 number
                normcoord = coords[ii] * powerOfTwo[levels_[ii]] - (double)aindex[ii]; // this should be 1 + ...
                //COMBIGRID_OUT_LEVEL3( verb , "FullG R , aindex[ii]:" << aindex[ii]);
              }
            }

            // we have to subtract one in case of no boundary points
            aindex[ii] = aindex[ii] - 1;
            // in some special cases this might become negative, so limit to zero (when level = 1 and no boundary point )
            aindex[ii] = (aindex[ii] < 0) ? 0 : aindex[ii];
          }

          //COMBIGRID_OUT_LEVEL3( verb , "FullGrid::eval ii:" << ii << " , coords[ii]:" << coords[ii] << " , aindex[ii]:" << aindex[ii]
          //                     << " , level[ii]:" << levels_[ii] << " , normcoord:" << normcoord );
          // evaluate the basis functions on the 1D coordinates
          intersect[2 * ii] = basis_->functionEval1( normcoord );
          intersect[2 * ii + 1] = basis_->functionEval2( normcoord );
        }

        nr = powerOfTwo[dim_];

        for (ii = 0; ii < nr; ii++) {
          baseVal = 1;//we will store here the coefficient for the next corner point of the cuboid
          i = 0;
          tmp_val = ii;
          vv = 0;

          // compute the "dim" dimensional basis function value(for one node) and the corresponding vector index
          for (jj = dim_ - 1 ; jj >= 0; jj--) {
            vv = ( nrPoints_[jj] > 1) ? (tmp_val & 1) : 0; // if we have only one point then use the same
            baseVal = baseVal * intersect[2 * jj + vv];
            i = i * nrPoints_[jj] + aindex[jj] + vv;
            tmp_val = tmp_val >> 1;
          }

          // if this is the first point then set directly the value, otherwise add the value
          //COMBIGRID_OUT_LEVEL3( verb , "Vect elem i:" << i << " , vect:" << fullgridVector_[i]);
          ret_val = (ii == 0) ? (baseVal * fullgridVector_[i]) : (ret_val + baseVal * fullgridVector_[i]);
        }

        //COMBIGRID_OUT_LEVEL3( verb , "FullGrid::eval RET VAL:" << ret_val);
        // just return the value
        return ret_val;
      }

      /** return the coordinates on the unit square
       * @param elemIndex [IN] index of the element in the vector
       * @param coords [OUT] the vector must be resized already*/
      void getCoords(int elemIndex , std::vector<double>& coords) const {
        // temporary variables
        //int verb = 6;
        int ind = 0;
        int tmp_add = 0;

        if (gridDomain_ == NULL) {
          for (int j = 0 ; j < dim_ ; j++) {
            ind       = elemIndex % (int)(nrPoints_[j]);
            elemIndex = elemIndex / (int)(nrPoints_[j]);
            // set the coordinate based on if we have boundary points
            tmp_add = (hasBoundaryPoints_[j] == true) ? (0) : (1);
            coords[j] = ((double)(ind + tmp_add)) * oneOverPowOfTwo[levels_[j]];
            //COMBIGRID_OUT_LEVEL3( verb , "FullGrid::getCoords j:" << j << " , coords[j]:" << coords[j]);
          }
        } else {
          // we have a valid Domain
          for (int j = 0 ; j < dim_ ; j++) {
            ind       = elemIndex % (int)(nrPoints_[j]);
            elemIndex = elemIndex / (int)(nrPoints_[j]);
            // set the coordinate based on if we have boundary points
            tmp_add = (hasBoundaryPoints_[j] == true) ? (0) : (1);
            (gridDomain_->get1DDomain(j)).transformUnitToReal( levels_[j] , ind + tmp_add , coords[j] );
            //COMBIGRID_OUT_LEVEL3( verb , "FullGrid::getCoords j:" << j << " , coords[j]:" << coords[j]);
          }
        }
      }


      /* returns the LI (level,index) notation for a given element in the full grid
       * @param elementIndex [IN] the linear index of the element
       * @param levels [OUT] the levels of the point in the LI notation
       * @param indexes [OUT] the indexes of the point in the LI notation */
      //void getLI_std(int elementIndex , std::vector<int>& levels , std::vector<int>& indexes) const {
      //  getLI( elementIndex , &(levels[0]) , &(indexes[0]) );
      //}

      /** returns the LI (level,index) notation for a given element in the full grid
       * @param elementIndex [IN] the linear index of the element
       * @param levels [OUT] the levels of the point in the LI notation
       * @param indexes [OUT] the indexes of the point in the LI notation */
      void getLI(int elementIndex ,  std::vector<int>& levels , std::vector<int>& indexes ) const {
        int k , startindex , tmp_val ;

        tmp_val = elementIndex;

        // first calculate intermediary indexes
        for ( k = 0 ; k < dim_ ; k++ ) {
          startindex = (hasBoundaryPoints_[k]) ? 0 : 1;
          indexes[k] = tmp_val % nrPoints_[k] + startindex;
          tmp_val = tmp_val / nrPoints_[k];
        }

        //The level and index of the element in the hashgridstorage are computed dividing by two the index and level in the fullgrid
        //until we obtain an impair number for the index, thus obtaining the level and index in the hierarchical basis (Aliz Nagy)
        // ...
        for ( k = 0 ; k < dim_ ; k++) {
          tmp_val = levels_[k];

          if ( indexes[k] != 0 ) {
            // todo: these operations can be optimized
            while ( indexes[k] % 2 == 0 ) {
              indexes[k] = indexes[k] / 2;
              tmp_val--;
            }
          } else {
            tmp_val = 0;
          }

          levels[k] = tmp_val;
        }
      }

      /** returns the vector index for one linear index
       * @param linIndex [IN] the linear index
       * @param axisIndex [OUT] the returned vector index */
      inline void getVectorIndex(const int linIndex, std::vector<int>& axisIndex) const {
        int tmp = linIndex;

        for (int i = dim_ - 1 ; i >= 0 ; i--) {
          axisIndex[i] = tmp / (this->getOffset(i));
          tmp = tmp % this->getOffset(i);
        }
      }

      /** returns the linear index for one linear index
       * @param axisIndex [IN] the vector index */
      inline int getLinearIndex(const std::vector<int>& axisIndex) const {
        int tmp = 0;

        for (int i = dim_ - 1 ; i >= 0 ; i--) {
          tmp = tmp +  offsets_[i] * axisIndex[i];
        }

        return tmp;
      }

      /** sets the domain of the full grid */
      inline void setDomain( GridDomain* gridDomain ) const {
        gridDomain_ = gridDomain;
      }

      /** returns the domain of the full grid */
      inline const GridDomain* getDomain() const {
        return gridDomain_;
      }

      /** returns the domain of the full grid */
      inline GridDomain* getDomain() {
        return gridDomain_;
      }

      /** */
      void addCustomDomain(std::vector<double> min, std::vector<double> max) {
        gridDomain_ = new GridDomain(dim_, min, max);
      }

      /** returns pointer to the basis function */
      inline const BasisFunctionBasis* getBasisFct() const {
        return basis_;
      }

      /** returns the dimension of the full grid */
      inline int getDimension() const {
        return dim_;
      }

      /** the getters for the full grid vector */
      inline std::vector<FG_ELEMENT>& getElementVector() {
        return fullgridVector_;
      }

      inline const std::vector<FG_ELEMENT>& getElementVector() const {
        return fullgridVector_;
      }

      /** return the offset in the full grid vector of the dimension */
      inline int getOffset(int i) const {
        return offsets_[i];
      }
      inline const std::vector<int>& getOffsets() const {
        return offsets_;
      }

      /** return the level vector */
      inline const std::vector<int>& getLevels() const {
        return levels_;
      }
      inline std::vector<int>& getLevels() {
        return levels_;
      }

      /** flag shows if the full grid is created */
      inline bool isGridCreated() const {
        return isFGcreated_;
      }

      /** returns the number of elements in the full grid */
      inline int getNrElements() const {
        return nrElements_;
      }

      /** number of points per dimension i */
      inline int length(int i) const {
        return nrPoints_[i];
      }

      /** return the vector for faster combination of the full grids <br>
       * this will be used in the Converter */
      inline std::vector<int>& getSGppIndex() const {
        return (*sgppIndex_);
      }

      /** vector of flags to show if the dimension has boundary points*/
      inline const std::vector<bool>& returnBoundaryFlags() const {
        return hasBoundaryPoints_;
      }

      /** copies the input vector to the full grid vector
       * @param in [IN] input vector*/
      void setElementVector(std::vector<FG_ELEMENT> in) {
        fullgridVector_ = in;
      }

    private:


      /** dimension of the full grid */
      int dim_;

      /** the size of the vector, nr of total elements */
      int nrElements_;

      /** flag to show if the grid is created */
      bool isFGcreated_;

      /** levels for each dimension */
      std::vector<int> levels_;

      /** number of points per axis*/
      std::vector<int> nrPoints_;

      /** flag to show if the dimension has boundary points*/
      std::vector<bool> hasBoundaryPoints_;

      /** the offsets in each direction*/
      std::vector<int> offsets_;

      /** the full grid vector, this contains the elements of the full grid */
      std::vector<FG_ELEMENT> fullgridVector_;

      /** pointer to the function basis*/
      const BasisFunctionBasis* basis_;

      /** the domain transformation */
      mutable GridDomain* gridDomain_;

      /** vector for the SGpp index, for recombination*/
      mutable std::vector<int>* sgppIndex_;

  };

}

#endif /* COMBIFULLGRID_HPP_ */
