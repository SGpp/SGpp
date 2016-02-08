/* ****************************************************************************
 * Copyright (C) 2011 Technische Universitaet Muenchen                         *
 * This file is part of the SG++ project. For conditions of distribution and   *
 * use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
 **************************************************************************** */
// @author Janos Benk (benk@in.tum.de)
#ifndef COMBIFULLGRID_HPP_
#define COMBIFULLGRID_HPP_

#include <boost/serialization/access.hpp>
#include "boost/serialization/vector.hpp"
#include "boost/serialization/complex.hpp"
#include "boost/archive/binary_oarchive.hpp"
#include "boost/archive/binary_iarchive.hpp"
#include "boost/archive/text_oarchive.hpp"
#include "boost/archive/text_iarchive.hpp"
#include <string>
#include <assert.h>
#include "sgpp/combigrid/utils/combigrid_ultils.hpp"
#include "sgpp/distributedcombigrid/utils/Types.hpp"
#include "sgpp/combigrid/basisfunction/CombiBasisFunctionBasis.hpp"
#include "sgpp/combigrid/basisfunction/CombiLinearBasisFunction.hpp"
#include "sgpp/combigrid/domain/CombiGridDomain.hpp"
#include "sgpp/distributedcombigrid/utils/LevelVector.hpp"

// switch on alternative assignment of level vector: the boundary points have
// level 1 and not level 0
//#define ALT_LEVEL_VECTOR

namespace combigrid {

template<typename FG_ELEMENT>
class SGrid;

// forward declarations, we cannot include Hierarchization.hpp
// before FullGrid.hpp is completely parsed
class Hierarchization;

// forward declarations necessary here because of SGrid and HierarchizedFullGrid
// including each other
// todo: separate declaration and definition. besides better readability, this
// problem would have been avoided. it cost me several hours to figure that out

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
  FullGrid(DimType dim, LevelType level, bool hasBdrPoints = true,
           const BasisFunctionBasis* basis = NULL);

  /** dimension adaptive Ctor */
  FullGrid(DimType dim, const LevelVector& levels, bool hasBdrPoints = true,
           const BasisFunctionBasis* basis = NULL);

  /** dimension adaptive Ctor */
  FullGrid(DimType dim, const LevelVector& levels,
           const std::vector<bool>& hasBdrPoints, const BasisFunctionBasis* basis =
             NULL);

  // load archived fg from file
  FullGrid(const char* filename, const BasisFunctionBasis* basis = NULL);

  /* create hierarchized fullgrid from SGrid */
  FullGrid(const LevelVector& levels, const SGrid<FG_ELEMENT>& sg);

  virtual
  ~FullGrid();

  /** allocates the memory for the element vector of the full grid <br>
   * Only after this needs the full grid considerable amount of memory */
  void
  createFullGrid();

  /** Deallocates the element vector, so the memory requirement should not be significant afer this*/
  void
  deleteFullGrid();

  /** evaluates the full grid on the specified coordinates
   * @param coords ND coordinates on the unit square [0,1]^D*/
  FG_ELEMENT
  eval(std::vector<real>& coords);

  /** return the coordinates on the unit square
   * @param elemIndex [IN] index of the element in the vector
   * @param coords [OUT] the vector must be resized already*/
  void
  getCoords(IndexType elemIndex, std::vector<real>& coords) const;

  /** returns the LI (level,index) notation for a given element in the full grid
   * @param elementIndex [IN] the linear index of the element
   * @param levels [OUT] the levels of the point in the LI notation
   * @param indexes [OUT] the indexes of the point in the LI notation */
  void
  getLI(IndexType elementIndex, LevelVector& levels,
        IndexVector& indexes) const;

  /** returns the vector index for one linear index
   * @param linIndex [IN] the linear index
   * @param axisIndex [OUT] the returned vector index */
  inline void
  getVectorIndex(const IndexType linIndex, IndexVector& axisIndex) const;

  /** returns the linear index for one linear index
   * @param axisIndex [IN] the vector index */
  inline IndexType
  getLinearIndex(const IndexVector& axisIndex) const;

  /** sets the domain of the full grid */
  inline void
  setDomain(GridDomain* gridDomain) const;

  /** returns the domain of the full grid */
  inline const GridDomain*
  getDomain() const;

  /** returns the domain of the full grid */
  inline GridDomain*
  getDomain();

  /** returns pointer to the basis function */
  inline const BasisFunctionBasis*
  getBasisFct() const;

  /** returns the dimension of the full grid */
  inline DimType
  getDimension() const;

  /** the getters for the full grid vector */
  inline std::vector<FG_ELEMENT>&
  getElementVector();

  inline const std::vector<FG_ELEMENT>&
  getElementVector() const;

  /** return the offset in the full grid vector of the dimension */
  inline IndexType
  getOffset(DimType i) const;

  inline const IndexVector&
  getOffsets() const;

  /** return the level vector */
  inline const LevelVector&
  getLevels() const;

  /** flag shows if the full grid is created */
  inline bool
  isGridCreated() const;

  /** function to check whether a fullgrid is hierarchized */
  inline bool
  isHierarchized() const;

  /** set the isHierarchized flag */
  inline void
  setHierarchized();

  /** returns the number of elements in the full grid */
  inline IndexType
  getNrElements() const;

  /** returns the number of elements in the full grid */
  inline IndexType
  getNrElementsNoBoundary() const;

  /** number of points per dimension i */
  inline IndexType
  length(DimType i) const;

  /** return the vector for faster combination of the full grids <br>
   * this will be used in the Converter */
  inline IndexVector&
  getSGppIndex() const;

  /** vector of flags to show if the dimension has boundary points*/
  inline const std::vector<bool>&
  returnBoundaryFlags() const;

  /** copies the input vector to the full grid vector
   * in most cases the add function would be more save
   * @param in [IN] input vector*/
  void
  setElementVector(const std::vector<FG_ELEMENT>& in);

  inline FG_ELEMENT*
  getData();

  inline const FG_ELEMENT*
  getData() const;

  void
  gridEval(FullGrid<FG_ELEMENT>& dst) const;

  void
  save(std::string& filename);

  /* add another full grid to the current full grid
   *
   * if the grids differ in discretization or boundary flags interpolation
   * based on the eval function is used
   */
  void
  add(FullGrid<FG_ELEMENT>& fg, real coeff);

  /* add another full grid to the current full grid
   *
   * special treatment for GENE grids which have fourier coefficients
   * in x-direction: no interpolation happens in x-direction. by using
   * a coordinate transformation we make sure the grids fit toghether in
   * the desired way
   */
  void
  addGENE(const FullGrid<FG_ELEMENT>& fg, real coeff);

  /* get lp norm of the grid. for p = 0 maximum norm is returned. */
  inline real
  getlpNorm(int p);

  /* normalize with lp norm and return norm */
  inline real
  normalizelp(int p);

  // return extents of grid
  inline const IndexVector& getSizes() const;

  inline void print(std::ostream& os) const;

 private:
  /** dimension of the full grid */
  DimType dim_;

  /** the size of the vector, nr of total elements */
  IndexType nrElements_;

  /** ADDED the size of the vector, nr of total elements without boundary */
  IndexType nrElementsNoBoundary_;

  /** flag to show if the grid is created */
  bool isFGcreated_;

  /** levels for each dimension */
  LevelVector levels_;

  /** number of points per axis*/
  IndexVector nrPoints_;

  /** flag to show if the dimension has boundary points*/
  std::vector<bool> hasBoundaryPoints_;

  /** the offsets in each direction*/
  IndexVector offsets_;

  /** the full grid vector, this contains the elements of the full grid */
  std::vector<FG_ELEMENT> fullgridVector_;

  /** pointer to the function basis*/
  const BasisFunctionBasis* basis_;

  /** the domain transformation */
  mutable GridDomain* gridDomain_;

  friend class boost::serialization::access;

  // serialize
  template<class Archive>
  void
  serialize(Archive& ar, const unsigned int version);

  // the hierarchization class should be the only class to modify the
  // isHierarchized flag
  friend class Hierarchization;

  /** indicates if the fullgrid is hierarchized, i.e. the stored values
   * are the coefficients of the hierarchical basis functions
   */
  bool isHierarchized_;

  void print1D(std::ostream& os) const;

  void print2D(std::ostream& os) const;

  void print3D(std::ostream& os) const;

};

} // namespace combigrid

namespace combigrid {

template<typename FG_ELEMENT>
FullGrid<FG_ELEMENT>::FullGrid(DimType dim, LevelType level, bool hasBdrPoints,
                               const BasisFunctionBasis* basis) {
  // set the basis function for the full grid
  if (basis == NULL)
    basis_ = LinearBasisFunction::getDefaultBasis();
  else
    basis_ = basis;

  gridDomain_ = NULL;
  isFGcreated_ = false;
  dim_ = dim;
  levels_.resize(dim, level);
  hasBoundaryPoints_.resize(dim, hasBdrPoints);
  nrElements_ = 1;

  // ADDED
  nrElementsNoBoundary_ = 1;

  offsets_.resize(dim_);
  nrPoints_.resize(dim_);

  for (DimType j = 0; j < dim_; j++) {
    nrPoints_[j] = (
                     (hasBoundaryPoints_[j] == true) ?
                     (powerOfTwo[levels_[j]] + 1) : (powerOfTwo[levels_[j]] - 1));
    offsets_[j] = nrElements_;
    nrElements_ = nrElements_ * nrPoints_[j];

    // ADDED
    nrElementsNoBoundary_ = nrElementsNoBoundary_
                            * (powerOfTwo[levels_[j]] - 1);
  }

  isHierarchized_ = false;
}

/** dimension adaptive Ctor */
template<typename FG_ELEMENT>
FullGrid<FG_ELEMENT>::FullGrid(DimType dim, const LevelVector& levels,
                               bool hasBdrPoints, const BasisFunctionBasis* basis) {
  // set the basis function for the full grid
  if (basis == NULL)
    basis_ = LinearBasisFunction::getDefaultBasis();
  else
    basis_ = basis;

  gridDomain_ = NULL;
  dim_ = dim;
  isFGcreated_ = false;
  levels_ = levels;
  hasBoundaryPoints_.resize(dim, hasBdrPoints);
  nrElements_ = 1;
  // ADDED
  nrElementsNoBoundary_ = 1;

  offsets_.resize(dim_);
  nrPoints_.resize(dim_);

  for (DimType j = 0; j < dim_; j++) {
    nrPoints_[j] = (
                     (hasBoundaryPoints_[j] == true) ?
                     (powerOfTwo[levels_[j]] + 1) : (powerOfTwo[levels_[j]] - 1));
    offsets_[j] = nrElements_;
    nrElements_ = nrElements_ * nrPoints_[j];

    // ADDED
    nrElementsNoBoundary_ = nrElementsNoBoundary_
                            * (powerOfTwo[levels_[j]] - 1);
  }

  isHierarchized_ = false;
}

/** dimension adaptive Ctor */
template<typename FG_ELEMENT>
FullGrid<FG_ELEMENT>::FullGrid(DimType dim, const LevelVector& levels,
                               const std::vector<bool>& hasBdrPoints, const BasisFunctionBasis* basis) {
  assert(levels.size() == dim);
  assert(hasBdrPoints.size() == dim);

  // set the basis function for the full grid
  if (basis == NULL)
    basis_ = LinearBasisFunction::getDefaultBasis();
  else
    basis_ = basis;

  gridDomain_ = NULL;
  dim_ = dim;
  isFGcreated_ = false;
  levels_ = levels;
  hasBoundaryPoints_ = hasBdrPoints;
  nrElements_ = 1;
  // ADDED
  nrElementsNoBoundary_ = 1;

  offsets_.resize(dim_);
  nrPoints_.resize(dim_);

  for (DimType j = 0; j < dim_; j++) {
    nrPoints_[j] = (
                     (hasBoundaryPoints_[j] == true) ?
                     (powerOfTwo[levels_[j]] + 1) : (powerOfTwo[levels_[j]] - 1));
    offsets_[j] = nrElements_;
    nrElements_ = nrElements_ * nrPoints_[j];
    nrElementsNoBoundary_ = nrElementsNoBoundary_
                            * (powerOfTwo[levels_[j]] - 1);
  }

  isHierarchized_ = false;
}

template<typename FG_ELEMENT>
FullGrid<FG_ELEMENT>::FullGrid(const char* filename,
                               const BasisFunctionBasis* basis) {
  // set the basis function for the full grid
  if (basis == NULL)
    basis_ = LinearBasisFunction::getDefaultBasis();
  else
    basis_ = basis;

  std::ifstream ifs(filename, std::ios::binary);
  //boost::archive::binary_iarchive ia(ifs);
  boost::archive::text_iarchive ia(ifs);
  ia >> *this;

  gridDomain_ = NULL;

  isHierarchized_ = false;

  isFGcreated_ = false;
}

/* create hierarchized fullgrid from SGrid */
template<typename FG_ELEMENT>
FullGrid<FG_ELEMENT>::FullGrid(const LevelVector& levels,
                               const SGrid<FG_ELEMENT>& sg) {
  SGrid<FG_ELEMENT>& sg2 = const_cast<SGrid<FG_ELEMENT>&>(sg);

  assert(levels.size() == sg2.getDim());

  // check if hfg contained in sg
  // in principle we dont need to do this and only copy the data of subspaces
  // which exists in both spaces
  assert(sg2.isContained(levels));

  // set the basis function for the full grid
  // at the moment we don't have any other basis for sg
  basis_ = LinearBasisFunction::getDefaultBasis();

  gridDomain_ = NULL;
  dim_ = sg2.getDim();
  isFGcreated_ = false;
  levels_ = levels;
  hasBoundaryPoints_ = sg2.getBoundaryVector();
  nrElements_ = 1;
  nrElementsNoBoundary_ = 1;

  offsets_.resize(dim_);
  nrPoints_.resize(dim_);

  for (DimType j = 0; j < dim_; j++) {
    nrPoints_[j] = (
                     (hasBoundaryPoints_[j] == true) ?
                     (powerOfTwo[levels_[j]] + 1) : (powerOfTwo[levels_[j]] - 1));
    offsets_[j] = nrElements_;
    nrElements_ = nrElements_ * nrPoints_[j];
    // ADDED
    nrElementsNoBoundary_ = nrElementsNoBoundary_
                            * (powerOfTwo[levels_[j]] - 1);
  }

  createFullGrid();

  // create iterator for each subspace of sgrid and set to begin
  typedef typename std::vector<FG_ELEMENT>::iterator SubspaceIterator;
  typename std::vector<SubspaceIterator> it_sub(sg2.getSize());

  for (size_t i = 0; i < it_sub.size(); ++i)
    it_sub[i] = sg2.getDataVector(i).begin();

  for (IndexType i = 0; i < getNrElements(); ++i) {
    // get level and index vector of element i
    LevelVector lvec(dim_);
    IndexVector ivec(dim_);
    getLI(i, lvec, ivec);

    // find index of sgrid space
    size_t j = sg2.getLevelIndex(lvec);

    // initialize subspace if necessary
    // this might invalidate the iterator, thus reassign iterator
    if (sg2.getDataSize(j) == 0) {
      sg2.initHierarchicalSpace(j, FG_ELEMENT(0));
      it_sub[j] = sg2.getDataVector(j).begin();
    }

    assert(it_sub[j] != sg2.getDataVector(j).end());

    // copy sg element to hfg
    fullgridVector_[i] = *it_sub[j];

    ++it_sub[j];
  }

  isHierarchized_ = true;
}

template<typename FG_ELEMENT>
FullGrid<FG_ELEMENT>::~FullGrid() {
}

/** allocates the memory for the element vector of the full grid <br>
 * Only after this needs the full grid considerable amount of memory */
template<typename FG_ELEMENT>
void FullGrid<FG_ELEMENT>::createFullGrid() {
  isFGcreated_ = true;
  fullgridVector_.resize(nrElements_, FG_ELEMENT(0));
}

/** Deallocates the element vector, so the memory requirement should not be significant afer this*/
template<typename FG_ELEMENT>
void FullGrid<FG_ELEMENT>::deleteFullGrid() {
  isFGcreated_ = false;
  //    fullgridVector_.flush();
  fullgridVector_.clear();
}

/** evaluates the full grid on the specified coordinates
 * @param coords ND coordinates on the unit square [0,1]^D*/
template<typename FG_ELEMENT>
FG_ELEMENT FullGrid<FG_ELEMENT>::eval(std::vector<real>& coords) {
  assert(!isHierarchized_ && "cannot eval hierarchized fg!");

  IndexType ii, i, tmp_val, vv, nr;
  IndexType jj;
  real baseVal;
  real normcoord;
  // this value will be reseted, but just to avoid compiler warnings
  FG_ELEMENT ret_val = fullgridVector_[0];
  real intersect[127];
  IndexType aindex[63];

  // if there is a transformation then transform to the unit coordinates
  if (gridDomain_ != NULL) {
    assert(false && "not implemented");

    /* todo: commented because transform not compatible with types
     for( ii = dim_ - 1; ii >= 0; ii-- ){
     ( gridDomain_->get1DDomain( ii ) ).transformRealToUnit(
     coords[ii], coords[ii], levels_[ii],
     ( hasBoundaryPoints_[ii] == false ) );
     }
     */
  }

  // the coordinates are on the unit square
  for (ii = dim_ - 1; ii >= 0; ii--) {
    // the value has to be between 0.0 and 1.0
    //coords[ii] = (coords[ii] > 1.0) ? 1.0 : coords[ii];
    //coords[ii] = (coords[ii] < 0.0) ? 0.0 : coords[ii];

    if ((coords[ii] < 0.0) || (coords[ii] > 1.0))
      return 0.0;

    // scale to the reference [0,1] interval the intersection point
    normcoord = coords[ii] * powerOfTwo[levels_[ii]];
    aindex[ii] = static_cast<IndexType>(std::floor(normcoord));

    if (hasBoundaryPoints_[ii] == true) {
      aindex[ii] = (aindex[ii] < 0) ? 0 : aindex[ii];
      aindex[ii] =
        (aindex[ii] >= (nrPoints_[ii] - 1)) ?
        (nrPoints_[ii] - 2) : aindex[ii];
      //calculate the coordinates
      normcoord = normcoord - static_cast<real>(aindex[ii]);
    } else {
      aindex[ii] = (aindex[ii] < 0) ? 0 : aindex[ii];
      aindex[ii] =
        (aindex[ii] >= (nrPoints_[ii])) ? (nrPoints_[ii]) : aindex[ii];
      normcoord = normcoord - static_cast<real>(aindex[ii]);

      if (aindex[ii] <= 0) {
        aindex[ii] = 1;
        // extrapolation to the left, must be a negative number
        normcoord = coords[ii] * powerOfTwo[levels_[ii]] - 1;
      } else {
        // make extrapolation at the last cell (where there is no boundary point)
        if (aindex[ii] >= nrPoints_[ii] - 1) {
          aindex[ii] = nrPoints_[ii] - 1;
          // extrapolation to the right this will be a >=1 number
          normcoord = coords[ii] * powerOfTwo[levels_[ii]]
                      - static_cast<real>(aindex[ii]); // this should be 1 + ...
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
    intersect[2 * ii] = basis_->functionEval1(normcoord);
    intersect[2 * ii + 1] = basis_->functionEval2(normcoord);
  }

  nr = powerOfTwo[dim_];

  for (ii = 0; ii < nr; ii++) {
    baseVal = 1; //we will store here the coefficient for the next corner point of the cuboid
    i = 0;
    tmp_val = ii;
    vv = 0;

    // compute the "dim" dimensional basis function value(for one node) and the corresponding vector index
    for (jj = dim_ - 1; jj >= 0; jj--) {
      vv = (nrPoints_[jj] > 1) ? (tmp_val & 1) :
           0; // if we have only one point then use the same

      // todo: find bug or reimplement eval
      // quick fix for the following bug: no boundary and l=1
      // always evals zero. i think the intersect values are not
      // computed correctly, but do not fully understand the code
      if (nrPoints_[jj] == 1) {
        assert(!hasBoundaryPoints_[jj]);
        assert(aindex[jj] == 0);
        assert(coords[jj] <= 1.0 && coords[jj] >= 0.0);

        if (coords[jj] > 0.5)
          baseVal = baseVal * basis_->functionEval1(coords[jj]);
        else
          baseVal = baseVal * basis_->functionEval2(coords[jj]);
      } else
        baseVal = baseVal * intersect[2 * jj + vv];

      i = i * nrPoints_[jj] + aindex[jj] + vv;
      tmp_val = tmp_val >> 1;
    }

    // if this is the first point then set directly the value, otherwise add the value
    //COMBIGRID_OUT_LEVEL3( verb , "Vect elem i:" << i << " , vect:" << fullgridVector_[i]);
    ret_val =
      (ii == 0) ?
      (baseVal * fullgridVector_[i]) :
      (ret_val + baseVal * fullgridVector_[i]);
  }

  //COMBIGRID_OUT_LEVEL3( 4 , "FullGrid::eval RET VAL:" << ret_val);
  // just return the value
  return ret_val;
}

/** return the coordinates on the unit square
 * @param elemIndex [IN] index of the element in the vector
 * @param coords [OUT] the vector must be resized already*/
template<typename FG_ELEMENT>
void FullGrid<FG_ELEMENT>::getCoords(IndexType elemIndex,
                                     std::vector<real>& coords) const {
  // temporary variables
  IndexType ind = 0;
  IndexType tmp_add = 0;

  if (gridDomain_ == NULL) {
    for (DimType j = 0; j < dim_; j++) {
      ind = elemIndex % nrPoints_[j];
      elemIndex = elemIndex / nrPoints_[j];
      // set the coordinate based on if we have boundary points
      tmp_add = (hasBoundaryPoints_[j] == true) ? (0) : (1);
      coords[j] = static_cast<real>(ind + tmp_add)
                  * oneOverPowOfTwo[levels_[j]];
    }
  } else {
    // we have a valid Domain
    assert(false && "not implemented");
    /* todo: commented because transformunittoreal not working with IndexType
     for( DimType j = 0; j < dim_; j++ ){
     ind = elemIndex % (IndexType) ( nrPoints_[j] );
     elemIndex = elemIndex / (IndexType) ( nrPoints_[j] );
     // set the coordinate based on if we have boundary points
     tmp_add = ( hasBoundaryPoints_[j] == true ) ? ( 0 ) : ( 1 );
     ( gridDomain_->get1DDomain( j ) ).transformUnitToReal( levels_[j],
     ind + tmp_add,
     coords[j] );
     }
     */
  }
}

template<typename FG_ELEMENT>
void FullGrid<FG_ELEMENT>::getLI(IndexType elementIndex, LevelVector& levels,
                                 IndexVector& indexes) const {
  IndexType startindex, tmp_val;

  tmp_val = elementIndex;

  // first calculate intermediary indexes
  for (DimType k = 0; k < dim_; k++) {
    startindex = (hasBoundaryPoints_[k]) ? 0 : 1;
    indexes[k] = tmp_val % nrPoints_[k] + startindex;
    tmp_val = tmp_val / nrPoints_[k];
  }

  //The level and index of the element in the hashgridstorage are computed dividing by two the index and level in the fullgrid
  //until we obtain an impair number for the index, thus obtaining the level and index in the hierarchical basis (Aliz Nagy)
  // ...
  for (DimType k = 0; k < dim_; k++) {
    tmp_val = levels_[k];

    if (indexes[k] != 0) {
      // todo: these operations can be optimized
      while (indexes[k] % 2 == 0) {
        indexes[k] = indexes[k] / 2;
        tmp_val--;
      }
    } else {
      tmp_val = 0;
    }

    levels[k] = tmp_val;
  }

  // a little hack to allow switching on an alternative assignment of
  // the level vectors: the boundary points have level 1 and not level 0
  // as in the standard case
#ifdef ALT_LEVEL_VECTOR

  for ( DimType k = 0; k < dim_; k++ )
    if ( levels[k] == 0 )
      levels[k] = 1;

#endif
}

/** returns the vector index for one linear index
 * @param linIndex [IN] the linear index
 * @param axisIndex [OUT] the returned vector index */
template<typename FG_ELEMENT>
inline void FullGrid<FG_ELEMENT>::getVectorIndex(const IndexType linIndex,
    IndexVector& axisIndex) const {
  IndexType tmp = linIndex;

  for (int i = static_cast<int>(dim_) - 1; i >= 0; i--) {
    axisIndex[i] = tmp / (this->getOffset(i));
    tmp = tmp % this->getOffset(i);
  }
}

/** returns the linear index for one linear index
 * @param axisIndex [IN] the vector index */
template<typename FG_ELEMENT>
inline IndexType FullGrid<FG_ELEMENT>::getLinearIndex(
  const IndexVector& axisIndex) const {
  IndexType tmp = 0;

  for (int i = static_cast<int>(dim_) - 1; i >= 0; i--) {
    tmp = tmp + offsets_[i] * axisIndex[i];
  }

  return tmp;
}

/** sets the domain of the full grid */
template<typename FG_ELEMENT>
inline void FullGrid<FG_ELEMENT>::setDomain(GridDomain* gridDomain) const {
  gridDomain_ = gridDomain;
}

/** returns the domain of the full grid */
template<typename FG_ELEMENT>
inline const GridDomain*
FullGrid<FG_ELEMENT>::getDomain() const {
  return gridDomain_;
}

/** returns the domain of the full grid */
template<typename FG_ELEMENT>
inline GridDomain*
FullGrid<FG_ELEMENT>::getDomain() {
  return gridDomain_;
}

/** returns pointer to the basis function */
template<typename FG_ELEMENT>
inline const BasisFunctionBasis*
FullGrid<FG_ELEMENT>::getBasisFct() const {
  return basis_;
}

/** returns the dimension of the full grid */
template<typename FG_ELEMENT>
inline DimType FullGrid<FG_ELEMENT>::getDimension() const {
  return dim_;
}

/** the getters for the full grid vector */
template<typename FG_ELEMENT>
inline std::vector<FG_ELEMENT>&
FullGrid<FG_ELEMENT>::getElementVector() {
  return fullgridVector_;
}

template<typename FG_ELEMENT>
inline const std::vector<FG_ELEMENT>&
FullGrid<FG_ELEMENT>::getElementVector() const {
  return fullgridVector_;
}

/** return the offset in the full grid vector of the dimension */
template<typename FG_ELEMENT>
inline IndexType FullGrid<FG_ELEMENT>::getOffset(DimType i) const {
  return offsets_[i];
}

template<typename FG_ELEMENT>
inline const IndexVector&
FullGrid<FG_ELEMENT>::getOffsets() const {
  return offsets_;
}

/** return the level vector */
template<typename FG_ELEMENT>
inline const LevelVector&
FullGrid<FG_ELEMENT>::getLevels() const {
  return levels_;
}

/** flag shows if the full grid is created */
template<typename FG_ELEMENT>
inline bool FullGrid<FG_ELEMENT>::isGridCreated() const {
  return isFGcreated_;
}

/** function to check whether a fullgrid is hierarchized */
template<typename FG_ELEMENT>
inline bool FullGrid<FG_ELEMENT>::isHierarchized() const {
  return isHierarchized_;
}

/** set the isHierarchized flag */
template<typename FG_ELEMENT>
inline void FullGrid<FG_ELEMENT>::setHierarchized() {
  isHierarchized_ = true;
}

/** returns the number of elements in the full grid */
template<typename FG_ELEMENT>
inline IndexType FullGrid<FG_ELEMENT>::getNrElements() const {
  return nrElements_;
}

/** returns the number of elements in the full grid */
template<typename FG_ELEMENT>
inline IndexType FullGrid<FG_ELEMENT>::getNrElementsNoBoundary() const {
  return nrElementsNoBoundary_;
}

/** number of points per dimension i */
template<typename FG_ELEMENT>
inline IndexType FullGrid<FG_ELEMENT>::length(DimType i) const {
  return nrPoints_[i];
}

/** vector of flags to show if the dimension has boundary points*/
template<typename FG_ELEMENT>
inline const std::vector<bool>&
FullGrid<FG_ELEMENT>::returnBoundaryFlags() const {
  return hasBoundaryPoints_;
}

/** copies the input vector to the full grid vector
 * in most cases the add function would be more save
 * @param in [IN] input vector*/
template<typename FG_ELEMENT>
void FullGrid<FG_ELEMENT>::setElementVector(const std::vector<FG_ELEMENT>& in) {
  assert(in.size() == static_cast<size_t>(nrElements_));
  fullgridVector_ = in;
  isFGcreated_ = true;
}

template<typename FG_ELEMENT>
inline FG_ELEMENT*
FullGrid<FG_ELEMENT>::getData() {
  return &fullgridVector_[0];
}

template<typename FG_ELEMENT>
inline const FG_ELEMENT*
FullGrid<FG_ELEMENT>::getData() const {
  return &fullgridVector_[0];
}

template<typename FG_ELEMENT>
void FullGrid<FG_ELEMENT>::gridEval(FullGrid<FG_ELEMENT>& dst) const {
  assert(!isHierarchized_ && "cannot eval hierarchized fg!");
  assert(getDimension() == dst.getDimension());
  assert(returnBoundaryFlags() == dst.returnBoundaryFlags());

  if (!isGridCreated())
    return;

  dst.createFullGrid();

  // loop over all elements of dst and eval fg at the gridpoints of dst
  for (IndexType i = 0; i < dst.getNrElements(); ++i) {
    std::vector<real> c(getDimension());
    dst.getCoords(i, c);

    dst.getData()[i] = eval(c);
  }
}

template<typename FG_ELEMENT>
void FullGrid<FG_ELEMENT>::save(std::string& filename) {
  std::ofstream ofs(filename.c_str(), std::ios_base::binary);
  //boost::archive::binary_oarchive oa(ofs);
  boost::archive::text_oarchive oa(ofs);
  oa << *this;
}

/* add another full grid to the current full grid
 *
 * if the grids differ in discretization or boundary flags interpolation
 * based on the eval function is used
 */
template<typename FG_ELEMENT>
void FullGrid<FG_ELEMENT>::add(FullGrid<FG_ELEMENT>& fg, real coeff) {
  assert(!isHierarchized_ && "can only add an unhierarchized fg!");

  // todo: a simple return should be enough here
  // however, probably there is an error somewhere else if this happens
  assert(fg.isGridCreated());

  assert(this->getDimension() == fg.getDimension());

  if (!this->isGridCreated())
    this->createFullGrid();

  // if boundary flags and levelvectors equal, we have the same grid
  // structure. thus we can directly add up the data. much faster.
  if (this->getLevels() == fg.getLevels()
      && this->returnBoundaryFlags() == fg.returnBoundaryFlags()) {
    const std::vector<FG_ELEMENT>& src = fg.getElementVector();
    std::vector<FG_ELEMENT>& dst = this->getElementVector();

    for (size_t i = 0; i < dst.size(); ++i) {
      dst[i] += coeff * src[i];
    }
  } else {
    std::vector<real> coords(this->getDimension());
    std::vector<FG_ELEMENT>& dst = this->getElementVector();

    for (IndexType i = 0; i < this->getNrElements(); ++i) {
      this->getCoords(i, coords);

      dst[i] += coeff * fg.eval(coords);
    }

  }
}

/* add another full grid to the current full grid
 *
 * special treatment for GENE grids which have fourier coefficients
 * in x-direction: no interpolation happens in x-direction. by using
 * a coordinate transformation we make sure the grids fit toghether in
 * the desired way
 */
template<typename FG_ELEMENT>
void FullGrid<FG_ELEMENT>::addGENE(const FullGrid<FG_ELEMENT>& fg, real coeff) {
  assert(!isHierarchized_ && "can only add an unhierarchized fg!");

  // todo: a simple return should be enough here
  // however, probably there is an error somewhere else if this happens
  assert(fg.isGridCreated());

  assert(this->getDimension() == fg.getDimension());

  if (!this->isGridCreated())
    this->createFullGrid();

  // if boundary flags and levelvectors equal, we have the same grid
  // structure. thus we can directly add up the data. much faster.
  if (this->getLevels() == fg.getLevels()
      && this->returnBoundaryFlags() == fg.returnBoundaryFlags()) {
    const std::vector<FG_ELEMENT>& src = fg.getElementVector();
    std::vector<FG_ELEMENT>& dst = this->getElementVector();

    for (IndexType i = 0; i < dst.size(); ++i) {
      dst[i] += coeff * src[i];
    }
  } else {
    std::vector<real> coords(this->getDimension());
    std::vector<FG_ELEMENT>& dst = this->getElementVector();

    const real lx_src = static_cast<real>(fg.getLevels()[0]);
    const real lx_dst = static_cast<real>(this->getLevels()[0]);
    const real a = std::pow(2.0, lx_src) / std::pow(2.0, lx_dst);
    const real b = 0.5 - 0.5 * a;
    const real aa = 1.0 / a;

    for (IndexType i = 0; i < this->getNrElements(); ++i) {
      this->getCoords(i, coords);

      // coordinate transformation in x-direction
      // x_src = 1/a * ( x_dst - b ) with
      // a = 2**lx_src / 2**lx_dst
      // b = 0.5 - 0.5*a
      coords[0] = aa * (coords[0] - b);

      dst[i] += coeff * fg.eval(coords);
    }
  }
}

// serialize
template<typename FG_ELEMENT>
template<class Archive>
void FullGrid<FG_ELEMENT>::serialize(Archive& ar, const unsigned int version) {
  ar& dim_;
  ar& nrElements_;
  ar& nrElementsNoBoundary_;
  ar& isFGcreated_;
  ar& levels_;
  ar& nrPoints_;
  ar& hasBoundaryPoints_;
  ar& offsets_;
  ar& fullgridVector_;
  ar& isHierarchized_;
}

template<typename FG_ELEMENT>
inline real FullGrid<FG_ELEMENT>::getlpNorm(int p) {
  assert(isGridCreated());
  assert(p > 0);

  // special case maximum norm
  if (p == 0) {
    std::vector<FG_ELEMENT>& data = getElementVector();
    real max = 0.0;

    for (size_t i = 0; i < data.size(); ++i) {
      if (std::abs(data[i]) > max)
        max = std::abs(data[i]);
    }

    return max;
  } else {
    real p_f = static_cast<real>(p);
    std::vector<FG_ELEMENT>& data = getElementVector();
    real res = 0.0;

    for (size_t i = 0; i < data.size(); ++i) {
      real abs = std::abs(data[i]);
      res += std::pow(abs, p_f);
    }

    return std::pow(res, 1.0 / p_f);
  }
}

template<typename FG_ELEMENT>
inline real FullGrid<FG_ELEMENT>::normalizelp(int p) {
  real norm = getlpNorm(p);

  for (IndexType i = 0; i < getNrElements(); ++i)
    getData()[i] = getData()[i] / norm;

  return norm;
}

template<typename FG_ELEMENT>
inline const IndexVector& FullGrid<FG_ELEMENT>::getSizes() const {
  return nrPoints_;
}

// 2d output
template<typename FG_ELEMENT>
void FullGrid<FG_ELEMENT>::print2D(std::ostream& os) const {
  assert(dim_ == 2);

  // loop over rows i -> dim2
  for (IndexType i = 0; i < nrPoints_[1]; ++i) {
    IndexType offset = offsets_[1] * i;

    for (IndexType j = 0; j < nrPoints_[0]; ++j) {
      os << fullgridVector_[offset + j] << "\t";
    }

    os << std::endl;
  }
}

// 1d output
template<typename FG_ELEMENT>
void FullGrid<FG_ELEMENT>::print1D(std::ostream& os) const {
  assert(dim_ == 1);

  for (IndexType j = 0; j < nrPoints_[0]; ++j) {
    os << fullgridVector_[j] << "\t";
  }

  os << std::endl;
}

// n-dim output. will print 2-dimensional slices of domain
template<typename FG_ELEMENT>
void FullGrid<FG_ELEMENT>::print3D(std::ostream& os) const {
  for (IndexType k = 0; k < nrPoints_[2]; ++k) {
    assert(dim_ == 3);

    IndexType offsetZ = offsets_[2] * k;

    // loop over rows i -> dim2
    for (IndexType i = 0; i < nrPoints_[1]; ++i) {
      IndexType offset = offsetZ + offsets_[1] * i;

      for (IndexType j = 0; j < nrPoints_[0]; ++j) {
        os << fullgridVector_[offset + j] << "\t";
      }

      os << std::endl;
    }
  }
}

template<typename FG_ELEMENT>
inline void FullGrid<FG_ELEMENT>::print(std::ostream& os) const {
  if (dim_ == 1)
    print1D(os);
  else if (dim_ == 2)
    print2D(os);
  else if (dim_ == 3)
    print3D(os);
  else
    assert(false && "Maximum dimension for printing is 2");
}

// output operator
template<typename FG_ELEMENT>
inline std::ostream& operator<<(std::ostream& os,
                                const FullGrid<FG_ELEMENT>& fg) {
  fg.print(os);

  return os;
}

}      //namespace
#endif /* COMBIFULLGRID_HPP_ */
