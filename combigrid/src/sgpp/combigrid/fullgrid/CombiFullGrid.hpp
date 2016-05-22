// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#ifndef COMBIFULLGRID_HPP_
#define COMBIFULLGRID_HPP_

#include <sgpp/combigrid/basisfunction/CombiBasisFunctionBasis.hpp>
#include <sgpp/combigrid/basisfunction/CombiLinearBasisFunction.hpp>
#include <sgpp/combigrid/domain/CombiGridDomain.hpp>
#include <sgpp/combigrid/domain/StretchingFactory.hpp>
#include <sgpp/combigrid/utils/combigrid_utils.hpp>

#include <vector>
#include <algorithm>

namespace combigrid {

// quick n dirty implementation of dirac delta function...
template <typename _Tp>
inline _Tp delta(const _Tp& a) {
  return (_Tp)(a == 0);
}

/** The full grid class which is the main building block of the combi grid <br>
 *  It is important that the grid will actually occupy memory when the
 * createFullGrid()
 *  method is called. <br>
 *  The index of a gridpoint in the full grid is given by the formula : <br>
 *  ind = i0 + i1*N0 + i2*N0*N1 + ... + id*N0*N1*N2*...*Nd, where i0,i1,i2,...
 * are the indexes in every dimension,
 *  and Nk is the number of gridpoints in direction k. <br>
 *  For more infos you can look at the "getVectorIndex" and "getLinearIndex"
 * functions and their implementation. <br>
 *  Nk=2^level[k]+1 for every k for the directions with boundary points and
 * Nk=2^level[k]-1 for the directions without boundary. <br>
 *  <br>
 *  The full grid can also be scaled which is done with a separate
 * "combigrid::GridDomain" object. <br>
 * */
template <typename FG_ELEMENT>
class FullGrid {
 public:
  /** simplest Ctor with homogeneous  levels */
  FullGrid(int dim, int level, bool hasBdrPoints = true, const BasisFunctionBasis* basis = NULL) {
    // set the basis function for the full grid
    if (basis == NULL)
      basis_ = LinearBasisFunction::getDefaultBasis();
    else
      basis_ = basis;

    sgppIndex_ = new std::vector<int>(0);

    isFGcreated_ = false;
    dim_ = dim;
    levels_.resize(dim, level);

    // init default gridDomain_
    CombiEquidistantStretching stretching;
    std::vector<double> min(dim, 0.0);
    std::vector<double> max(dim, 1.0);
    gridDomain_ = new GridDomain(dim, levels_, min, max, stretching);

    hasBoundaryPoints_.resize(dim, hasBdrPoints);
    nrElements_ = 1;
    offsets_.resize(dim_);
    nrPoints_.resize(dim_);

    for (int j = 0; j < dim_; j++) {
      nrPoints_[j] = ((hasBoundaryPoints_[j] == true) ? (powerOfTwo[levels_[j]] + 1)
                                                      : (powerOfTwo[levels_[j]] - 1));
      offsets_[j] = nrElements_;
      nrElements_ = nrElements_ * nrPoints_[j];
    }
  }

  /** dimension adaptive Ctor */
  FullGrid(int dim, const std::vector<int>& levels, bool hasBdrPoints = true,
           const BasisFunctionBasis* basis = NULL) {
    // set the basis function for the full grid
    if (basis == NULL)
      basis_ = LinearBasisFunction::getDefaultBasis();
    else
      basis_ = basis;

    sgppIndex_ = new std::vector<int>(0);
    dim_ = dim;
    isFGcreated_ = false;
    levels_ = levels;
    // init default gridDomain_
    CombiEquidistantStretching stretching;
    std::vector<double> min(dim, 0.0);
    std::vector<double> max(dim, 1.0);
    gridDomain_ = new GridDomain(dim, levels_, min, max, stretching);

    hasBoundaryPoints_.resize(dim, hasBdrPoints);
    nrElements_ = 1;
    offsets_.resize(dim_);
    nrPoints_.resize(dim_);

    for (int j = 0; j < dim_; j++) {
      nrPoints_[j] = ((hasBoundaryPoints_[j] == true) ? (powerOfTwo[levels_[j]] + 1)
                                                      : (powerOfTwo[levels_[j]] - 1));
      offsets_[j] = nrElements_;
      nrElements_ = nrElements_ * nrPoints_[j];
    }
  }

  /** dimension adaptive Ctor */
  FullGrid(int dim, const std::vector<int>& levels, const std::vector<bool>& hasBdrPoints,
           const BasisFunctionBasis* basis = NULL) {
    // set the basis function for the full grid
    if (basis == NULL)
      basis_ = LinearBasisFunction::getDefaultBasis();
    else
      basis_ = basis;

    sgppIndex_ = new std::vector<int>(0);
    dim_ = dim;
    isFGcreated_ = false;
    levels_ = levels;
    // init default gridDomain_
    CombiEquidistantStretching stretching;
    std::vector<double> min(dim, 0.0);
    std::vector<double> max(dim, 1.0);
    gridDomain_ = new GridDomain(dim, levels_, min, max, stretching);

    hasBoundaryPoints_ = hasBdrPoints;
    nrElements_ = 1;
    offsets_.resize(dim_);
    nrPoints_.resize(dim_);

    for (int j = 0; j < dim_; j++) {
      nrPoints_[j] = ((hasBoundaryPoints_[j] == true) ? (powerOfTwo[levels_[j]] + 1)
                                                      : (powerOfTwo[levels_[j]] - 1));
      offsets_[j] = nrElements_;
      nrElements_ = nrElements_ * nrPoints_[j];
    }
  }

  virtual ~FullGrid() {
    delete sgppIndex_;

    if (gridDomain_ != NULL) delete gridDomain_;
  }

  /** allocates the memory for the element vector of the full grid <br>
   * Only after this needs the full grid considerable amount of memory */
  void createFullGrid() {
    if (!isFGcreated_) {
      isFGcreated_ = true;  // add a check if the fullgrid is created?
      fullgridVector_.resize(nrElements_);
    }
  }

  /** Deallocates the element vector, so the memory requirement should not be
   * significant afer this*/
  void deleteFullGrid() {
    isFGcreated_ = false;
    //    fullgridVector_.flush();
    fullgridVector_.clear();
    delete sgppIndex_;
    sgppIndex_ = new std::vector<int>(0);
  }

  /** evaluates the full grid on the specified coordinates
   * @param coords N-D coordinates on the unit square [0,1]^D*/
  FG_ELEMENT eval(std::vector<double>& coords) const {
    int ii, i, tmp_val, vv, nr;
    int jj;
    FG_ELEMENT baseVal;
    double normcoord;

    // this value will be reseted, but just to avoid compiler warnings
    FG_ELEMENT ret_val = fullgridVector_[0];
    double intersect[127];  // alternatievly [127]
    int aindex[63];         // alternatievly [63]

    // if there is a transformation then transform to the unit coordinates
    if (gridDomain_ != NULL) {
      for (ii = dim_ - 1; ii >= 0; ii--) {
        (gridDomain_->get1DDomain(ii))
            .transformRealToUnit(coords[ii], coords[ii], levels_[ii],
                                 (hasBoundaryPoints_[ii] == false));
      }
    }

    // the coordinates are on the unit square
    for (ii = dim_ - 1; ii >= 0; ii--) {
      // the value has to be between 0.0 and 1.0
      // coords[ii] = (coords[ii] > 1.0) ? 1.0 : coords[ii];
      // coords[ii] = (coords[ii] < 0.0) ? 0.0 : coords[ii];
      // if ((coords[ii]<0.0)||(coords[ii]>1.0)) return 0.0;
      // scale to the reference [0,1] interval the intersection point
      normcoord = coords[ii] * powerOfTwo[levels_[ii]];
      aindex[ii] = static_cast<int>(std::floor(normcoord));

      if (hasBoundaryPoints_[ii] == true) {
        aindex[ii] = (aindex[ii] < 0) ? 0 : aindex[ii];
        aindex[ii] = (aindex[ii] >= static_cast<int>(nrPoints_[ii] - 1))
                         ? static_cast<int>(nrPoints_[ii] - 2)
                         : aindex[ii];
        // calculate the coordinates
        normcoord = normcoord - static_cast<double>(aindex[ii]);
      } else {
        aindex[ii] = (aindex[ii] < 0) ? 0 : aindex[ii];
        aindex[ii] = (aindex[ii] >= static_cast<int>(nrPoints_[ii]))
                         ? static_cast<int>(nrPoints_[ii])
                         : aindex[ii];
        normcoord = normcoord - static_cast<double>(aindex[ii]);

        if (aindex[ii] <= 0) {
          aindex[ii] = 1;
          // extrapolation to the left, must be a negative number
          normcoord = coords[ii] * powerOfTwo[levels_[ii]] - 1;
          // COMBIGRID_OUT_LEVEL3( verb , "FullG L , aindex[ii]:" <<
          // aindex[ii]);
        } else {
          // make extrapolation at the last cell (where there is no boundary
          // point)
          if (aindex[ii] >= static_cast<int>(nrPoints_[ii]) - 1) {
            aindex[ii] = static_cast<int>(nrPoints_[ii]) - 1;
            // extrapolation to the right this will be a >=1 number
            normcoord = coords[ii] * powerOfTwo[levels_[ii]] -
                        static_cast<double>(aindex[ii]);  // this should be 1 + ...
            // COMBIGRID_OUT_LEVEL3( verb , "FullG R , aindex[ii]:" <<
            // aindex[ii]);
          }
        }

        // we have to subtract one in case of no boundary points
        aindex[ii] = aindex[ii] - 1;
        // in some special cases this might become negative, so limit to zero
        // (when level = 1 and no boundary point )
        aindex[ii] = (aindex[ii] < 0) ? 0 : aindex[ii];
      }

      // COMBIGRID_OUT_LEVEL3( verb , "FullGrid::eval ii:" << ii << " ,
      // coords[ii]:" << coords[ii] << " , aindex[ii]:" << aindex[ii]
      //                     << " , level[ii]:" << levels_[ii] << " ,
      //                     normcoord:" << normcoord );
      // evaluate the basis functions on the 1D coordinates
      intersect[2 * ii] = basis_->functionEval1(normcoord);
      intersect[2 * ii + 1] = basis_->functionEval2(normcoord);
    }

    nr = powerOfTwo[dim_];

    for (ii = 0; ii < nr; ii++) {
      baseVal = 1;  // we will store here the coefficient for the next corner
                    // point of the cuboid
      i = 0;
      tmp_val = ii;
      vv = 0;

      // compute the "dim" dimensional basis function value(for one node) and
      // the corresponding vector index
      for (jj = dim_ - 1; jj >= 0; jj--) {
        vv =
            (nrPoints_[jj] > 1) ? (tmp_val & 1) : 0;  // if we have only one point then use the same
        baseVal = baseVal * (FG_ELEMENT)intersect[2 * jj + vv];
        i = i * nrPoints_[jj] + aindex[jj] + vv;
        tmp_val = tmp_val >> 1;
      }

      // if this is the first point then set directly the value, otherwise add
      // the value
      // COMBIGRID_OUT_LEVEL3( verb , "Vect elem i:" << i << " , vect:" <<
      // fullgridVector_[i]);
      ret_val =
          (ii == 0) ? (baseVal * fullgridVector_[i]) : (ret_val + baseVal * fullgridVector_[i]);
    }

    // COMBIGRID_OUT_LEVEL3( verb , "FullGrid::eval RET VAL:" << ret_val);
    // just return the value
    return ret_val;
  }

  /**
   *  PeTz:Implementation of the prolongation/reduction operations for a
   *d-dimensional FullGrid!
   *  The function:
   *
   *     transform_gird(std::vector<int> new_levels ) : FullGrid<FG_ELEMENT>* ;
   *
   *  returns an extension (G_new) of the current grid (G) with levels vector
   *new_levels and grid values which
   *  are linearly interpolated from G
   *
   * @param new_levels an integer vector of size "dim" specifying the desired
   *levels for eahc of the dimensions of the new grid.
   * @return returns a FullGrid<FG_ELEMENT>* with levels vector equal to
   *"new_levels" and linearly interpolated grid values
   *
   *
   *  Computational complexity is O( N 2^dim) where N is the total number of
   *points in the new grid!
   */

  FullGrid<FG_ELEMENT>* transform_grid(std::vector<int> new_levels) {
    if (new_levels.size() != levels_.size()) {
      COMBIGRID_OUT_ERR(
          "Newly requested grid dimensions incompatible with existing grid's "
          "dimensions! Please reconsider your inputs and"
          "try again. Returning NULL!",
          __FILE__, __LINE__)
      return NULL;
    }

    // check for non-positive values of the new_levels input vector
    for (int i = 0; i < new_levels.size(); i++) {
      if (new_levels[i] < 0) {
        COMBIGRID_OUT_ERR(
            "Not all newly requested grid levels are non-negative.! Please "
            "reconsider your inputs and"
            "try again. Returning NULL!",
            __FILE__, __LINE__)
        return NULL;
      }
    }

    // if input is all right create new grid instance and interpolate
    // elements....
    FullGrid<FG_ELEMENT>* new_grid =
        new FullGrid<FG_ELEMENT>(dim_, new_levels, hasBoundaryPoints_, basis_);
    new_grid->createFullGrid();  // allocate mem for the new fullGrid

    // now set the new_grid domain
    if (gridDomain_ != NULL) {
      // in case of inhomogeneous stretching...
      std::vector<AbstractStretchingMaker*> stretchings(dim_, NULL);

      for (int d = 0; d < dim_; d++)
        stretchings[d] =
            combigrid::createStretchingMaker(gridDomain_->get1DDomain(d).getStretchingType());

      std::vector<double> min = gridDomain_->getMin();
      std::vector<double> max = gridDomain_->getMax();
      GridDomain* new_domain = new GridDomain(dim_, new_levels, min, max, stretchings);
      new_grid->setDomain(new_domain);
    }

    int nr_elems = new_grid->getNrElements();
    std::vector<double> k_coords(dim_, 0.0);  // current point coordinates
    std::vector<int> k_idx(dim_, 0);

    std::vector<int> bb_idx(dim_, 0);              // bounding box indices (only left)
    std::vector<double> bb_coords(2 * dim_, 0.0);  // bounding box coords (left and right)

    // get the next element of the new grid
    // speed-up somehow :> ???
    for (int j_k = 0; j_k < nr_elems; j_k++) {
      // take the coordinates of the e-th element of the new grid
      // O(dim)

      new_grid->getVectorIndex(j_k, k_idx);
      new_grid->getStandardCoordinates(j_k, k_coords);

      // take the bounding box of e O(dim)
      bounding_box(k_idx, new_levels, bb_idx, bb_coords);
      // calculate and store the interpolation coeffs for F[e] in each separate
      // dimension
      std::vector<double> interpol_coeffs(2 * dim_, 0.0);

      // O(dim)
      for (int d = 0; d < dim_; d++) {
        double x_d0 = bb_coords[d * 2];
        double x_d1 = bb_coords[d * 2 + 1];
        interpol_coeffs[d * 2] = calculate_interpolation_coeff(x_d0, x_d1, k_coords[d], -1);
        interpol_coeffs[d * 2 + 1] = calculate_interpolation_coeff(x_d0, x_d1, k_coords[d], +1);
      }

      // now perform the interpolation sum!
      double F_e = 0.0;  // the grid value at point "e"

      //      O(2^dim)
      for (int j = 0; j < powerOfTwo[dim_]; j++) {
        // stores the interpolation for the j-th vertex of the bounding box
        double a_j = 1.0;
        // we need this touple to extract the global linear coordinate of the
        // j-the vertex of the bounding box...
        std::vector<int> f_glob_idx(dim_, 0);

        for (int d = 0; d < dim_; d++) {
          int b_d_j = (j & (1 << d)) >> d;  // get the d-th bit of j with a
                                            // 000010000 mask (1 is at d-th
                                            // position).
          a_j *= interpol_coeffs[d * 2 + b_d_j];
          f_glob_idx[d] = bb_idx[d] + b_d_j;
        }

        int f_glob_j = getLinearIndex(f_glob_idx);
        F_e += a_j * fullgridVector_[f_glob_j];
      }

      (new_grid->getElementVector())[j_k] = F_e;
    }

    return new_grid;
  }

  /**
   *  PeTz: A quick and dirty function mapping the j_k-th index of a 1D grid of
   *level k( denoted G_k)
   *
   *  to a corresponding element of 1D gird of level l (denoted G_l),
   *
   *  The function handles two different situations:
   *  a) when k>=l, i.e. we are dealing with restriction from G_k -> G_l, then
   *  we have (due to nestedness) that the the j_k-th point of G_k,
   *  and the j_k/2^(k-l) the point from G_l have the same coordinates,
   *  provided that:
   *    mod(j_k,2^(k-l)) == 0
   *
   *  Even though the above condition is sometimes violated , this
   *"transform_idx" function
   *  always returns the rounded value of j_k/2^(k-l), in effect mapping the
   *j_k-th point to the left boundary of
   *  its surrounding interval in G_l:
   *
   *    G_k (k =2 )   ->   G_l (k = 1)
   *    X--X--X--X--X   ->  x-----x-----x
   *
   *
   *
   *  b) when k<l , i.e. we are dealing with prolongation along this dimension,
   *  then we have that the j_k-th point of G_k  and the
   *  j_k*2^(k-l)-th point of G_l, have the same coordinates.
   *
   *    G_k (k=1 )   ->  G_l (k = 2)
   *    X-----x-----x  ->  X--X--X--X--X
   *
   */

  inline int bb_left(int j_k, int l, int k) const {
    int dL = k - l;
    return (dL >= 0 ? (j_k / powerOfTwo[dL]) : j_k * powerOfTwo[-dL]) - j_k / powerOfTwo[k];
    // at the end we subtract j_k/powerOfTwo[k] to ensure that points on the
    // right boundary (i.e. when j_k = 2^k
    // are mapped onto the point 2^l-1, instead of 2^l!!! this is done because
    // bb_left is supposed to give us
    // the index of the left boundary of the bounding box of the coordinate
    // x_{j_k}...
  }

  /** return the coordinates on the unit square
   * @param elemIndex [IN] index of the element in the vector
   * @param coords [OUT] the vector must be resized already*/
  void getStandardCoordinates(int elemIndex, std::vector<double>& coords) const {
    // temporary variables
    // int verb = 6;
    int ind = 0;
    int tmp_add = 0;

    if (gridDomain_ == NULL) {  // default we are getting coordinates in the [0;1]^dim cube...
      for (int j = 0; j < dim_; j++) {
        ind = elemIndex % static_cast<int>(nrPoints_[j]);
        elemIndex = elemIndex / static_cast<int>(nrPoints_[j]);
        // set the coordinate based on if we have boundary points
        tmp_add = (hasBoundaryPoints_[j] == true) ? (0) : (1);
        coords[j] = (static_cast<double>(ind + tmp_add)) * oneOverPowOfTwo[levels_[j]];
      }
    } else {
      // we have a valid Domain
      for (int j = 0; j < dim_; j++) {
        ind = elemIndex % static_cast<int>(nrPoints_[j]);
        elemIndex = elemIndex / static_cast<int>(nrPoints_[j]);
        // set the coordinate based on if we have boundary points
        tmp_add = (hasBoundaryPoints_[j] == true) ? (0) : (1);
        (gridDomain_->get1DDomain(j)).transformUnitToReal(levels_[j], ind + tmp_add, coords[j]);
        // COMBIGRID_OUT_LEVEL3( verb , "FullGrid::getStandardCoordinates j:" << j << " ,
        // coords[j]:" << coords[j]);
      }
    }
  }

  /* returns the LI (level,index) notation for a given element in the full grid
   * @param elementIndex [IN] the linear index of the element
   * @param levels [OUT] the levels of the point in the LI notation
   * @param indexes [OUT] the indexes of the point in the LI notation */
  // void getLI_std(int elementIndex , std::vector<int>& levels ,
  // std::vector<int>& indexes) const {
  //  getLI( elementIndex , &(levels[0]) , &(indexes[0]) );
  //}
  /** returns the LI (level,index) notation for a given element in the full grid
   * @param elementIndex [IN] the linear index of the element
   * @param levels [OUT] the levels of the point in the LI notation
   * @param indexes [OUT] the indexes of the point in the LI notation */
  void getLI(int elementIndex, std::vector<int>& levels, std::vector<int>& indexes) const {
    int k, startindex, tmp_val;

    tmp_val = elementIndex;

    // first calculate intermediary indexes
    for (k = 0; k < dim_; k++) {
      startindex = (hasBoundaryPoints_[k]) ? 0 : 1;
      indexes[k] = tmp_val % nrPoints_[k] + startindex;
      tmp_val = tmp_val / nrPoints_[k];
    }

    // The level and index of the element in the hashgridstorage are computed
    // dividing by two the index and level in the fullgrid
    // until we obtain an impair number for the index, thus obtaining the level
    // and index in the hierarchical basis (Aliz Nagy)
    // ...
    for (k = 0; k < dim_; k++) {
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
  }

  /** returns the vector index for one linear index! Please note that
   *irrespective of whether or not
   * the full grid data storage includes the boundary points, the vector indices
   *ALWAYS assume that the
   * the abscissa point lying on the left boundary has index 0 and the one lying
   *on the right boundary
   * has index 2^level + 1.
   *
   * Naturally the getVectorIndex method follows that convention. I.e. for a
   *full grid of say 2 dimensions
   * and for e.g. 3 points per dimension,
   *  the linear index enumeration when bdries are included will look as
   *follows:
   *      _____
   *  2  |6 7 8|
   *  1  |3 4 5|
   *  0  |0 1 2|
   *
   *     0 1 2
   *   for the above construction, a call to getVectorIndex(4,vector) shall
   *return the vector (1,1)
   *
   *  if we exclude the boundaries the grid looks as follows
   *  2 * * *
   *  1 *|0|*
   *  0 * * *
   *    0 1 2
   *
   *  for this case a call to getVectorIndex should also return
   *  (1,1)...
   *
   *
   * @param linIndex [IN] the linear index
   * @param axisIndex [OUT] the returned vector index */
  inline void getVectorIndex(const int linIndex, std::vector<int>& axisIndex) const {
    int tmp = linIndex;

    for (int i = dim_ - 1; i >= 0; i--) {
      int bdry = hasBoundaryPoints_[i] ? 0 : 1;
      axisIndex[i] = tmp / (this->getOffset(i)) + bdry;
      tmp = tmp % this->getOffset(i);
    }
  }

  /** returns the linear index for one linear index
   * @param axisIndex [IN] the vector index */
  inline int getLinearIndex(const std::vector<int>& axisIndex) const {
    int tmp = 0;

    for (int i = dim_ - 1; i >= 0; i--) {
      tmp = tmp + offsets_[i] * axisIndex[i];
    }

    return tmp;
  }

  /** sets the domain of the full grid */
  inline void setDomain(GridDomain* gridDomain) { gridDomain_ = gridDomain; }

  /** returns the domain of the full grid */
  inline const GridDomain* getDomain() const { return gridDomain_; }

  /** returns the domain of the full grid */
  inline GridDomain* getDomain() { return gridDomain_; }

  /** returns pointer to the basis function */
  inline const BasisFunctionBasis* getBasisFct() const { return basis_; }

  /** returns the dimension of the full grid */
  inline int getDimension() const { return dim_; }

  /** the getters for the full grid vector */
  inline std::vector<FG_ELEMENT>& getElementVector() { return fullgridVector_; }

  inline const std::vector<FG_ELEMENT>& getElementVector() const { return fullgridVector_; }

  /** return the offset in the full grid vector of the dimension */
  inline int getOffset(int i) const { return offsets_[i]; }
  inline const std::vector<int>& getOffsets() const { return offsets_; }

  /** return the level vector */
  inline const std::vector<int>& getLevels() const { return levels_; }
  inline std::vector<int>& getLevels() { return levels_; }

  /** flag shows if the full grid is created */
  inline bool isGridCreated() const { return isFGcreated_; }

  /** returns the number of elements in the full grid */
  inline int getNrElements() const { return nrElements_; }

  /** number of points per dimension i */
  inline int length(int i) const { return nrPoints_[i]; }

  /** return the vector for faster combination of the full grids <br>
   * this will be used in the Converter */
  inline std::vector<int>& getSGppIndex() const { return (*sgppIndex_); }

  /** vector of flags to show if the dimension has boundary points*/
  inline const std::vector<bool>& returnBoundaryFlags() const { return hasBoundaryPoints_; }

  /** copies the input vector to the full grid vector
   * @param in [IN] input vector*/
  void setElementVector(std::vector<FG_ELEMENT> in) { fullgridVector_ = in; }

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

  /**
   * NOTE!!! ALREADY EXISTS IN CombiLinearBasisFunction.hpp
   *
   * PeTz:
   *
   * a simple  function returning the 1-D interpolation coefficient at point x
   *within the interval [x_0;x_1]
   * @param x_0 - the coordinate of the left neighbour of x
   * @param x_1 - the coordinate of the right neighbour of x
   * @param x - the coordinates of our sampling (interpolation) point
   * @param direction - an integer specifying the boundary point we need the
   *interpolation coefficient of. if direction is a negative number
   * , the method returns the interpolation coefficient for the left and if
   *positive, for right boundary.
   * @return a double value in the range [0.0;1.0] specifying the 1D linear
   *interpolation coefficient in the respective direction
   *
   * Needless to say, for the interpolation to make sense we need to ensure that
   * x_0 <= x <= x_1 and that direction != 0
   *
   * In case of one of those condition being violated, the method prints a
   *warning and returns 0.0.  !!!!
   *
   */

  double calculate_interpolation_coeff(double x_0, double x_1, double x, int direction) {
    // error checking
    if (direction == 0) {
      COMBIGRID_OUT_WRN("interpolation direciton parameter = 0. Returning 0.0!!!", __FILE__,
                        __LINE__)
      return 0.0;
    }

    if (x_0 >= x_1 || x < x_0 || x > x_1) {
      COMBIGRID_OUT_WRN(
          " Invalid interpolation parameters! Please reconsider your input and "
          "try again."
          "Returning 0.0!!!",
          __FILE__, __LINE__)
      return 0.0;
    }

    // now calculate the interpolation coefficient according to the linear
    // interpolation formula...
    double coef = 0.0;

    if (direction < 0) {
      coef = (x_1 - x) / (x_1 - x_0);
    } else {
      coef = (x - x_0) / (x_1 - x_0);
    }

    return coef;
  }

  /**
   * PeTz:
   *
   * this method finds the bounding box of the coordinate x, expressed as the
   *d-dimensional cube 4
   * [x_i1;x_{i1+1}] x [x_i2;x_{i2+1}] ... [x_id;x_{id+1}]
   * where i1, i2 ... id are indices of the current grid's points.
   *
   * The method returns {i1,i2 ... id} stored in the out-vector "bb_idx" as well
   *as the physical
   * coordinates {x_i1, x_i2 ... x_id} stored in the out-vector "bb_coord".
   *
   * @param k_idx - an int vector (of size dim) specifying the multi-index of
   *the sample point
   * @param new_levels - an int vector , also of size dim, specifying the new
   *grid's levels vector
   * @param bb_idx - (out param) an int vector (of size dim) containing the
   *computed bounding box indices.
   * @param bb_coord - (out param) a double vector containing the  2*d points'
   *coordinates specifying the bounding box of x!
   *
   *
   */

  void bounding_box(std::vector<int> k_idx, std::vector<int> new_levels, std::vector<int>& bb_idx,
                    std::vector<double>& bb_coord) const {
    if (k_idx.size() != dim_ || bb_idx.size() != dim_) {
      COMBIGRID_OUT_WRN(
          "Cannot find current point's bounding box! Invalid input "
          "parameters. ",
          __FILE__, __LINE__)
      return;
    }

    bb_coord.resize(2 * bb_idx.size(), 0.0);

    if (gridDomain_ == NULL) {
      /**
       *  No grid domain means that the grid consists of equally distributed
       *points in the interval [0;1]
       *
       */
      for (int d = 0; d < dim_; d++) {
        // calculate the interval width in the original level-l grid including
        // the boundaries !!!
        // here l = levels_[d];
        double dx = 1.0 / (powerOfTwo[levels_[d]]);
        // take into consideration if there are boundaries or not
        double offset = hasBoundaryPoints_[d] ? 0 : dx;

        bb_idx[d] = bb_left(k_idx[d], this->levels_[d], new_levels[d]);
        // and the interval's coordinates themselves
        bb_coord[2 * d] = bb_idx[d] * dx + offset;
        bb_coord[2 * d + 1] = (bb_idx[d] + 1) * dx + offset;
      }
    } else {
      /**
       *  Having grid domain means that we can extract the d_axis arrays from
       *each 1D domain separately...
       *
       */
      for (int d = 0; d < dim_; d++) {
        // as simple as ...
        std::vector<double> d_axis = gridDomain_->get1DDomain(d).axisScaling();
        bb_idx[d] = bb_left(k_idx[d], this->levels_[d], new_levels[d]);
        bb_coord[2 * d] = d_axis[bb_idx[d]];
        bb_coord[2 * d + 1] = d_axis[bb_idx[d] + 1];
      }
    }
  }

  /***
   * PeTz: A standard binary search algorithm searching through the ordered
   *double array "grid_line"
   * for the position of x in the ordering, i.e. effectively finds the index of
   *the smallest interval containing x.
   *
   * @param x - the point we are trying to find the surrounding interval of
   * @param grid_line - the grid axis we are looking in
   * @return the index of the largest element of grid_line, which is smaller or
   *equal to x.
   *
   * O( log(N)) operations ...
   *
   *
   */

  int binary_search(double x, std::vector<double> grid_line) const {
    if (grid_line.size() == 0) {
      COMBIGRID_OUT_ERR("Cannot search inside an empty array! Aborting!", __FILE__, __LINE__);
      return 0;
    }

    if (x < grid_line[0] || x > grid_line[grid_line.size() - 1]) {
      COMBIGRID_OUT_ERR("Input element outside array bounds! Aborting!", __FILE__, __LINE__);
      return 0;
    }

    int start = 0;
    int end = grid_line.size() - 1;
    int m = (start + end) / 2;

    while (!(x >= grid_line[m] && x <= grid_line[m + 1])) {
      if (x < grid_line[m]) {
        end = m;
      } else {
        start = m;
      }

      m = (start + end) / 2;
    }

    return m;
  }
};
}  // namespace combigrid

#endif /* COMBIFULLGRID_HPP_ */
