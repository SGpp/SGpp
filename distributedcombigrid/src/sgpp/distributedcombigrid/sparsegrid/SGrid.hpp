/*
 * SGrid.hpp
 *
 *  Created on: May 14, 2013
 *      Author: heenemo
 */

#ifndef SGRID_HPP_
#define SGRID_HPP_

#include <assert.h>
#include <algorithm>
#include <cmath>
#include <cstdlib>      // std::rand, std::srand
#include <ctime>        // std::time
#include <iostream>
#include <vector>

#include "sgpp/distributedcombigrid/utils/IndexVector.hpp"
#include "sgpp/distributedcombigrid/utils/LevelVector.hpp"
#include "sgpp/distributedcombigrid/utils/Types.hpp"
#include "sgpp/distributedcombigrid/fullgrid/FullGrid.hpp"

namespace combigrid {

template<typename FG_ELEMENT>
class SGrid {

 public:
  /** create sparse grid of dimension d and level n with or without boundary
   *  in all dimensions
   */
  SGrid(DimType dim, LevelType n, bool boundary = false);

  /** create sparse grid of dimension d and specify for each dimension the
   *  the maximum discretization level and whether there is a boundary or not
   */
  SGrid(DimType dim, const LevelVector& nmax, const LevelVector& lmin,
        const std::vector<bool>& boundary);

  virtual ~SGrid();

  void print(std::ostream& os) const;

  void initHierarchicalSpace(const LevelVector& l, const FG_ELEMENT& val);

  void initAll(FG_ELEMENT val);

  void initCombiSpace(const LevelVector& l, const FG_ELEMENT& val);

  void add(const SGrid<FG_ELEMENT>& sg);

  void add(const FullGrid<FG_ELEMENT>& hfg, real coeff);

  // TODO: not sure whether this is a good idea, since it runs the risk of
  // accessing data which is not initialized
  // however, this is analogous to the std::vector::data of C++11
  // the advantage is that we are decoupled from the actual way how data is
  // organised.
  inline FG_ELEMENT* getData(LevelVector& l);

  inline FG_ELEMENT* getData(size_t i);

  inline std::vector<FG_ELEMENT>& getDataVector(size_t i);

  inline size_t getSize() const;

  inline size_t getDataSize(size_t i) const;

  int getInitList(std::vector<int>& list) const;

  size_t getLevelIndex(const LevelVector& l) const;

  inline bool isContained(const LevelVector& l, size_t& index) const;

  inline bool isContained(const LevelVector& l) const;

  inline size_t getDim() const;

  inline const LevelVector& getNMax() const;

  inline const LevelVector& getNMin() const;

  inline const LevelVector& getLevelVector(size_t i) const;

  inline void initHierarchicalSpace(size_t i, const FG_ELEMENT& val);

  inline size_t getCombinedDataSize() const;

  //todo: remove
  void permute();

  inline const std::vector<bool>& getBoundaryVector() const;

  inline const IndexVector& getSubspaceSizes(IndexType i) const;

 private:
  void createLevels(DimType dim, const LevelVector& nmax,
                    const LevelVector& lmin);

  void createLevelsRec(size_t dim, size_t n, size_t d, LevelVector& l,
                       const LevelVector& nmax);

  void setSizes();

  DimType dim_;

  LevelVector nmax_;

  LevelVector lmin_;

  std::vector<LevelVector> levels_;

  std::vector<std::vector<FG_ELEMENT> > data_;

  size_t combinedDataSize_;

  friend class CombiComTree;

  std::vector<bool> boundary_;

  std::vector<IndexVector> sizes_;
};

} // namespace

namespace combigrid {

template<typename FG_ELEMENT>
std::ostream& operator<<(std::ostream& os, const SGrid<FG_ELEMENT>& sg);

template<typename FG_ELEMENT>
inline void SGrid<FG_ELEMENT>::initHierarchicalSpace(size_t i,
    const FG_ELEMENT& val) {
  size_t num(1);
  const LevelVector& l = getLevelVector(i);

  // loop over all dimensions
  for (size_t j = 0; j < dim_; ++j) {
    if (l[j] == 1 && boundary_[j]) {
      num *= (size_t(std::pow(2.0, real(levels_[i][j] - 1))) + size_t(2));
    } else {
      num *= size_t(std::pow(2.0, real(levels_[i][j] - 1)));
    }
  }

  if (data_[i].size() == 0) {
    data_[i].resize(num, val);
    combinedDataSize_ += num;
  }
}

template<typename FG_ELEMENT>
inline FG_ELEMENT* SGrid<FG_ELEMENT>::getData(LevelVector& l) {
  return &data_[getLevelIndex(l)][0];
}

template<typename FG_ELEMENT>
inline FG_ELEMENT* SGrid<FG_ELEMENT>::getData(size_t i) {
  return &data_[i][0];
}

template<typename FG_ELEMENT>
std::vector<FG_ELEMENT>& SGrid<FG_ELEMENT>::getDataVector(size_t i) {
  return data_[i];
}

template<typename FG_ELEMENT>
inline size_t SGrid<FG_ELEMENT>::getSize() const {
  return levels_.size();
}

template<typename FG_ELEMENT>
inline size_t SGrid<FG_ELEMENT>::getDataSize(size_t i) const {
  return data_[i].size();
}

template<typename FG_ELEMENT>
inline DimType SGrid<FG_ELEMENT>::getDim() const {
  return dim_;
}

template<typename FG_ELEMENT>
inline const LevelVector& SGrid<FG_ELEMENT>::getNMax() const {
  return nmax_;
}

template<typename FG_ELEMENT>
inline const LevelVector& SGrid<FG_ELEMENT>::getNMin() const {
  return lmin_;
}

template<typename FG_ELEMENT>
inline const LevelVector& SGrid<FG_ELEMENT>::getLevelVector(size_t i) const {
  return levels_[i];
}

template<typename FG_ELEMENT>
inline size_t SGrid<FG_ELEMENT>::getCombinedDataSize() const {
  return combinedDataSize_;
}

// at construction create only levels, no data
template<typename FG_ELEMENT>
SGrid<FG_ELEMENT>::SGrid(DimType dim, LevelType n, bool boundary) :
  dim_(dim), combinedDataSize_(0) {
  assert(dim > 0);
  assert(n > 0);

  nmax_.resize(dim, n);
  lmin_.resize(dim, 1);
  boundary_.resize(dim, boundary);

  createLevels(dim, nmax_, lmin_);

  data_.resize(levels_.size());

  setSizes();
}

// at construction create only levels, no data
template<typename FG_ELEMENT>
SGrid<FG_ELEMENT>::SGrid(DimType dim, const LevelVector& nmax,
                         const LevelVector& lmin, const std::vector<bool>& boundary) :
  dim_(dim), combinedDataSize_(0) {
  assert(dim > 0);

  assert(nmax.size() == dim);

  for (size_t i = 0; i < nmax.size(); ++i)
    assert(nmax[i] > 0);

  assert(lmin.size() == dim);

  for (size_t i = 0; i < lmin.size(); ++i)
    assert(lmin[i] > 0);

  assert(boundary.size() == dim);

  nmax_ = nmax;
  lmin_ = lmin;
  boundary_ = boundary;

  createLevels(dim, nmax_, lmin_);

  data_.resize(levels_.size());

  setSizes();
}

template<typename FG_ELEMENT>
void SGrid<FG_ELEMENT>::print(std::ostream& os) const {
  for (size_t i = 0; i < levels_.size(); ++i) {
    for (size_t j = 0; j < levels_[i].size(); ++j)
      os << levels_[i][j] << " ";

    os << data_[i].size() << "\n";
  }
}

template<typename FG_ELEMENT>
std::ostream& operator<<(std::ostream& os, const SGrid<FG_ELEMENT>& sg) {
  sg.print(os);
  return os;
}

template<typename FG_ELEMENT>
void SGrid<FG_ELEMENT>::initHierarchicalSpace(const LevelVector& l,
    const FG_ELEMENT& val) {
  initHierarchicalSpace(getLevelIndex(l), val);
}

template<typename FG_ELEMENT>
void SGrid<FG_ELEMENT>::initAll(FG_ELEMENT val) {
  for (size_t i = 0; i < levels_.size(); ++i) {
    initHierarchicalSpace(i, val);
  }
}

/** initialize hierarchical spaces which are necessary for sparse grid
 *  representation of combi-fullgrid with level l
 *
 */
template<typename FG_ELEMENT>
void SGrid<FG_ELEMENT>::initCombiSpace(const LevelVector& l,
                                       const FG_ELEMENT& val) {
  assert(l.size() == dim_);

  // check if level vector contained in levels
  // implicitly done by getindex function
  getLevelIndex(l);

  for (size_t i = 0; i < levels_.size(); ++i) {
    LevelVector li(getLevelVector(i));

    if (li <= l)
      initHierarchicalSpace(i, val);
  }
}

/** add sparse grid sg to current sparse grid
 *
 */
template<typename FG_ELEMENT>
void SGrid<FG_ELEMENT>::add(const SGrid<FG_ELEMENT>& sg) {
  // todo: sanity check if sgrids are equal

  for (size_t i = 0; i < levels_.size(); ++i) {
    // case 1: subspace i not contained in both sparse grids
    if (this->data_[i].size() == 0 && sg.data_[i].size() == 0)
      continue;

    // subspace i only contained in "left" operand
    // nothing to do
    if (this->data_[i].size() != 0 && sg.data_[i].size() == 0)
      continue;

    // case 3: subspace i only contained in "right" operand
    // has to be created before addition
    if (this->data_[i].size() == 0) {
      initHierarchicalSpace(i, 0.0);
    }

    // add subspaces
    assert(data_[i].size() == sg.data_[i].size());

    for (size_t j = 0; j < data_[i].size(); ++j)
      data_[i][j] += sg.data_[i][j];
  }
}

/** add hierarchized full grid to sparse grid
 *
 * the hierarchized full grid must be completely contained in the space of
 * the sparse grid
 *
 */
template<typename FG_ELEMENT>
void SGrid<FG_ELEMENT>::add(const FullGrid<FG_ELEMENT>& hfg, real coeff) {
  assert(hfg.isGridCreated());
  assert(hfg.isHierarchized());

  // check if hfg contained in sg
  // in principle we dont need to do this and only copy the data of subspaces
  // which exists in both spaces
  LevelVector l_hfg(hfg.getLevels().begin(), hfg.getLevels().end());
  assert(isContained(l_hfg));

  // check if boundary flags equal
  assert(hfg.returnBoundaryFlags() == boundary_);

  // create iterator for each subspace of sgrid and set to begin
  typedef typename std::vector<FG_ELEMENT>::iterator SubspaceIterator;
  typename std::vector<SubspaceIterator> it_sub(getSize());

  for (size_t i = 0; i < it_sub.size(); ++i)
    it_sub[i] = data_[i].begin();

  for (size_t i = 0; i < hfg.getNrElements(); ++i) {
    // get level and index vector of element i
    LevelVector lvec(dim_);
    IndexVector ivec(dim_);
    hfg.getLI(i, lvec, ivec);

    // find index of sgrid space
    size_t j = getLevelIndex(lvec);

    // initialize subspace if necessary
    // this might invalidate the iterator, thus reassign iterator
    if (getDataSize(j) == 0) {
      initHierarchicalSpace(j, FG_ELEMENT(0));
      it_sub[j] = data_[j].begin();
    }

    assert(it_sub[j] != data_[j].end());

    // add hfg element to sgrid element
    *it_sub[j] += coeff * hfg.getData()[i];

    ++it_sub[j];
  }
}

template<typename FG_ELEMENT>
size_t SGrid<FG_ELEMENT>::getLevelIndex(const LevelVector& l) const {
  // get index of l
  bool found = false;
  size_t i;

  for (i = 0; i < levels_.size(); ++i) {
    if (levels_[i] == l) {
      found = true;
      break;
    }
  }

  if (!found)
    std::cout << l << " not included in sgrid" << std::endl;

  assert(found);

  return i;
}

template<typename FG_ELEMENT>
bool SGrid<FG_ELEMENT>::isContained(const LevelVector& l, size_t& index) const {
  // get index of l
  bool found = false;

  for (size_t i = 0; i < levels_.size(); ++i) {
    if (levels_[i] == l) {
      found = true;
      index = i;
      break;
    }
  }

  return found;
}

template<typename FG_ELEMENT>
bool SGrid<FG_ELEMENT>::isContained(const LevelVector& l) const {
  size_t tmp;
  return isContained(l, tmp);
}

template<typename FG_ELEMENT>
int SGrid<FG_ELEMENT>::getInitList(std::vector<int>& list) const {
  // TODO: check list size????
  int count(0);

  for (size_t i = 0; i < levels_.size(); ++i) {
    if (data_[i].size() > 0) {
      list[i] = 1;
      ++count;
    } else
      list[i] = 0;
  }

  return count;
}

template<typename FG_ELEMENT>
SGrid<FG_ELEMENT>::~SGrid() {
  // TODO Auto-generated destructor stub
}

// start recursion by setting dim=d=dimensionality of the vector space
template<typename FG_ELEMENT>
void SGrid<FG_ELEMENT>::createLevelsRec(size_t dim, size_t n, size_t d,
                                        LevelVector& l, const LevelVector& nmax) {
  // sum rightmost entries of level vector
  LevelType lsum(0);

  for (size_t i = dim; i < l.size(); ++i)
    lsum += l[i];

  for (LevelType ldim = 1; ldim <= LevelType(n) + LevelType(d) - 1 - lsum;
       ++ldim) {
    l[dim - 1] = ldim;

    if (dim == 1) {
      if (l <= nmax) {
        levels_.push_back(l);
        //std::cout << l << std::endl;
      }
    } else {
      createLevelsRec(dim - 1, n, d, l, nmax);
    }
  }
}

template<typename FG_ELEMENT>
void SGrid<FG_ELEMENT>::createLevels(DimType dim, const LevelVector& nmax,
                                     const LevelVector& lmin) {
  assert(nmax.size() == dim);
  assert(lmin.size() == dim);

  // compute c which fulfills nmax - c*1  >= lmin

  LevelVector ltmp(nmax);
  LevelType c = 0;

  while (ltmp > lmin) {
    ++c;

    for (size_t i = 0; i < dim; ++i) {
      ltmp[i] = nmax[i] - c;

      if (ltmp[i] < 1)
        ltmp[i] = 1;
    }
  }

  LevelVector rlmin(dim);

  for (size_t i = 0; i < rlmin.size(); ++i) {
    rlmin[i] = nmax[i] - c;
  }

  LevelType n = sum(rlmin) + c - dim + 1;

  LevelVector l(dim);
  createLevelsRec(dim, n, dim, l, nmax);
}

template<typename FG_ELEMENT>
void SGrid<FG_ELEMENT>::permute() {
  std::srand(unsigned(std::time(0)));
  std::random_shuffle(levels_.begin(), levels_.end());
}

template<typename FG_ELEMENT>
inline const std::vector<bool>& SGrid<FG_ELEMENT>::getBoundaryVector() const {
  return boundary_;
}

template<typename FG_ELEMENT>
void SGrid<FG_ELEMENT>::setSizes() {
  sizes_.resize(levels_.size(), IndexVector(dim_));

  for (size_t i = 0; i < levels_.size(); ++i) {
    IndexVector& sizes = sizes_[i];
    const LevelVector& l = getLevelVector(i);

    // loop over all dimensions
    for (size_t j = 0; j < dim_; ++j) {
      if (l[j] == 1 && boundary_[j]) {
        sizes[j] = (size_t(std::pow(2.0, real(l[j] - 1))) + size_t(2));
      } else {
        sizes[j] = size_t(std::pow(2.0, real(l[j] - 1)));
      }
    }
  }
}

template<typename FG_ELEMENT>
const IndexVector& SGrid<FG_ELEMENT>::getSubspaceSizes(IndexType i) const {
  return sizes_[i];
}

} /* namespace combigrid */

#endif /* SGRID_HPP_ */
