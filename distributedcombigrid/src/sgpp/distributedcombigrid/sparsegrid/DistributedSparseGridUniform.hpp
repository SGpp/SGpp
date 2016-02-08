/*
 * DistributedSparseGrid.h
 *
 *  Created on: Oct 19, 2015
 *      Author: heenemo
 */

#ifndef SRC_SGPP_COMBIGRID_SPARSEGRID_DISTRIBUTEDSPARSEGRIDUNIFORM_HPP_
#define SRC_SGPP_COMBIGRID_SPARSEGRID_DISTRIBUTEDSPARSEGRIDUNIFORM_HPP_

#include <assert.h>

#include "sgpp/distributedcombigrid/utils/Types.hpp"

using namespace combigrid;

/*
 * Instead of having private static functions, I put these functions in an
 * unnamed namespace. So, they are not accessible from outside the file, as well.
 * In the general case, this would have the advantage, that we can change
 * the declaration of these functions without touching the declaration of the
 * class. So we avoid recompilation of all files that use the class.
 */
namespace {

template<typename FG_ELEMENT>
struct SubspaceSGU {
  LevelVector level_;

  IndexVector sizes_;

  size_t dataSize_;

  std::vector<FG_ELEMENT> data_;
};

/*
 template<typename FG_ELEMENT>
 bool mycomp( const SubspaceSGU<FG_ELEMENT>& ss1, const SubspaceSGU<FG_ELEMENT>& ss2 ){
 return ( ss1.dataSize_ > ss2.dataSize_ );
 }*/

} // end anonymous namespace

namespace combigrid {

template<typename FG_ELEMENT>
class DistributedSparseGridUniform {

 public:
  /** create sparse grid of dimension d and specify for each dimension the
   *  the maximum discretization level and whether there is a boundary or not
   */
  DistributedSparseGridUniform(DimType dim, const LevelVector& lmax,
                               const LevelVector& lmin, const std::vector<bool>& boundary,
                               CommunicatorType comm, size_t procsPerNode = 0);

  virtual ~DistributedSparseGridUniform();

  void print(std::ostream& os) const;

  // return level vector of subspace i
  inline const LevelVector& getLevelVector(size_t i) const;

  // return index of subspace i
  inline IndexType getIndex(const LevelVector& l) const;

  inline const std::vector<bool>& getBoundaryVector() const;

  // get pointer to first element in subspace with l
  inline FG_ELEMENT* getData(const LevelVector& l);

  // get pointer to first element in subspace i
  inline FG_ELEMENT* getData(size_t i);

  // get reference to data vector of subspace i.
  inline std::vector<FG_ELEMENT>& getDataVector(size_t i);

  // get reference to data vector of subspace with l.
  inline std::vector<FG_ELEMENT>& getDataVector(const LevelVector& l);

  inline size_t getDim() const;

  inline const LevelVector& getNMax() const;

  inline const LevelVector& getNMin() const;

  // return the number of subspaces
  inline size_t getNumSubspaces() const;

  // return the sizes for each dimension for subspace i
  inline const IndexVector& getSubspaceSizes(size_t i) const;

  // return the sizes for each dimension for subspace with l
  inline const IndexVector& getSubspaceSizes(const LevelVector& l) const;

  // return the number of elements of subspace i.
  // this number is independent of whether the subspace is initialized on this
  // process or not.
  inline size_t getSubspaceSize(size_t i) const;

  // return the number of elements of subspace i.
  // this number is independent of whether the subspace is initialized on this
  // process or not.
  inline size_t getSubspaceSize(const LevelVector& l) const;

  // check if a subspace with l is contained in the sparse grid
  // unlike getIndex this will not throw an assert in case l is not contained
  bool isContained(const LevelVector& l) const;

  inline size_t getDataSize(size_t i) const;

  inline size_t getDataSize(const LevelVector& l) const;

  inline CommunicatorType getCommunicator() const;

  inline int getCommunicatorSize() const;

 private:
  void createLevels(DimType dim, const LevelVector& nmax,
                    const LevelVector& lmin);

  void createLevelsRec(size_t dim, size_t n, size_t d, LevelVector& l,
                       const LevelVector& nmax);

  void setSizes();

  void calcProcAssignment(int procsPerNode);

  DimType dim_;

  LevelVector nmax_;

  LevelVector lmin_;

  std::vector<LevelVector> levels_;

  std::vector<bool> boundary_;

  CommunicatorType comm_;

  std::vector<RankType> subspaceToProc_;

  RankType rank_;

  int commSize_;

  std::vector<SubspaceSGU<FG_ELEMENT> > subspaces_;
};

} // namespace

namespace combigrid {

// at construction create only levels, no data
template<typename FG_ELEMENT>
DistributedSparseGridUniform<FG_ELEMENT>::DistributedSparseGridUniform(
  DimType dim, const LevelVector& lmax, const LevelVector& lmin,
  const std::vector<bool>& boundary, CommunicatorType comm,
  size_t procsPerNode) :
  dim_(dim) {
  assert(dim > 0);

  assert(lmax.size() == dim);

  for (size_t i = 0; i < lmax.size(); ++i)
    assert(lmax[i] > 0);

  assert(lmin.size() == dim);

  for (size_t i = 0; i < lmin.size(); ++i)
    assert(lmin[i] > 0);

  assert(boundary.size() == dim);

  MPI_Comm_rank(comm, &rank_);
  MPI_Comm_size(comm, &commSize_);
  comm_ = comm;

  nmax_ = lmax;
  lmin_ = lmin;
  boundary_ = boundary;

  std::cout << "nmax = " << nmax_ << std::endl;

  createLevels(dim, nmax_, lmin_);

  subspaces_.resize(levels_.size());

  for (size_t i = 0; i < levels_.size(); ++i)
    subspaces_[i].level_ = levels_[i];

  setSizes();

  /*
   // sort subspaces by datasize in descending order and update levels_ accordingly
   std::sort( subspaces_.begin(), subspaces_.end(), mycomp<FG_ELEMENT> );
   for( size_t i=0; i<levels_.size(); ++i )
   levels_[i] = subspaces_[i].level_;

   subspaceToProc_.resize( subspaces_.size() );
   calcProcAssignment( int(procsPerNode) );
   */
}

template<typename FG_ELEMENT>
void DistributedSparseGridUniform<FG_ELEMENT>::print(std::ostream& os) const {
  for (size_t i = 0; i < subspaces_.size(); ++i) {
    std::cout << i << " " << subspaces_[i].level_ << " " << subspaces_[i].sizes_
              << " " << subspaces_[i].dataSize_ << " " << "r " << subspaceToProc_[i]
              << std::endl;
  }
}

template<typename FG_ELEMENT>
std::ostream& operator<<(std::ostream& os,
                         const DistributedSparseGridUniform<FG_ELEMENT>& sg) {
  sg.print(os);
  return os;
}

template<typename FG_ELEMENT>
DistributedSparseGridUniform<FG_ELEMENT>::~DistributedSparseGridUniform() {
}

// start recursion by setting dim=d=dimensionality of the vector space
template<typename FG_ELEMENT>
void DistributedSparseGridUniform<FG_ELEMENT>::createLevelsRec(size_t dim,
    size_t n, size_t d, LevelVector& l, const LevelVector& nmax) {
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
void DistributedSparseGridUniform<FG_ELEMENT>::createLevels(DimType dim,
    const LevelVector& nmax, const LevelVector& lmin) {
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
void DistributedSparseGridUniform<FG_ELEMENT>::setSizes() {
  for (size_t i = 0; i < subspaces_.size(); ++i) {
    IndexVector& sizes = subspaces_[i].sizes_;
    sizes.resize(dim_);
    const LevelVector& l = subspaces_[i].level_;

    // loop over all dimensions
    for (size_t j = 0; j < dim_; ++j) {
      if (l[j] == 1 && boundary_[j]) {
        sizes[j] = (size_t(std::pow(2.0, real(l[j] - 1))) + size_t(2));
      } else {
        sizes[j] = size_t(std::pow(2.0, real(l[j] - 1)));
      }
    }

    IndexType tmp(1);

    for (auto s : sizes)
      tmp *= s;

    subspaces_[i].dataSize_ = size_t(tmp);
  }
}

template<typename FG_ELEMENT>
inline const LevelVector&
DistributedSparseGridUniform<FG_ELEMENT>::getLevelVector(size_t i) const {
  return levels_[i];
}

template<typename FG_ELEMENT>
void DistributedSparseGridUniform<FG_ELEMENT>::calcProcAssignment(
  int procsPerNode) {
  if (procsPerNode == 0) {
    for (size_t i = 0; i < subspaces_.size(); ++i)
      subspaceToProc_[i] = int(i) % commSize_;
  }
  /*
   else{
   // check if commsize a multiple of procs per node
   assert( (commSize_ % procsPerNode) == 0
   && "number of procs in comm must be multiple of procsPerNode" );
   int numNodes = commSize_ / procsPerNode;
   for( int i=0; i<int( subspaces_.size() ); ++i ){
   int nodeID = i % numNodes;
   int procInNodeID = ( i / numNodes ) % procsPerNode;
   subspaceToProc_[i] = nodeID * procsPerNode + procInNodeID;
   }
   }*/

  // use this to assign subspaces to first proc in group
  else {
    // todo: ceil
    int numNodes = static_cast<int>(std::ceil(
                                      real(commSize_) / real(procsPerNode)));

    for (int i = 0; i < int(subspaces_.size()); ++i) {
      int nodeID = i % numNodes;
      subspaceToProc_[i] = nodeID * procsPerNode;
    }
  }
}

/* get index of subspace with l. returns -1 if not included */
template<typename FG_ELEMENT>
IndexType DistributedSparseGridUniform<FG_ELEMENT>::getIndex(
  const LevelVector& l) const {
  // get index of l
  for (IndexType i = 0; i < IndexType(levels_.size()); ++i) {
    if (levels_[i] == l) {
      return i;
    }
  }

  return -1;
}

template<typename FG_ELEMENT>
inline const std::vector<bool>&
DistributedSparseGridUniform<FG_ELEMENT>::getBoundaryVector() const {
  return boundary_;
}

template<typename FG_ELEMENT>
inline FG_ELEMENT*
DistributedSparseGridUniform<FG_ELEMENT>::getData(const LevelVector& l) {
  IndexType i = getIndex(l);

  if (i < 0) {
    std::cout << "l = " << l << " not included in distributed sparse grid"
              << std::endl;
    assert(false);
  }

  return &subspaces_[i].data_[0];
}

template<typename FG_ELEMENT>
inline FG_ELEMENT*
DistributedSparseGridUniform<FG_ELEMENT>::getData(size_t i) {
  return &subspaces_[i].data_[0];
}

template<typename FG_ELEMENT>
inline std::vector<FG_ELEMENT>&
DistributedSparseGridUniform<FG_ELEMENT>::getDataVector(size_t i) {
  return subspaces_[i].data_;
}

template<typename FG_ELEMENT>
inline std::vector<FG_ELEMENT>&
DistributedSparseGridUniform<FG_ELEMENT>::getDataVector(const LevelVector& l) {
  IndexType i = getIndex(l);

  if (i < 0) {
    std::cout << "l = " << l << " not included in distributed sparse grid"
              << std::endl;
    assert(false);
  }

  return subspaces_[i].data_;
}

template<typename FG_ELEMENT>
inline DimType DistributedSparseGridUniform<FG_ELEMENT>::getDim() const {
  return dim_;
}

template<typename FG_ELEMENT>
inline const LevelVector& DistributedSparseGridUniform<FG_ELEMENT>::getNMax()
const {
  return nmax_;
}

template<typename FG_ELEMENT>
inline const LevelVector& DistributedSparseGridUniform<FG_ELEMENT>::getNMin()
const {
  return lmin_;
}

template<typename FG_ELEMENT>
inline size_t DistributedSparseGridUniform<FG_ELEMENT>::getNumSubspaces()
const {
  return subspaces_.size();
}

template<typename FG_ELEMENT>
bool DistributedSparseGridUniform<FG_ELEMENT>::isContained(
  const LevelVector& l) const {
  // get index of l
  bool found = false;

  for (size_t i = 0; i < levels_.size(); ++i) {
    if (levels_[i] == l) {
      found = true;
      break;
    }
  }

  return found;
}

template<typename FG_ELEMENT>
const IndexVector&
DistributedSparseGridUniform<FG_ELEMENT>::getSubspaceSizes(size_t i) const {
  return subspaces_[i].sizes_;
}

template<typename FG_ELEMENT>
const IndexVector&
DistributedSparseGridUniform<FG_ELEMENT>::getSubspaceSizes(
  const LevelVector& l) const {
  IndexType i = getIndex(l);

  if (i < 0) {
    std::cout << "l = " << l << " not included in distributed sparse grid"
              << std::endl;
    assert(false);
  }

  return subspaces_[i].sizes_;
}

template<typename FG_ELEMENT>
size_t DistributedSparseGridUniform<FG_ELEMENT>::getSubspaceSize(
  size_t i) const {
  return subspaces_[i].dataSize_;
}

template<typename FG_ELEMENT>
size_t DistributedSparseGridUniform<FG_ELEMENT>::getSubspaceSize(
  const LevelVector& l) const {
  IndexType i = getIndex(l);

  if (i < 0) {
    std::cout << "l = " << l << " not included in distributed sparse grid"
              << std::endl;
    assert(false);
  }

  return subspaces_[i].dataSize_;
}

template<typename FG_ELEMENT>
size_t DistributedSparseGridUniform<FG_ELEMENT>::getDataSize(size_t i) const {
  return subspaces_[i].data_.size();
}

template<typename FG_ELEMENT>
size_t DistributedSparseGridUniform<FG_ELEMENT>::getDataSize(
  const LevelVector& l) const {
  IndexType i = getIndex(l);

  if (i < 0) {
    std::cout << "l = " << l << " not included in distributed sparse grid"
              << std::endl;
    assert(false);
  }

  return subspaces_[i].data_.size();
}

template<typename FG_ELEMENT>
CommunicatorType DistributedSparseGridUniform<FG_ELEMENT>::getCommunicator()
const {
  return comm_;
}

template<typename FG_ELEMENT>
inline int DistributedSparseGridUniform<FG_ELEMENT>::getCommunicatorSize()
const {
  return commSize_;
}

} /* namespace combigrid */

#endif /* SRC_SGPP_COMBIGRID_SPARSEGRID_DISTRIBUTEDSPARSEGRID_HPP_ */
