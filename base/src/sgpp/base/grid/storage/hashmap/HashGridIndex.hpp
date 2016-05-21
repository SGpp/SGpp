// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#ifndef HASHGRIDINDEX_HPP
#define HASHGRIDINDEX_HPP

#include <sgpp/base/datatypes/DataVector.hpp>

#include <sgpp/base/grid/common/BoundingBox.hpp>
#include <sgpp/base/grid/common/Stretching.hpp>

#include <sgpp/globaldef.hpp>

#include <sys/types.h>

#include <iostream>
#include <sstream>
#include <string>
#include <cmath>
#include <algorithm>
#include <map>

namespace sgpp {
namespace base {

/**
 * This Class represents one Gridpoint.
 *
 * A Gridpoint is given by its
 * ansatzfunctions that are not zero in every dimension. Instances
 * of this class are members in the hashmap that represents the
 * whole grid.
 */
class HashGridIndex {
 public:
  /// level type
  typedef uint32_t level_type;
  /// index type
  typedef uint32_t index_type;

  /// how the coordinates of the points are calculated
  enum class PointDistribution { Normal, ClenshawCurtis };

  /**
   * Constructor of a n-Dim gridpoint
   *
   * @param dimension the dimension of the gridpoint
   */
  explicit HashGridIndex(size_t dimension);

  /**
   * Standard-Constructor
   */
  HashGridIndex();

  /**
   * Copy-Constructor
   *
   * @param o constant reference to HashGridIndex object
   */
  HashGridIndex(const HashGridIndex& o);

  /**
   * Serialisation-Constructor
   *
   * @param istream instream object the contains the information about the gridpoint
   * @param version the serialization version of the file
   */
  HashGridIndex(std::istream& istream, int version);

  /**
   * Destructor
   */
  ~HashGridIndex();

  /**
   * Serialize this Gridpoint e.g. for a storage or checkpointing
   *
   * @param ostream outstream object to which the gridpoint's information is written
   * @param version the serialization version of the file
   */
  void serialize(std::ostream& ostream, int version);

  /**
   * Gets the dimension of the gridpoint
   *
   * @return the dimension of the gridpoint
   */
  size_t getDimension() const;

  /**
   * Sets level <i>l</i> and index <i>i</i> in dimension <i>d</i> and rehashs the HashGridIndex
   * object
   *
   * @param d the dimension in which the ansatzfunction is set
   * @param l the level of the ansatzfunction
   * @param i the index of the ansatzfunction
   */
  inline void set(size_t d, level_type l, index_type i) {
    level[d] = l;
    index[d] = i;
    rehash();
  }

  /**
   * Sets level <i>l</i> and index <i>i</i> in dimension <i>d</i> and the Leaf property and rehashs
   * the HashGridIndex object
   *
   * @param d the dimension in which the ansatzfunction is set
   * @param l the level of the ansatzfunction
   * @param i the index of the ansatzfunction
   * @param isLeaf specifies if this gridpoint has any childrens in any dimension
   */
  inline void set(size_t d, level_type l, index_type i, bool isLeaf) {
    level[d] = l;
    index[d] = i;
    leaf = isLeaf;
    rehash();
  }

  /**
   * Sets level <i>l</i> and index <i>i</i> in dimension <i>d</i> and doesn't rehash the
   * HashGridIndex object
   *
   * @param d the dimension in which the ansatzfunction is set
   * @param l the level of the ansatzfunction
   * @param i the index of the ansatzfunction
   */
  inline void push(size_t d, level_type l, index_type i) {
    level[d] = l;
    index[d] = i;
  }

  /**
   * Sets level <i>l</i> and index <i>i</i> in dimension <i>d</i> and the Leaf property and doesn't
   * rehash the HashGridIndex object
   *
   * @param d the dimension in which the ansatzfunction is set
   * @param l the level of the ansatzfunction
   * @param i the index of the ansatzfunction
   * @param isLeaf specifies if this gridpoint has any childrens in any dimension
   */
  inline void push(size_t d, level_type l, index_type i, bool isLeaf) {
    level[d] = l;
    index[d] = i;
    leaf = isLeaf;
  }

  /**
   * gets level <i>l</i> and index <i>i</i> in dimension <i>d</i> by reference parameters
   *
   * @param d the dimension in which the ansatz function should be read
   * @param l reference parameter for the level of the ansatz function
   * @param i reference parameter for the index of the ansatz function
   */
  inline void get(size_t d, level_type& l, index_type& i) const {
    l = level[d];
    i = index[d];
  }

  /**
   * gets level <i>l</i> in dimension <i>d</i>
   *
   * @param d the dimension in which the ansatz function should be read
   * @return level
   */
  inline level_type getLevel(size_t d) const { return level[d]; }

  /**
   * gets index <i>i</i> in dimension <i>d</i>
   *
   * @param d the dimension in which the ansatz function should be read
   * @return index
   */
  inline index_type getIndex(size_t d) const { return index[d]; }

  /**
   * Gets the point distribution of the grid point.
   *
   * @return point distribution
   */
  PointDistribution getPointDistribution() const;

  /**
   * Sets the point distribution of the grid point.
   *
   * @param distr new point distribution
   */
  void setPointDistribution(PointDistribution distr);

  /**
   * Set the leaf property; a grid point is called a leaf, if it has <b>not a single</b> child.
   *
   * @param isLeaf specifies if the current index is a leaf (i.e. has <b>no</b> child nodes) or not
   */
  void setLeaf(bool isLeaf);

  /**
   * Checks if this grid point has <b>not a single</b> child in any dimension.
   *
   * @return Returns true if this grid point has <b>no</b> children, otherwise false
   */
  bool isLeaf();

  /**
   * determines the coordinate in a given dimension
   *
   * @param d the dimension in which the coordinate should be calculated
   *
   * @return the coordinate in the given dimension
   */
  double getCoord(size_t d) const;

  /**
   * determines the coordinate in a given dimension
   *
   * @param d the dimension in which the coordinate should be calculated
   * @param q the intervals width in this dimension
   * @param t the offset in this dimension
   *
   * @return the coordinate in the given dimension
   */
  double getCoordBB(size_t d, double q, double t) const;

  /**
   * determines the coordinate in a given dimension
   *
   * @param d the dimension in which the coordinate should be calculated
   * @param stretch the stretching the index belongs to
   *
   * @return the coordinate in the given dimension
   */
  double getCoordStretching(size_t d, Stretching* stretch);

  /**
   * determines if the grid point is an inner grid point
   *
   * @return true if the grid point is an inner grid point
   */
  bool isInnerPoint() const;

  /**
   * rehashs the current gridpoint and sets hInv
   */
  void rehash();

  /**
   * gets the hash value of the current instance
   *
   * @return the hash value of the instance
   */
  size_t getHash() const;

  /**
   * checks whether this gridpoints is identical to another one
   *
   * Running under WINDOW this is defined the way around, MSDN 2009:
   * A binary predicate f(x,y) is a function object that has two
   * argument objects x and y and a return value of true or false.
   * An ordering imposed on a hash_map is a strict weak ordering
   * if the binary predicate is irreflexive, antisymmetric,
   * and transitive and if equivalence is transitive, where
   * two objects x and y are defined to be equivalent
   * when both f(x,y) and f(y,x) are false -> equalsSGLRBHash
   *
   * @param rhs reference the another HashGridIndex instance
   *
   * @return true if the gridpoints are identical otherwise false
   */
  bool equals(const HashGridIndex& rhs) const;

  /**
   * A wrapper for operator=
   *
   * @param rhs a reference to a HashGridIndex that contains the values that should be copied
   *
   * @return returns a reference HashGridIndex
   */
  HashGridIndex& assign(const HashGridIndex& rhs);

  /**
   * operator to assign the current grid point with the values of another one
   *
   * @param rhs a reference to a HashGridIndex that contains the values that should be copied
   *
   * @return returns a reference HashGridIndex
   */
  HashGridIndex& operator=(const HashGridIndex& rhs);

  /**
   * Generates a string with level and index of the gridpoint.
   * The format is <tt>[l1, i1, l2, i2, ..., ld, id]</tt>.
   * Needed for Java compatibility.
   *
   * @returns string into which the gridpoint is written
   */
  std::string toString() const;

  /**
   * Generates a string with level and index of the gridpoint.
   * The format is <tt>[l1, i1, l2, i2, ..., ld, id]</tt>.
   *
   * @param stream reference to a output stream
   */
  void toString(std::ostream& stream) const;

  /**
   * Sets the entries of DataVector p to the coordinates of the gridpoint
   *
   * @param p the (result) DataVector p that should be overwritten
   */
  void getCoords(DataVector& p) const;

  /**
   * Sets the entries of DataVector p to the coordinates of the gridpoint with bounding box
   *
   * @param p the (result) DataVector p that should be overwritten
   * @param BB reference to BoundingBox Object, that stores all boundaries for all dimensions
   */
  void getCoordsBB(DataVector& p, BoundingBox& BB) const;

  /**
   * Sets the entries of DataVector p to the coordinates of the gridpoint with stretching
   *
   * @param p the (result) DataVector p that should be overwritten
   * @param stretch reference to Stretching Object, that stores grid points in all dimensions
   */
  void getCoordsStretching(DataVector& p, Stretching& stretch) const;

  /**
   * Generates a string with all coordinates of the grid point.
   * The accuracy is up to 6 digits, i.e. beginning with level 8 there are rounding errors.
   *
   * @return returns a string with the coordinates of the grid point separated by whitespace
   */
  std::string getCoordsString() const;

  /**
   * Generates a string with all coordinates of the grid point with bounding box
   * The accuracy is up to 6 digits, i.e. beginning with level 8 there are rounding errors.
   *
   * This version scales the coordinates with q and t
   *
   * @param BB reference to BoundingBox Object, that stores all boundaries for all dimensionst
   *
   * @return returns a string with the coordinates of the grid point separated by whitespace
   */
  std::string getCoordsStringBB(BoundingBox& BB) const;

  /**
   * Generates a string with all coordinates of the grid point with bounding box
   * The accuracy is up to 6 digits, i.e. beginning with level 8 there are rounding errors.
   *
   * This version scales the coordinates with q and t
   *
   * @param stretch reference to Stretching Object, that stores all boundaries for all dimensions
   *
   * @return returns a string with the coordinates of the grid point separated by whitespace
   */
  std::string getCoordsStringStretching(Stretching& stretch) const;

  /**
   * Returns the sum of the one-dimensional levels, i.e., @f$ |\vec{l}|_1 @f$.
   *
   * @return the sum of the one-dimensional levels
   */
  level_type getLevelSum() const;

  /**
   * Returns the maximum of the one-dimensional levels, i.e., @f$ |\vec{l}|_\infty @f$.
   *
   * @return the maximum of the one-dimensional levels
   */
  level_type getLevelMax() const;

  /**
   * Returns the minimum of the one-dimensional levels.
   *
   * @return the minimum of the one-dimensional levels
   */
  level_type getLevelMin() const;

  /**
   * Sets the index to the left level zero parent
   *
   * @param dim the dimension in which the modification is taken place
   */
  inline void getLeftLevelZero(size_t dim) {
    set(dim, 0, 0);
  }

  /**
   * Sets the index to the right level zero parent
   *
   * @param dim the dimension in which the modification is taken place
   */
  inline void getRightLevelZero(size_t dim) {
    set(dim, 0, 1);
  }

  /**
   * Sets the index to the left child
   *
   * @param dim the dimension in which the modification is taken place
   */
  inline void getLeftChild(size_t dim) {
    level_type l;
    index_type i;
    get(dim, l, i);
    set(dim, l + 1, 2 * i - 1);
  }

  /**
   * Sets the index to the right child
   *
   * @param dim the dimension in which the modification is taken place
   */
  inline void getRightChild(size_t dim) {
    level_type l;
    index_type i;
    get(dim, l, i);
    set(dim, l + 1, 2 * i + 1);
  }

  /**
   * Resets the index to the top level in direction d
   *
   * @param d the dimension in which the modification is taken place
   */
  inline void getRoot(size_t d) {
    set(d, 1, 1);
  }

 protected:
  typedef std::map<std::string, PointDistribution> pointDistributionMap;
  typedef std::map<PointDistribution, std::string> pointDistributionVerboseMap;

 private:
  /// the dimension of the gridpoint
  size_t dimension;
  /// pointer to array that stores the ansatzfunctions' level
  level_type* level;
  /// pointer to array that stores the ansatzfunctions' indices
  index_type* index;
  /// pointer to array that stores the mesh widths (1 << level[d] for each dimension)
  index_type* hInv;
  /// distribution of the grid point (Normal or ClenshawCurtis)
  PointDistribution distr;
  /// stores if this gridpoint is a leaf
  bool leaf;
  /// stores the hashvalue of the gridpoint
  size_t hash;

  static pointDistributionMap& typeMap();
  static pointDistributionVerboseMap& typeVerboseMap();

  friend struct HashGridIndexPointerHashFunctor;
  friend struct HashGridIndexPointerEqualityFunctor;
  friend struct HashGridIndexHashFunctor;
  friend struct HashGridIndexEqualityFunctor;
};

struct HashGridIndexPointerHashFunctor {
  size_t operator()(const HashGridIndex* index) const { return index->getHash(); }
};

struct HashGridIndexPointerEqualityFunctor {
  size_t operator()(const HashGridIndex* s1, const HashGridIndex* s2) const {
    return s1->equals(*s2);
  }
};

struct HashGridIndexHashFunctor {
  size_t operator()(const HashGridIndex& index) const { return index.getHash(); }
};

struct HashGridIndexEqualityFunctor {
  size_t operator()(const HashGridIndex& s1, const HashGridIndex& s2) const {
    return s1.equals(s2);
  }
};

}  // namespace base
}  // namespace sgpp

#endif /* HASHGRIDINDEX_HPP */
