// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <sgpp/base/grid/storage/hashmap/HashGridPoint.hpp>
#include <sgpp/base/tools/ClenshawCurtisTable.hpp>

#include <sys/types.h>

#include <iostream>
#include <sstream>
#include <string>
#include <cmath>
#include <algorithm>
#include <utility>
#include <map>

namespace sgpp {
namespace base {

HashGridPoint::HashGridPoint(size_t dimension)
    : dimension(dimension), level(nullptr), index(nullptr), hInv(nullptr), hash(0) {
  level = new level_type[dimension];
  index = new index_type[dimension];
  hInv = new index_type[dimension];
  leaf = false;
}

HashGridPoint::HashGridPoint()
    : dimension(0), level(nullptr), index(nullptr), hInv(nullptr), hash(0) {
  leaf = false;
}

HashGridPoint::HashGridPoint(const HashGridPoint& o)
    : dimension(o.dimension), level(nullptr), index(nullptr), hInv(nullptr), hash(0) {
  level = new level_type[dimension];
  index = new index_type[dimension];
  hInv = new index_type[dimension];
  leaf = false;

  for (size_t d = 0; d < dimension; d++) {
    level[d] = o.level[d];
    index[d] = o.index[d];
  }

  leaf = o.leaf;
  rehash();
}

HashGridPoint::HashGridPoint(std::istream& istream, int version)
    : dimension(0), level(nullptr), index(nullptr), hInv(nullptr), hash(0) {
  size_t temp_leaf;

  istream >> dimension;

  level = new level_type[dimension];
  index = new index_type[dimension];
  hInv = new index_type[dimension];
  leaf = false;

  for (size_t d = 0; d < dimension; d++) {
    istream >> level[d];
    istream >> index[d];
  }

  if (version >= 2 && version != 4) {
    // read leaf option
    istream >> temp_leaf;

    if (temp_leaf == 0) {
      leaf = false;
    } else {
      leaf = true;
    }
  } else {
    leaf = false;
  }

  // PointDistribution was introduced in version 6 and removed in version 9;
  // this information is now stored in Stretching globally for all grid points
  if ((version >= 6) && (version <= 8)) {
    size_t temp_distr;
    istream >> temp_distr;
  }

  rehash();
}

/**
 * Destructor
 */
HashGridPoint::~HashGridPoint() {
  if (level) {
    delete[] level;
  }

  if (index) {
    delete[] index;
  }

  if (hInv) {
    delete[] hInv;
  }
}

void HashGridPoint::serialize(std::ostream& ostream, int version) {
  ostream << dimension << std::endl;

  for (size_t d = 0; d < dimension; d++) {
    ostream << level[d] << " ";
    ostream << index[d] << " ";
  }

  ostream << std::endl;

  ostream << leaf << std::endl;

  // PointDistribution was introduced in version 6 and removed in version 9;
  // this information is now stored in Stretching globally for all grid points
  if ((version >= 6) && (version <= 8)) {
    ostream << "Normal" << std::endl;
  }
}

size_t HashGridPoint::getDimension() const { return dimension; }

void HashGridPoint::setLeaf(bool isLeaf) { leaf = isLeaf; }

bool HashGridPoint::isLeaf() { return leaf; }

bool HashGridPoint::hasOverlappingSupport(size_t dim, HashGridPoint& gp) {
  size_t leveli = level[dim], indexi = index[dim];
  size_t levelj = gp.getLevel(dim), indexj = gp.getIndex(dim);

  if (leveli == levelj) return indexi == indexj;

  if (leveli < levelj)
    return isHierarchicalAncestor(dim, gp);
  else
    return gp.isHierarchicalAncestor(dim, *this);
}

bool HashGridPoint::isHierarchicalAncestor(size_t dim, HashGridPoint& gp) {
  size_t leveli = level[dim], indexi = index[dim];
  size_t levelj = gp.getLevel(dim), indexj = gp.getIndex(dim);

  return (levelj >= leveli) && (indexi == ((indexj >> (levelj - leveli)) | 1));
}

bool HashGridPoint::hasOverlappingSupport(HashGridPoint& gp) {
  size_t idim = 0;

  while (idim < dimension && hasOverlappingSupport(idim, gp)) ++idim;

  // check whether the supports are overlapping in all dimensions
  return idim == dimension;
}

bool HashGridPoint::isHierarchicalAncestor(HashGridPoint& gp) {
  size_t idim = 0;

  while (idim < dimension && isHierarchicalAncestor(idim, gp)) ++idim;

  // check whether the supports are overlapping in all dimensions
  return idim == dimension;
}

void HashGridPoint::getStandardCoordinates(DataVector& coordinates) const {
  coordinates.resize(dimension);

  for (size_t d = 0; d < dimension; d++) {
    coordinates.set(d, getStandardCoordinate(d));
  }
}

bool HashGridPoint::isInnerPoint() const {
  for (size_t d = 0; d < dimension; d++) {
    if (level[d] == 0) {
      return false;
    }
  }

  return true;
}

void HashGridPoint::rehash() {
  size_t hash = 0xdeadbeef;

  for (size_t d = 0; d < dimension; d++) {
    hInv[d] = static_cast<index_type>(1) << level[d];
    hash = hInv[d] + index[d] + hash * 65599;
  }

  this->hash = hash;
}

size_t HashGridPoint::getHash() const { return hash; }

bool HashGridPoint::equals(const HashGridPoint& rhs) const {
  for (size_t d = 0; d < dimension; d++) {
    if (level[d] != rhs.level[d]) {
      return false;
    }
  }

  for (size_t d = 0; d < dimension; d++) {
    if (index[d] != rhs.index[d]) {
      return false;
    }
  }

  return true;
}

HashGridPoint& HashGridPoint::assign(const HashGridPoint& rhs) { return this->operator=(rhs); }

HashGridPoint& HashGridPoint::operator=(const HashGridPoint& rhs) {
  if (this == &rhs) {
    return *this;
  }

  if (dimension != rhs.dimension) {
    if (level) {
      delete[] level;
    }

    if (index) {
      delete[] index;
    }

    if (hInv) {
      delete[] hInv;
    }

    dimension = rhs.dimension;

    level = new level_type[dimension];
    index = new index_type[dimension];
    hInv = new index_type[dimension];
  }

  for (size_t d = 0; d < dimension; d++) {
    level[d] = rhs.level[d];
    index[d] = rhs.index[d];
  }

  leaf = rhs.leaf;

  rehash();
  return *this;
}

std::string HashGridPoint::toString() const {
  std::ostringstream ostream;
  toString(ostream);

  return ostream.str();
}

void HashGridPoint::toString(std::ostream& stream) const {
  stream << "[";

  for (size_t i = 0; i < dimension; i++) {
    if (i != 0) {
      stream << ",";
    }

    stream << " " << this->level[i];
    stream << ", " << this->index[i];
  }

  stream << " ]";
}

HashGridPoint::level_type HashGridPoint::getLevelSum() const {
  HashGridPoint::level_type levelsum = 0;

  for (size_t d = 0; d < dimension; d++) {
    levelsum += level[d];
  }

  return levelsum;
}

HashGridPoint::level_type HashGridPoint::getLevelMax() const {
  HashGridPoint::level_type levelmax = level[0];

  for (size_t d = 1; d < dimension; d++) {
    levelmax = std::max(levelmax, level[d]);
  }

  return levelmax;
}

HashGridPoint::level_type HashGridPoint::getLevelMin() const {
  HashGridPoint::level_type levelmin = level[0];

  for (size_t d = 1; d < dimension; d++) {
    levelmin = std::min(levelmin, level[d]);
  }

  return levelmin;
}

}  // namespace base
}  // namespace sgpp
