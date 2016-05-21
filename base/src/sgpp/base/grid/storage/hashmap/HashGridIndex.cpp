// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <sgpp/base/grid/storage/hashmap/HashGridIndex.hpp>
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

HashGridIndex::HashGridIndex(size_t dimension)
    : dimension(dimension), level(NULL), index(NULL), distr(PointDistribution::Normal),
      hash_value(0) {
  level = new level_type[dimension];
  index = new index_type[dimension];
  Leaf = false;dimension
}

HashGridIndex::HashGridIndex()
    : dimension(0), level(NULL), index(NULL), distr(PointDistribution::Normal), hash_value(0) {
  Leaf = false;
}

HashGridIndex::HashGridIndex(const HashGridIndex& o)
    : dimension(o.dimension), level(NULL), index(NULL), distr(PointDistribution::Normal),
      hash_value(0) {
  level = new level_type[dimension];
  index = new index_type[dimension];
  Leaf = false;

  for (size_t d = 0; d < dimension; d++) {
    level[d] = o.level[d];
    index[d] = o.index[d];
  }

  distr = o.distr;
  Leaf = o.Leaf;
  rehash();
}

HashGridIndex::HashGridIndex(std::istream& istream, int version)
    : dimension(0), level(NULL), index(NULL), hash_value(0) {
  size_t temp_leaf;

  istream >> dimension;

  level = new level_type[dimension];
  index = new index_type[dimension];
  Leaf = false;

  for (size_t d = 0; d < dimension; d++) {
    istream >> level[d];
    istream >> index[d];
  }

  if (version >= 2 && version != 4) {
    // read leaf option
    istream >> temp_leaf;

    if (temp_leaf == 0) {
      Leaf = false;
    } else {
      Leaf = true;
    }
  } else {
    Leaf = false;
  }

  if (version == 6) {
    size_t temp_distr;
    istream >> temp_distr;
    distr = static_cast<PointDistribution>(temp_distr);
  } else if (version >= 7) {
    std::string temp_distr;
    istream >> temp_distr;

    if (typeMap().count(temp_distr) > 0) {
      distr = typeMap()[temp_distr];
    } else {
      distr = PointDistribution::Normal;
    }
  } else {
    distr = PointDistribution::Normal;
  }

  rehash();
}

/**
 * Destructor
 */
HashGridIndex::~HashGridIndex() {
  if (level) {
    delete[] level;
  }

  if (index) {
    delete[] index;
  }
}

void HashGridIndex::serialize(std::ostream& ostream, int version) {
  ostream << dimension << std::endl;

  for (size_t d = 0; d < dimension; d++) {
    ostream << level[d] << " ";
    ostream << index[d] << " ";
  }

  ostream << std::endl;

  ostream << Leaf << std::endl;
  if (version > 5) {
    ostream << typeVerboseMap()[distr] << std::endl;
  }
}

size_t HashGridIndex::getDimension() const { return dimension; }

HashGridIndex::PointDistribution HashGridIndex::getPointDistribution() const { return distr; }

void HashGridIndex::setPointDistribution(HashGridIndex::PointDistribution distr) {
  this->distr = distr;
}

void HashGridIndex::setLeaf(bool isLeaf) { Leaf = isLeaf; }

bool HashGridIndex::isLeaf() { return Leaf; }

double HashGridIndex::getCoord(size_t d) const {
  if (distr == PointDistribution::Normal) {
    // cast 1 to index_type to ensure that 1 << level[d] doesn't overflow
    return static_cast<double>(index[d]) /
           static_cast<double>(static_cast<index_type>(1) << level[d]);
  } else {
    return ClenshawCurtisTable::getInstance().getPoint(level[d], index[d]);
  }
}

double HashGridIndex::getCoordBB(size_t d, double q, double t) const { return q * getCoord(d) + t; }

double HashGridIndex::getCoordStretching(size_t d, Stretching* stretch) {
  return stretch->getCoordinates(level[d], index[d], d);
}

bool HashGridIndex::isInnerPoint() const {
  for (size_t d = 0; d < dimension; d++) {
    if (level[d] == 0) return false;
  }

  return true;
}

HashGridIndex* HashGridIndex::getPointer() { return this; }

void HashGridIndex::rehash() {
  size_t hash = 0xdeadbeef;

  for (size_t d = 0; d < dimension; d++) {
    hash = (1 << level[d]) + index[d] + hash * 65599;
  }

  hash_value = hash;
}

size_t HashGridIndex::hash() const { return hash_value; }

bool HashGridIndex::equals(const HashGridIndex& rhs) const {
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

HashGridIndex& HashGridIndex::assign(const HashGridIndex& rhs) { return this->operator=(rhs); }

HashGridIndex& HashGridIndex::operator=(const HashGridIndex& rhs) {
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

    dimension = rhs.dimension;

    level = new level_type[dimension];
    index = new index_type[dimension];
  }

  for (size_t d = 0; d < dimension; d++) {
    level[d] = rhs.level[d];
    index[d] = rhs.index[d];
  }

  distr = rhs.distr;
  Leaf = rhs.Leaf;

  rehash();
  return *this;
}

std::string HashGridIndex::toString() const {
  std::ostringstream ostream;
  toString(ostream);

  return ostream.str();
}

void HashGridIndex::toString(std::ostream& stream) const {
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

void HashGridIndex::getCoords(DataVector& p) const {
  for (size_t d = 0; d < dimension; d++) {
    p.set(d, getCoord(d));
  }
}

void HashGridIndex::getCoordsBB(DataVector& p, BoundingBox& BB) const {
  for (size_t d = 0; d < dimension; d++) {
    p.set(d, BB.getIntervalWidth(d) * getCoord(d) + BB.getIntervalOffset(d));
  }
}

void HashGridIndex::getCoordsStretching(DataVector& p, Stretching& stretch) const {
  for (size_t d = 0; d < dimension; d++) {
    if (level[d] == 0) {
      p.set(d, stretch.getIntervalWidth(d) * static_cast<double>(index[d]) +
                   stretch.getIntervalOffset(d));
    } else {
      p.set(d, stretch.getCoordinates(level[d], index[d], d));
    }
  }
}

std::string HashGridIndex::getCoordsString() const {
  std::stringstream return_stream;

  // switch on scientific notation:
  // return_stream << std::scientific;

  for (size_t d = 0; d < dimension; d++) {
    if (level[d] == 0) {
      return_stream << index[d];
    } else {
      return_stream << std::scientific << getCoord(d);
    }

    if (d < dimension - 1) {
      return_stream << " ";
    }
  }

  return return_stream.str();
}

std::string HashGridIndex::getCoordsStringBB(BoundingBox& BB) const {
  std::stringstream return_stream;

  for (size_t d = 0; d < dimension; d++) {
    return_stream << std::scientific
                  << BB.getIntervalWidth(d) * getCoord(d) + BB.getIntervalOffset(d);

    if (d < dimension - 1) {
      return_stream << " ";
    }
  }

  return return_stream.str();
}

std::string HashGridIndex::getCoordsStringStretching(Stretching& stretch) const {
  std::stringstream return_stream;

  for (size_t d = 0; d < dimension; d++) {
    return_stream << std::scientific << stretch.getCoordinates(level[d], index[d], d);

    if (d < dimension - 1) {
      return_stream << " ";
    }
  }

  return return_stream.str();
}

HashGridIndex::level_type HashGridIndex::getLevelSum() const {
  HashGridIndex::level_type levelsum = 0;

  for (size_t d = 0; d < dimension; d++) {
    levelsum += level[d];
  }

  return levelsum;
}

HashGridIndex::level_type HashGridIndex::getLevelMax() const {
  HashGridIndex::level_type levelmax = level[0];

  for (size_t d = 1; d < dimension; d++) {
    levelmax = std::max(levelmax, level[d]);
  }

  return levelmax;
}

HashGridIndex::level_type HashGridIndex::getLevelMin() const {
  HashGridIndex::level_type levelmin = level[0];

  for (size_t d = 1; d < dimension; d++) {
    levelmin = std::min(levelmin, level[d]);
  }

  return levelmin;
}

std::map<std::string, base::HashGridIndex::PointDistribution>& HashGridIndex::typeMap() {
  // This is only executed once!
  static pointDistributionMap* tMap = new pointDistributionMap();

  if (tMap->size() == 0) {
/*
 * Insert strings here.
 */
#ifdef _WIN32
    tMap->insert(std::pair<std::string, base::HashGridIndex::PointDistribution>(
        "Normal", HashGridIndex::PointDistribution::Normal));
    tMap->insert(std::pair<std::string, base::HashGridIndex::PointDistribution>(
        "ClenshawCurtis", HashGridIndex::PointDistribution::ClenshawCurtis));
#else
    tMap->insert(std::make_pair("Normal", HashGridIndex::PointDistribution::Normal));
    tMap->insert(
        std::make_pair("ClenshawCurtis", HashGridIndex::PointDistribution::ClenshawCurtis));
#endif
  }

  return *tMap;
}

std::map<base::HashGridIndex::PointDistribution, std::string>& HashGridIndex::typeVerboseMap() {
  // This is only executed once!
  static pointDistributionVerboseMap* verboseMap = new pointDistributionVerboseMap();

  if (verboseMap->size() == 0) {
/*
 * Insert strings here.
 */
#ifdef _WIN32
    verboseMap->insert(std::pair<base::HashGridIndex::PointDistribution, std::string>(
        HashGridIndex::PointDistribution::Normal, "Normal"));
    verboseMap->insert(std::pair<base::HashGridIndex::PointDistribution, std::string>(
        HashGridIndex::PointDistribution::ClenshawCurtis, "ClenshawCurtis"));
#else
    verboseMap->insert(std::make_pair(HashGridIndex::PointDistribution::Normal, "Normal"));
    verboseMap->insert(
        std::make_pair(HashGridIndex::PointDistribution::ClenshawCurtis, "ClenshawCurtis"));
#endif
  }

  return *verboseMap;
}

}  // namespace base
}  // namespace sgpp
