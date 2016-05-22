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

HashGridPoint::HashGridPoint(size_t dimension) :
    dimension(dimension), level(nullptr), index(nullptr), hInv(nullptr),
    distr(PointDistribution::Normal), hash(0) {
  level = new level_type[dimension];
  index = new index_type[dimension];
  hInv = new index_type[dimension];
  leaf = false;
}

HashGridPoint::HashGridPoint() :
    dimension(0), level(nullptr), index(nullptr), hInv(nullptr),
    distr(PointDistribution::Normal), hash(0) {
  leaf = false;
}

HashGridPoint::HashGridPoint(const HashGridPoint& o) :
    dimension(o.dimension), level(nullptr), index(nullptr), hInv(nullptr),
    distr(PointDistribution::Normal), hash(0) {
  level = new level_type[dimension];
  index = new index_type[dimension];
  hInv = new index_type[dimension];
  leaf = false;

  for (size_t d = 0; d < dimension; d++) {
    level[d] = o.level[d];
    index[d] = o.index[d];
  }

  distr = o.distr;
  leaf = o.leaf;
  rehash();
}

HashGridPoint::HashGridPoint(std::istream& istream, int version) :
    dimension(0), level(nullptr), index(nullptr), hInv(nullptr), hash(0) {
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
  if (version > 5) {
    ostream << typeVerboseMap()[distr] << std::endl;
  }
}

size_t HashGridPoint::getDimension() const {
  return dimension;
}

HashGridPoint::PointDistribution HashGridPoint::getPointDistribution() const {
  return distr;
}

void HashGridPoint::setPointDistribution(HashGridPoint::PointDistribution distr) {
  this->distr = distr;
}

void HashGridPoint::setLeaf(bool isLeaf) {
  leaf = isLeaf;
}

bool HashGridPoint::isLeaf() {
  return leaf;
}

double HashGridPoint::getCoord(size_t d) const {
  if (distr == PointDistribution::Normal) {
    // cast 1 to index_type to ensure that 1 << level[d] doesn't overflow
    return static_cast<double>(index[d]) / static_cast<double>(hInv[d]);
  } else {
    return ClenshawCurtisTable::getInstance().getPoint(level[d], index[d]);
  }
}

double HashGridPoint::getCoordBB(size_t d, double q, double t) const {
  return q * getCoord(d) + t;
}

double HashGridPoint::getCoordStretching(size_t d, Stretching* stretch) {
  return stretch->getCoordinates(level[d], index[d], d);
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

size_t HashGridPoint::getHash() const {
  return hash;
}

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

HashGridPoint& HashGridPoint::assign(const HashGridPoint& rhs) {
  return this->operator=(rhs);
}

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

  distr = rhs.distr;
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

void HashGridPoint::getCoords(DataVector& p) const {
  for (size_t d = 0; d < dimension; d++) {
    p.set(d, getCoord(d));
  }
}

void HashGridPoint::getCoordsBB(DataVector& p, BoundingBox& BB) const {
  for (size_t d = 0; d < dimension; d++) {
    p.set(d, BB.getIntervalWidth(d) * getCoord(d) + BB.getIntervalOffset(d));
  }
}

void HashGridPoint::getCoordsStretching(DataVector& p, Stretching& stretch) const {
  for (size_t d = 0; d < dimension; d++) {
    if (level[d] == 0) {
      p.set(d, stretch.getIntervalWidth(d) * static_cast<double>(index[d]) +
                   stretch.getIntervalOffset(d));
    } else {
      p.set(d, stretch.getCoordinates(level[d], index[d], d));
    }
  }
}

std::string HashGridPoint::getCoordsString() const {
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

std::string HashGridPoint::getCoordsStringBB(BoundingBox& BB) const {
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

std::string HashGridPoint::getCoordsStringStretching(Stretching& stretch) const {
  std::stringstream return_stream;

  for (size_t d = 0; d < dimension; d++) {
    return_stream << std::scientific << stretch.getCoordinates(level[d], index[d], d);

    if (d < dimension - 1) {
      return_stream << " ";
    }
  }

  return return_stream.str();
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

std::map<std::string, base::HashGridPoint::PointDistribution>& HashGridPoint::typeMap() {
  // This is only executed once!
  static pointDistributionMap* tMap = new pointDistributionMap();

  if (tMap->size() == 0) {
/*
 * Insert strings here.
 */
#ifdef _WIN32
    tMap->insert(std::pair<std::string, base::HashGridPoint::PointDistribution>(
        "Normal", HashGridPoint::PointDistribution::Normal));
    tMap->insert(std::pair<std::string, base::HashGridPoint::PointDistribution>(
        "ClenshawCurtis", HashGridPoint::PointDistribution::ClenshawCurtis));
#else
    tMap->insert(std::make_pair("Normal", HashGridPoint::PointDistribution::Normal));
    tMap->insert(
        std::make_pair("ClenshawCurtis", HashGridPoint::PointDistribution::ClenshawCurtis));
#endif
  }

  return *tMap;
}

std::map<base::HashGridPoint::PointDistribution, std::string>& HashGridPoint::typeVerboseMap() {
  // This is only executed once!
  static pointDistributionVerboseMap* verboseMap = new pointDistributionVerboseMap();

  if (verboseMap->size() == 0) {
/*
 * Insert strings here.
 */
#ifdef _WIN32
    verboseMap->insert(std::pair<base::HashGridPoint::PointDistribution, std::string>(
        HashGridPoint::PointDistribution::Normal, "Normal"));
    verboseMap->insert(std::pair<base::HashGridPoint::PointDistribution, std::string>(
        HashGridPoint::PointDistribution::ClenshawCurtis, "ClenshawCurtis"));
#else
    verboseMap->insert(std::make_pair(HashGridPoint::PointDistribution::Normal, "Normal"));
    verboseMap->insert(
        std::make_pair(HashGridPoint::PointDistribution::ClenshawCurtis, "ClenshawCurtis"));
#endif
  }

  return *verboseMap;
}

}  // namespace base
}  // namespace sgpp
