/* ****************************************************************************
 * Copyright (C) 2009 Technische Universitaet Muenchen                         *
 * This file is part of the SG++ project. For conditions of distribution and   *
 * use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
 **************************************************************************** */
// @author Dirk Pflueger (dirk.pflueger@in.tum.de), Jörg Blank (blankj@in.tum.de), Alexander Heinecke (Alexander.Heinecke@mytum.de)

#include <iostream>
#include <sstream>
#include <string>
#include <sys/types.h>
#include <cmath>
#include <algorithm>

#include <sgpp/base/grid/storage/hashmap/HashGridIndex.hpp>


namespace SGPP {
namespace base {


  HashGridIndex::HashGridIndex(size_t dim) :
      DIM(dim), level(NULL), index(NULL), hash_value(0) {
    level = new level_type[dim];
    index = new index_type[dim];
    Leaf = false;
  }

  HashGridIndex::HashGridIndex() :
      DIM(0), level(NULL), index(NULL), hash_value(0) {
    Leaf = false;
  }

  HashGridIndex::HashGridIndex(const HashGridIndex* o) :
      DIM(o->DIM), level(NULL), index(NULL), hash_value(0) {
    level = new level_type[DIM];
    index = new index_type[DIM];
    Leaf = false;

    for (size_t d = 0; d < DIM; d++) {
      level[d] = o->level[d];
      index[d] = o->index[d];
    }

    Leaf = o->Leaf;
    rehash();
  }

  HashGridIndex::HashGridIndex(std::istream& istream, int version) :
      DIM(0), level(NULL), index(NULL), hash_value(0) {
    size_t temp_leaf;

    istream >> DIM;

    level = new level_type[DIM];
    index = new index_type[DIM];
    Leaf = false;

    for (size_t d = 0; d < DIM; d++) {
      istream >> level[d];
      istream >> index[d];
    };

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

  void
HashGridIndex::serialize(std::ostream& ostream) {
    ostream << DIM << std::endl;

    for (size_t d = 0; d < DIM; d++) {
      ostream << level[d] << " ";
      ostream << index[d] << " ";
    }

    ostream << std::endl;

    ostream << Leaf << std::endl;
  }

  size_t
  HashGridIndex::dim() const {
    return DIM;
  }

  void
HashGridIndex::set(size_t d, HashGridIndex::level_type l, HashGridIndex::index_type i) {
    level[d] = l;
    index[d] = i;
    rehash();
  }

  void
HashGridIndex::set(size_t d, HashGridIndex::level_type l, HashGridIndex::index_type i, bool isLeaf) {
    level[d] = l;
    index[d] = i;
    Leaf = isLeaf;
    rehash();
  }

  void
HashGridIndex::push(size_t d, HashGridIndex::level_type l, HashGridIndex::index_type i) {
    level[d] = l;
    index[d] = i;
  }

  void
HashGridIndex::push(size_t d, HashGridIndex::level_type l, HashGridIndex::index_type i, bool isLeaf) {
    level[d] = l;
    index[d] = i;
    Leaf = isLeaf;
  }

  void
HashGridIndex::get(size_t d, HashGridIndex::level_type& l, HashGridIndex::index_type& i) const {
    l = level[d];
    i = index[d];
  }

  int
  HashGridIndex::getLevel(size_t d) const {
    return level[d];
  }

  int
  HashGridIndex::getIndex(size_t d) const {
    return index[d];
  }

  void
HashGridIndex::setLeaf(bool isLeaf) {
    Leaf = isLeaf;
  }

  bool
  HashGridIndex::isLeaf() {
    return Leaf;
  }

  double
  HashGridIndex::abs(size_t d) const {
    return index[d] * pow(2.0, -static_cast<double>(level[d]));
  }

  double
  HashGridIndex::getCoordBB(size_t d, double q, double t) const {
    return q * (index[d] * pow(2.0, -static_cast<double>(level[d]))) + t;
  }

  double
  HashGridIndex::getCoordStretching(size_t d, Stretching* stretch) {
    return stretch->getCoordinates(level[d], index[d], d);
  }

  bool
  HashGridIndex::isInnerPoint() {
    for (size_t d = 0; d < DIM; d++) {
      if (level[d] == 0)
        return false;
    }

    return true;
  }

  HashGridIndex*
  HashGridIndex::getPointer() {
    return this;
  }

  void
HashGridIndex::rehash() {
    size_t hash = 0xdeadbeef;

    for (size_t d = 0; d < DIM; d++) {
      hash = (1 << level[d]) + index[d] + hash * 65599;
    }

    hash_value = hash;
  }

  size_t
  HashGridIndex::hash() const {
    return hash_value;
  }

  bool
  HashGridIndex::equals(const HashGridIndex& rhs) const {
    for (size_t d = 0; d < DIM; d++) {
      if (level[d] != rhs.level[d]) {
        return false;
      }
    }

    for (size_t d = 0; d < DIM; d++) {
      if (index[d] != rhs.index[d]) {
        return false;
      }
    }

    return true;
  }

  HashGridIndex&
  HashGridIndex::assign(
      const HashGridIndex& rhs) {
    return this->operator=(rhs);
  }

  HashGridIndex&
  HashGridIndex::operator=(
      const HashGridIndex& rhs) {
    if (this == &rhs) {
      return *this;
    }

    if (DIM != rhs.DIM) {

      if (level) {
        delete[] level;
      }

      if (index) {
        delete[] index;
      }

      DIM = rhs.DIM;

      level = new level_type[DIM];
      index = new index_type[DIM];
    }

    for (size_t d = 0; d < DIM; d++) {
      level[d] = rhs.level[d];
      index[d] = rhs.index[d];
    }

    Leaf = rhs.Leaf;

    rehash();
    return *this;
  }

  std::string
  HashGridIndex::toString() {
    std::ostringstream ostream;
    toString(ostream);

    return ostream.str();
  }

  void
HashGridIndex::toString(std::ostream& stream) {
    stream << "[";

    for (size_t i = 0; i < DIM; i++) {
      if (i != 0) {
        stream << ",";
      }

      stream << " " << this->level[i];
      stream << ", " << this->index[i];
    }

    stream << " ]";
  }

  void
HashGridIndex::getCoords(DataVector& p) {
    for (size_t i = 0; i < DIM; i++) {
      if (level[i] == 0) {
        p.set(i, index[i]);
      } else {
        p.set(i, pow(0.5, static_cast<double>(level[i])) * index[i]);
      }
    }
  }

  void
HashGridIndex::getCoordsBB(DataVector& p, BoundingBox& BB) {
    for (size_t i = 0; i < DIM; i++) {
      if (level[i] == 0) {
        p.set(i, (BB.getIntervalWidth(i) * index[i]) + BB.getIntervalOffset(i));
      } else {
        p.set(i,
            (BB.getIntervalWidth(i)
                * (pow(0.5, static_cast<double>(level[i])) * index[i]))
                + BB.getIntervalOffset(i));
      }
    }
  }

  void
HashGridIndex::getCoordsStretching(DataVector& p, Stretching& stretch) {
    for (size_t i = 0; i < DIM; i++) {
      if (level[i] == 0) {
        p.set(i,
            (stretch.getIntervalWidth(i) * index[i])
                + stretch.getIntervalOffset(i));
      } else {
        p.set(i, stretch.getCoordinates(level[i], index[i], i));
      }
    }
  }

  std::string
  HashGridIndex::getCoordsString() {
    std::stringstream return_stream;

    // switch on scientific notation:
    //return_stream << std::scientific;

    for (size_t i = 0; i < DIM; i++) {
      if (level[i] == 0) {
        return_stream << index[i];
      } else {
        return_stream << std::scientific
            << (pow(0.5, static_cast<double>(level[i])) * index[i]);
      }

      if (i < DIM - 1) {
        return_stream << " ";
      }
    }

    return return_stream.str();
  }

  std::string
  HashGridIndex::getCoordsStringBB(BoundingBox& BB) {
    std::stringstream return_stream;

    for (size_t i = 0; i < DIM; i++) {
      if (level[i] == 0) {
        return_stream << std::scientific
            << (BB.getIntervalWidth(i) * index[i]) + BB.getIntervalOffset(i);
      } else {
        return_stream << std::scientific
            << (BB.getIntervalWidth(i)
                * (pow(0.5, static_cast<double>(level[i])) * index[i]))
                + BB.getIntervalOffset(i);
      }

      if (i < DIM - 1) {
        return_stream << " ";
      }
    }

    return return_stream.str();
  }

  std::string
  HashGridIndex::getCoordsStringStretching(Stretching& stretch) {
    std::stringstream return_stream;

    for (size_t i = 0; i < DIM; i++) {
      return_stream << std::scientific
          << stretch.getCoordinates(level[i], index[i], i);

      if (i < DIM - 1) {
        return_stream << " ";
      }
    }

    return return_stream.str();
  }

  HashGridIndex::level_type
  HashGridIndex::getLevelSum() {
    HashGridIndex::level_type levelsum = 0;

    for (size_t i = 0; i < DIM; i++) {
      levelsum += level[i];
    }

    return levelsum;
  }

  HashGridIndex::level_type
  HashGridIndex::getLevelMax() {
    HashGridIndex::level_type levelmax = level[0];

    for (size_t i = 1; i < DIM; i++) {
      levelmax = std::max(levelmax, level[i]);
    }

    return levelmax;
  }

  HashGridIndex::level_type
  HashGridIndex::getLevelMin() {
    HashGridIndex::level_type levelmin = level[0];

    for (size_t i = 1; i < DIM; i++) {
      levelmin = std::min(levelmin, level[i]);
    }

    return levelmin;
  }


}
}

