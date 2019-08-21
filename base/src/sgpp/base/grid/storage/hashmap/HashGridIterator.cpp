// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <sgpp/base/grid/storage/hashmap/HashGridIterator.hpp>
#include <sgpp/base/exception/generation_exception.hpp>
#include <sgpp/base/grid/storage/hashmap/SerializationVersion.hpp>

#include <memory>
#include <string>
#include <sstream>
#include <exception>

namespace sgpp {
namespace base {

HashGridIterator::HashGridIterator(HashGridStorage& storage) :
  storage(storage), index(storage.getDimension()) {
  for (size_t i = 0; i < storage.getDimension(); i++) {
    index.push(i, 1, 1);
  }

  index.rehash();
  this->seq_ = storage.getSequenceNumber(index);
}


HashGridIterator::HashGridIterator(HashGridIterator& copy) :
  storage(copy.storage), index(copy.storage.getDimension()) {
  index_type::level_type l;
  index_type::index_type i;

  for (size_t dim = 0; dim < storage.getDimension(); dim++) {
    copy.get(dim, l, i);
    index.push(dim, l, i);
  }

  index.rehash();
  this->seq_ = storage.getSequenceNumber(index);
}

HashGridIterator::~HashGridIterator() {
}


void
HashGridIterator::resetToLevelZero() {
  for (size_t i = 0; i < storage.getDimension(); i++) {
    index.push(i, 0, 0);
  }

  index.rehash();
  this->seq_ = storage.getSequenceNumber(index);
}

void
HashGridIterator::resetToLeftLevelZero(size_t dim) {
  index.set(dim, 0, 0);
  this->seq_ = storage.getSequenceNumber(index);
}

void
HashGridIterator::resetToRightLevelZero(size_t dim) {
  index.set(dim, 0, 1);
  this->seq_ = storage.getSequenceNumber(index);
}

void
HashGridIterator::resetToLevelOne(size_t d) {
  index.set(d, 1, 1);
  this->seq_ = storage.getSequenceNumber(index);
}

void
HashGridIterator::leftChild(size_t dim) {
  index_type::level_type l;
  index_type::index_type i;
  index.get(dim, l, i);
  index.set(dim, l + 1, 2 * i - 1);
  this->seq_ = storage.getSequenceNumber(index);
}

void
HashGridIterator::rightChild(size_t dim) {
  index_type::level_type l;
  index_type::index_type i;
  index.get(dim, l, i);
  index.set(dim, l + 1, 2 * i + 1);
  this->seq_ = storage.getSequenceNumber(index);
}

void
HashGridIterator::up(size_t d) {
  index_type::level_type l;
  index_type::index_type i;
  index.get(d, l, i);

  i /= 2;
  i += i % 2 == 0 ? 1 : 0;

  index.set(d, l - 1, i);
  this->seq_ = storage.getSequenceNumber(index);
}

void
HashGridIterator::stepLeft(size_t d) {
  index_type::level_type l;
  index_type::index_type i;
  index.get(d, l, i);
  index.set(d, l, i - 2);
  this->seq_ = storage.getSequenceNumber(index);
}

void
HashGridIterator::stepRight(size_t d) {
  index_type::level_type l;
  index_type::index_type i;
  index.get(d, l, i);
  index.set(d, l, i + 2);
  this->seq_ = storage.getSequenceNumber(index);
}

bool
HashGridIterator::isInnerPoint() const {
  return index.isInnerPoint();
}

bool
HashGridIterator::hint() const {
  return storage.getPoint(this->seq_).isLeaf();
}

bool
HashGridIterator::hintLeft(size_t d) {
  index_type::level_type l;
  index_type::index_type i;
  bool hasIndex = true;

  index.get(d, l, i);
  index.set(d, l + 1, 2 * i - 1);

  hasIndex = storage.isContaining(index);

  index.set(d, l, i);

  return hasIndex;
}

bool
HashGridIterator::hintRight(size_t d) {
  index_type::level_type l;
  index_type::index_type i;
  bool hasIndex = true;

  index.get(d, l, i);
  index.set(d, l + 1, 2 * i + 1);

  hasIndex = storage.isContaining(index);

  index.set(d, l, i);

  return hasIndex;
}

size_t
HashGridIterator::seq() const {
  return seq_;
}

HashGridIterator::level_t
HashGridIterator::getGridDepth(size_t dim) {
  index_type::level_type depth = 1;
  index_type::level_type orig_level, cur_level;
  index_type::index_type orig_index, cur_index;

  index.get(dim, orig_level, orig_index);

  while (true) {
    if (this->hintLeft(dim)) {
      depth++;
      this->leftChild(dim);
    } else if (this->hintRight(dim)) {
      depth++;
      this->rightChild(dim);
    } else {
      index.get(dim, cur_level, cur_index);

      bool hasFound = false;  // Was a next index found?

      // Ok, we have no more childs left.
      // Now we slide from left to right in the dim on
      // the same level, to see, if there are adaptive refinements
      for (size_t i = cur_index + 2; i < static_cast<unsigned int>(1 << (depth));
           i = i + 2) {
        this->set(dim, cur_level, static_cast<unsigned int>(i));

        // does this index exist?
        if (!storage.isInvalidSequenceNumber(this->seq())) {
          if (this->hintLeft(dim)) {
            depth++;
            this->leftChild(dim);
            hasFound = true;
            break;
          } else if (this->hintRight(dim)) {
            depth++;
            this->rightChild(dim);
            hasFound = true;
            break;
          }
        }
      }

      if (!hasFound) {
        break;
      }
    }
  }

  this->set(dim, orig_level, orig_index);
  return depth;
}

std::string
HashGridIterator::toString() {
  return index.toString();
}

}  // namespace base
}  // namespace sgpp

