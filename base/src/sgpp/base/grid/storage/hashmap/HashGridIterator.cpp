/* ****************************************************************************
 * Copyright (C) 2009 Technische Universitaet Muenchen                         *
 * This file is part of the SG++ project. For conditions of distribution and   *
 * use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
 **************************************************************************** */
// @author Jörg Blank (blankj@in.tum.de), Alexander Heinecke (Alexander.Heinecke@mytum.de), Dirk Pflueger (pflueged@in.tum.de)

#include <memory>
#include <string>
#include <sstream>
#include <exception>

#include <sgpp/base/grid/storage/hashmap/HashGridIterator.hpp>

#include <sgpp/base/exception/generation_exception.hpp>
#include <sgpp/base/grid/storage/hashmap/SerializationVersion.hpp>

namespace SGPP {
namespace base {

  HashGridIterator::HashGridIterator(HashGridStorage* storage) :
      storage(storage), index(storage->dim()) {
    for (size_t i = 0; i < storage->dim(); i++) {
      index.push(i, 1, 1);
    }

    index.rehash();
    this->seq_ = storage->seq(&index);
  }


HashGridIterator::HashGridIterator(HashGridIterator& copy) :
      storage(copy.storage), index(copy.storage->dim()) {
    index_type::level_type l;
    index_type::index_type i;

    for (size_t dim = 0; dim < storage->dim(); dim++) {
      copy.get(dim, l, i);
      index.push(dim, l, i);
    }

    index.rehash();
    this->seq_ = storage->seq(&index);
  }

  HashGridIterator::~HashGridIterator() {
  }


  void
HashGridIterator::resetToLevelZero() {
    for (size_t i = 0; i < storage->dim(); i++) {
      index.push(i, 0, 0);
    }

    index.rehash();
    this->seq_ = storage->seq(&index);
  }

  void
HashGridIterator::left_levelzero(size_t dim) {
    index.set(dim, 0, 0);
    this->seq_ = storage->seq(&index);
  }

  void
HashGridIterator::right_levelzero(size_t dim) {
    index.set(dim, 0, 1);
    this->seq_ = storage->seq(&index);
  }

  void
HashGridIterator::left_child(size_t dim) {
    index_type::level_type l;
    index_type::index_type i;
    index.get(dim, l, i);
    index.set(dim, l + 1, 2 * i - 1);
    this->seq_ = storage->seq(&index);
  }

  void
HashGridIterator::right_child(size_t dim) {
    index_type::level_type l;
    index_type::index_type i;
    index.get(dim, l, i);
    index.set(dim, l + 1, 2 * i + 1);
    this->seq_ = storage->seq(&index);
  }

  void
HashGridIterator::top(size_t d) {
    index.set(d, 1, 1);
    this->seq_ = storage->seq(&index);
  }

  void
HashGridIterator::up(size_t d) {
    index_type::level_type l;
    index_type::index_type i;
    index.get(d, l, i);

    i /= 2;
    i += i % 2 == 0 ? 1 : 0;

    index.set(d, l - 1, i);
    this->seq_ = storage->seq(&index);
  }

  void
HashGridIterator::step_left(size_t d) {
    index_type::level_type l;
    index_type::index_type i;
    index.get(d, l, i);
    index.set(d, l, i - 2);
    this->seq_ = storage->seq(&index);

  }

  void
HashGridIterator::step_right(size_t d) {
    index_type::level_type l;
    index_type::index_type i;
    index.get(d, l, i);
    index.set(d, l, i + 2);
    this->seq_ = storage->seq(&index);

  }

  bool
  HashGridIterator::isInnerPoint() {
    return index.isInnerPoint();
  }

  bool
  HashGridIterator::hint() const {
    return storage->get(this->seq_)->isLeaf();
  }

  bool
  HashGridIterator::hint_left(size_t d) {
    index_type::level_type l;
    index_type::index_type i;
    bool hasIndex = true;

    index.get(d, l, i);
    index.set(d, l + 1, 2 * i - 1);

    HashGridIndex* my_Index = index.getPointer();
    hasIndex = storage->has_key(my_Index);

    index.set(d, l, i);

    return hasIndex;
  }

  bool
  HashGridIterator::hint_right(size_t d) {
    index_type::level_type l;
    index_type::index_type i;
    bool hasIndex = true;

    index.get(d, l, i);
    index.set(d, l + 1, 2 * i + 1);

    HashGridIndex* my_Index = index.getPointer();
    hasIndex = storage->has_key(my_Index);

    index.set(d, l, i);

    return hasIndex;
  }

  void
HashGridIterator::get(size_t d, index_type::level_type& l,
      index_type::index_type& i) const {
    index.get(d, l, i);
  }

  void
HashGridIterator::set(size_t d, index_type::level_type l,
      index_type::index_type i) {
    index.set(d, l, i);
    this->seq_ = storage->seq(&index);
  }

  void
HashGridIterator::push(size_t d, index_type::level_type l,
      index_type::index_type i) {
    index.push(d, l, i);
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
      if (this->hint_left(dim)) {
        depth++;
        this->left_child(dim);
      } else if (this->hint_right(dim)) {
        depth++;
        this->right_child(dim);
      } else {

        index.get(dim, cur_level, cur_index);

        bool hasFound = false; //Was a next index found?

        //Ok, we have no more childs left. Now we slide from left to right in the dim on
        //the same level, to see, if there are adaptive refinements
        for (size_t i = cur_index + 2; i < (unsigned int) (1 << (depth));
            i = i + 2) {
          this->set(dim, cur_level, *reinterpret_cast<unsigned int*>(&i));

          //does this index exist?
          if (!storage->end(this->seq())) {
            if (this->hint_left(dim)) {
              depth++;
              this->left_child(dim);
              hasFound = true;
              break;
            } else if (this->hint_right(dim)) {
              depth++;
              this->right_child(dim);
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

} // namespace base
} // namespace SGPP

