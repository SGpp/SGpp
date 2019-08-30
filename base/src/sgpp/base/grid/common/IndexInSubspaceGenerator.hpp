// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <sgpp/globaldef.hpp>
#include <sgpp/base/grid/LevelIndexTypes.hpp>

#include <algorithm>  // std::nth_element
#include <cmath>
#include <iostream>  // std::cout
#include <iterator>
#include <memory>
#include <numeric>
#include <queue>
#include <stdexcept>
#include <utility>  // std::pair
#include <vector>  // std::vector

#ifndef SRC_SGPP_BASE_GRID_COMMON_INDEXINSUBSPACEGENERATOR_HPP_
#define SRC_SGPP_BASE_GRID_COMMON_INDEXINSUBSPACEGENERATOR_HPP_

namespace sgpp {
namespace base {


/**
 * Container for the index_vectors of a subspace.
 * The subclass iterator implements the STL forward iterator
 *
 * The indices are generated on the fly and hence only one iterator instance can
 * be possible (second call to begin() will empty the working queue)
 */
class IndexInSubspaceGenerator {
 public:
  typedef std::vector<index_t> value_type;
  typedef std::shared_ptr<value_type> pointer_type;
  typedef std::pair<pointer_type, size_t> queue_value_type;

  /**
   * Constructor
   * @param level_vector std::vector<int> with the level vector of the subspace
   */
  explicit IndexInSubspaceGenerator(const value_type& level_vector);


  /**
   * Destructor
   */
  ~IndexInSubspaceGenerator() {
    this->queue_.empty();
  }


  /**
   * Iterator class compatible with STL forward iterator (no const iterator)
   */
  class iterator
    : std::iterator<std::forward_iterator_tag, value_type> {
   public:
    explicit iterator(IndexInSubspaceGenerator* p = nullptr) : ptr_(p) {}
    // implicit copy constructor, copy assignment and destructor

    reference operator* () {
      return *(ptr_->val_);
    }

    iterator& operator++ () {
      if (!ptr_->queue_.empty()) {
        ptr_ = ptr_->next_();
      } else {
        ptr_ = nullptr;
      }

      return *this;
    }

    iterator operator++(int) {
      iterator tmp = *this;
      ++*this;
      return tmp;
    }

    bool operator== (const iterator& other) const {
      return ptr_ == other.ptr_;
    }
    bool operator!= (const iterator& other) const {
      return ptr_ != other.ptr_;
    }

   private:
    IndexInSubspaceGenerator* ptr_;
  };



  /**
   * returns an iterator to the beginning
   *
   * @return iterator to the beginning
   */
  iterator begin();

  /**
   * returns an iterator to the end;
   *
   * @return iterator to the end
   */
  iterator end()   {
    return iterator(nullptr);
  }

 private:
  value_type level_vector;
  value_type max_index_vector;
  pointer_type val_;
  size_t dim_;
  std::queue <queue_value_type> queue_;

  /**
   * Assures that the vectors are of the same size and vec1 is not larger than vec2
   * @param vec1 first vector to compare
   * @param vec2 second vector to compare
   * @return true if vec1 components-wise smaller or equal to vec2
   */
  bool compareVectors(value_type& vec1, value_type& vec2);


  /**
   * Move to the next element in the queue
   * @return pointer to self
   */
  IndexInSubspaceGenerator* next_();
};

}  // namespace base
}  // namespace sgpp

#endif /* SRC_SGPP_BASE_GRID_COMMON_INDEXINSUBSPACEGENERATOR_HPP_ */
