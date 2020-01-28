// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <sgpp/globaldef.hpp>

#include <iostream>  // std::cout
#include <algorithm>  // std::nth_element
#include <vector>  // std::vector
#include <utility>  // std::pair
#include <iterator>
#include <queue>
#include <numeric>

#ifndef SRC_SGPP_BASE_GRID_COMMON_SUBSPACEGENERATOR_HPP_
#define SRC_SGPP_BASE_GRID_COMMON_SUBSPACEGENERATOR_HPP_

namespace sgpp {
namespace base {

class SubspaceGenerator {
 public:
  typedef std::vector<unsigned int> value_type;
  typedef value_type* pointer_type;

  typedef std::pair<pointer_type, unsigned int> queue_value_type;

  explicit SubspaceGenerator(unsigned int dim, unsigned int max_level);

  class iterator
    : std::iterator<std::forward_iterator_tag, pointer_type> {
   public:
    explicit iterator(SubspaceGenerator* p = nullptr) : ptr_(p) {}
    // implicit copy constructor, copy assignment and destructor

    reference operator* () {
      return ptr_->val_;
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
    SubspaceGenerator* ptr_;
  };




  iterator begin() {
    return iterator(this);
  }
  iterator end()   {
    return iterator(nullptr);
  }

 private:
  pointer_type val_;
  unsigned int dim_;
  unsigned int max_sum_;
  std::queue <queue_value_type> queue_;

  SubspaceGenerator* next_();
};

}  // namespace base
}  // namespace sgpp

#endif /* SRC_SGPP_BASE_GRID_COMMON_SUBSPACEGENERATOR_HPP_ */
