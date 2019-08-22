// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <sgpp/globaldef.hpp>
#include <sgpp/base/grid/common/IndexInSubspaceGenerator.hpp>

namespace sgpp {
namespace base {

IndexInSubspaceGenerator::IndexInSubspaceGenerator(const value_type&
    level_vector) :
  level_vector(level_vector),
  max_index_vector(level_vector.size()), val_(nullptr), dim_(level_vector.size()) {
  for (size_t d = 0; d < dim_; d++) {
    if (level_vector[d] > 0) {
      max_index_vector[d] = (2 << (level_vector[d] - 1)) - 1;
    } else {
      max_index_vector[d] = 1;
    }
  }
}

IndexInSubspaceGenerator::iterator IndexInSubspaceGenerator::begin() {
  val_ = std::make_shared<value_type>(dim_);

  for (size_t d = 0; d < dim_; d++) {
    if (level_vector[d] > 0) {
      (*val_)[d] = 1;
    } else {
      (*val_)[d] = 0;
    }
  }

  queue_value_type pair(val_, 0);
  queue_.push(pair);
  this->next_();
  return iterator(this);
}


bool IndexInSubspaceGenerator::compareVectors(value_type& vec1,
    value_type& vec2) {
  if (vec1.size() != vec2.size()) {
    throw std::length_error("Vector size mismatch");
  }

  bool res = true;

  for (size_t i = 0; i < vec1.size(); i++) {
    if (vec1[i] > vec2[i]) {
      res = false;
      break;
    }
  }

  return res;
}

IndexInSubspaceGenerator* IndexInSubspaceGenerator::next_() {
  queue_value_type pair = queue_.front();
  val_ = pair.first;
  const size_t startDim = pair.second;

  for (size_t d = startDim; d < this->dim_; d++) {
    pointer_type new_index = pointer_type(new value_type(*val_));

    if (level_vector[d] > 0) {
      (*new_index)[d] += 2;
    } else {
      (*new_index)[d]++;
    }

    if (compareVectors(*new_index, max_index_vector)) {
      queue_.push(std::make_pair(new_index, d));
    }
  }

  queue_.pop();
  return this;
}

}  // namespace base
}  // namespace sgpp
