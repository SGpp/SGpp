// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <sgpp/base/grid/common/SubspaceGenerator.hpp>
#include <sgpp/globaldef.hpp>

namespace sgpp {
namespace base {

SubspaceGenerator::SubspaceGenerator(unsigned int dim,
                                     unsigned int max_level) :
  val_(nullptr), dim_(dim), max_sum_(max_level + dim - 1) {
  pointer_type root_subspace = new value_type(dim, 1);
  queue_value_type pair(root_subspace, 0);
  queue_.push(pair);
  this->next_();
}



SubspaceGenerator* SubspaceGenerator::next_() {
  queue_value_type pair = queue_.front();
  val_ = pair.first;
  unsigned int start_dim = pair.second;
  unsigned int sum = std::accumulate(val_->begin(), val_->end(), 0);

  if (sum < max_sum_) {
    for (unsigned int d = start_dim; d < this->dim_; d++) {
      pointer_type new_subspace = new value_type(*val_);
      (*new_subspace)[d] += 1;
      queue_.push(std::make_pair(new_subspace, d));
    }
  }

  queue_.pop();
  return this;
}

}  // namespace base
}  // namespace sgpp
