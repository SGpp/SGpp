

#include "IndexInSubspaceGenerator.hpp"

namespace sg {
  namespace base {

    IndexInSubspaceGenerator::IndexInSubspaceGenerator (const value_type& level_vector) :
      level_vector(level_vector),
      max_index_vector(level_vector.size()), val_(NULL), dim_(level_vector.size()) {
      for (size_t d = 0; d < dim_; d++) {
        max_index_vector[d] = (2 << (level_vector[d] - 1)) - 1;
      }

    }

    IndexInSubspaceGenerator::iterator IndexInSubspaceGenerator::begin() {
      val_ = std::make_shared<value_type>(dim_, 1);
      queue_value_type pair(val_, 0);
      queue_.push(pair);
      this->next_();
      return iterator(this);
    }


    bool IndexInSubspaceGenerator::compareVectors(value_type& vec1, value_type& vec2) {
      if (vec1.size() != vec2.size()) {
        throw new std::length_error("Vector size mismatch");
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
      unsigned int start_dim = pair.second;

      for (unsigned int d = start_dim; d < this->dim_; d++) {
        pointer_type new_index = pointer_type(new value_type(*val_));
        (*new_index)[d] += 2;

        if (compareVectors(*new_index, max_index_vector)) {
          queue_.push(std::make_pair(new_index, d));
        }
      }

      queue_.pop();
      return this;
    }




  } /* namespace base */
} /* namespace sg */
