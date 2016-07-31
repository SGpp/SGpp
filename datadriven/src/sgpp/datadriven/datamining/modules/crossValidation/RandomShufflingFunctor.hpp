/*
 * RandomShufflingFunctor.hpp
 *
 *  Created on: 31.07.2016
 *      Author: Michael Lettrich
 */
#pragma once

#include <sgpp/datadriven/datamining/modules/crossValidation/ShufflingFunctor.hpp>

#include <algorithm>

namespace sgpp {
namespace datadriven {

class RandomShufflingFunctor : public ShufflingFunctor {
 public:
  RandomShufflingFunctor() : ShufflingFunctor(){};
  virtual ~RandomShufflingFunctor(){};
  virtual void shuffle(std::vector<size_t>& indices) const {
    std::shuffle(indices.begin(), indices.end(), generator);
  };
};

} /* namespace datadriven */
} /* namespace sgpp */
