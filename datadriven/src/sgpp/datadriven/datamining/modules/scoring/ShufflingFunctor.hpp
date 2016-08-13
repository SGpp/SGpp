/* Copyright (C) 2008-today The SG++ project
 * This file is part of the SG++ project. For conditions of distribution and
 * use, please see the copyright notice provided with SG++ or at
 * sgpp.sparsegrids.org
 *
 * ShufflingFunctor.hpp
 *
 *  Created on: 31.07.2016
 *      Author: Michael Lettrich
 */

#pragma once

#include <memory>
#include <random>
#include <vector>

namespace sgpp {
namespace datadriven {

class ShufflingFunctor {
 public:
  ShufflingFunctor() {
    std::random_device rd;
    seed = rd();
    seed = 1;
    generator = std::mt19937(rd());
  };
  virtual ~ShufflingFunctor(){};

  virtual void shuffle(std::vector<size_t>& indices) = 0;
  int64_t getSeed() const { return seed; }
  void setSeed(int64_t seed) {
    this->seed = seed;
    generator = std::mt19937(seed);
  }

 protected:
  int64_t seed;
  std::mt19937 generator;
};

} /* namespace datadriven */
} /* namespace sgpp */
