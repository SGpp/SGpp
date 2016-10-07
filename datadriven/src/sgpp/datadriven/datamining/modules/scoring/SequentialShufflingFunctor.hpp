/*
 * Copyright (C) 2008-today The SG++ project
 * This file is part of the SG++ project. For conditions of distribution and
 * use, please see the copyright notice provided with SG++ or at
 * sgpp.sparsegrids.org
 *
 * SequentialShufflingFunctor.hpp
 *
 * Created on: Oct 7, 2016
 *     Author: Michael Lettrich
 */

#pragma once

#include <sgpp/datadriven/datamining/modules/scoring/ShufflingFunctor.hpp>

namespace sgpp {
namespace datadriven {

class SequentialShufflingFunctor : public ShufflingFunctor {
 public:
  SequentialShufflingFunctor() : ShufflingFunctor(){};
  virtual ~SequentialShufflingFunctor(){};
  virtual void shuffle(std::vector<size_t>& indices){
      // doesn't do any permutation to provided indices, since we're using the entries sequentially
  };
};

} /* namespace datadriven */
} /* namespace sgpp */
