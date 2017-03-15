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

#include <vector>

namespace sgpp {
namespace datadriven {

/**
 * Sequential Shuffling does not permute indices at all and thus keeps their order unchanged. To be
 * used on  datasets that are already arranged in the desired manner.
 */
class SequentialShufflingFunctor : public ShufflingFunctor {
 public:
  ShufflingFunctor* clone() const override;

  /**
   * Does not permute indices and just returns the same permutation passed as a parameter. To be
   * used on datasets that are already arranged as desired.
   * @param indices vector of indices to permute.
   */
  void shuffle(std::vector<size_t>& indices) override;
};

} /* namespace datadriven */
} /* namespace sgpp */
