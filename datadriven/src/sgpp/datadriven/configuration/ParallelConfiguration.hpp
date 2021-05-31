// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#pragma once

#include <sgpp/globaldef.hpp>

namespace sgpp {
namespace datadriven {

/**
 * Struct that stores all the configuration information for parallelization with ScaLAPACK.
 */
struct ParallelConfiguration {
  // disable by default, enable if config is found. Does not have to be set in the config file.
  bool scalapackEnabled_ = false;

  int processRows_ = -1;
  int processCols_ = -1;
  size_t rowBlockSize_ = 64;
  size_t columnBlockSize_ = 64;
};

}  // namespace datadriven
}  // namespace sgpp
