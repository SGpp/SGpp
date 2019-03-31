/*
 * Copyright (C) 2008-today The SG++ project
 * This file is part of the SG++ project. For conditions of distribution and
 * use, please see the copyright notice provided with SG++ or at
 * sgpp.sparsegrids.org
 *
 * ParallelConfiguration.hpp
 *
 * Created on: Mar 13, 2019
 *     Author: Jan Schopohl
 */

#pragma once

#include <sgpp/globaldef.hpp>

namespace sgpp {
namespace datadriven {

/**
 * Struct that stores all the configuration information
 * for parallelization with ScaLAPACK.
 */
struct ParallelConfiguration {
  int processRows_ = -1;
  int processCols_ = -1;
  size_t rowBlockSize_ = 128;
  size_t columnBlockSize_ = 128;
};

}  // namespace datadriven
}  // namespace sgpp
