// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#pragma once

#include <sgpp/globaldef.hpp>

namespace sgpp {
namespace datadriven {

/**
 * Struct that stores all the configuration information
 * for parallelization with ScaLAPACK.
 */
struct ParallelConfiguration {
  // disable by default, enable if config is found. Does not have to be set in
  // the config file.
  bool scalapackEnabled_ = false;

  int processRows_ = -1;
  int processCols_ = -1;
  size_t rowBlockSize_ = 64;
  size_t columnBlockSize_ = 64;

  /*
  // Debug method to neatly print internal data
  void dumpToStream(std::ostream& stream_out = std::cout) const {
    stream_out << "scalapackEnabled: \t\t" << std::boolalpha
               << scalapackEnabled_ << std::endl;
    stream_out << "processRows: \t\t" << processRows_ << std::endl;
    stream_out << "processCols: \t\t" << processCols_ << std::endl;
    stream_out << "rowBlockSize: \t\t" << rowBlockSize_ << std::endl;
    stream_out << "columnBlockSize: \t\t" << columnBlockSize_ << std::endl;
  }
  */
};

}  // namespace datadriven
}  // namespace sgpp
