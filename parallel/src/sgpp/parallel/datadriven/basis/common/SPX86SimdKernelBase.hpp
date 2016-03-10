// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#ifndef SPX86SIMDKERNELBASE_HPP
#define SPX86SIMDKERNELBASE_HPP

#include <sgpp/parallel/datadriven/basis/common/CommonX86SimdKernelBase.hpp>

#include <sgpp/globaldef.hpp>

namespace sgpp {
namespace parallel {

class SPX86SimdKernelBase {
 public:
  static inline size_t getChunkGridPoints() { return 12; }
  static inline size_t getChunkDataPoints() {
    return 48;  // must be divisible by 48
  }
  static inline void resetKernel() {}
};

}  // namespace parallel
}  // namespace sgpp

#endif  // SPX86SIMDKERNELBASE_HPP
