// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#ifndef X86SIMDKERNELBASE_HPP
#define X86SIMDKERNELBASE_HPP

#include <sgpp/parallel/datadriven/basis/common/CommonX86SimdKernelBase.hpp>

#include <sgpp/globaldef.hpp>


namespace SGPP {
  namespace parallel {

    class X86SimdKernelBase {
      public:
        static inline size_t getChunkGridPoints() {
          return 12;
        }
        static inline size_t getChunkDataPoints() {
          return 24; // must be divisible by 24
        }
        static inline void resetKernel() {}
    };

    class X86SimdKernelBase1 {
      public:
        static inline size_t getChunkGridPoints() {
          return 1;
        }
        static inline size_t getChunkDataPoints() {
          return 1;
        }
        static inline void resetKernel() {}
    };

  }
}

#endif // X86SIMDKERNELBASE_HPP