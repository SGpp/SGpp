/* ****************************************************************************
* Copyright (C) 2010 Technische Universitaet Muenchen                         *
* This file is part of the SG++ project. For conditions of distribution and   *
* use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
**************************************************************************** */
// @author Alexander Heinecke (Alexander.Heinecke@mytum.de)
// @author Roman Karlstetter (karlstetter@mytum.de)

#ifndef X86SIMDKERNELBASE_HPP
#define X86SIMDKERNELBASE_HPP

#include "CommonX86SimdKernelBase.hpp"

namespace sg {
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
