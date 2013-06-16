/* ****************************************************************************
* Copyright (C) 2010-2013 Technische Universitaet Muenchen                    *
* This file is part of the SG++ project. For conditions of distribution and   *
* use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
**************************************************************************** */
// @author Alexander Heinecke (Alexander.Heinecke@mytum.de)
// @author Roman Karlstetter (karlstetter@mytum.de)

#ifndef SPX86SIMDKERNELBASE_HPP
#define SPX86SIMDKERNELBASE_HPP

#include "CommonX86SimdKernelBase.hpp"

namespace sg {
  namespace parallel {

    class SPX86SimdKernelBase {
      public:
        static inline size_t getChunkGridPoints() {
          return 12;
        }
        static inline size_t getChunkDataPoints() {
          return 48; //must be divisible by 48
        }
        static inline void resetKernel() {}
    };

  }
}

#endif // SPX86SIMDKERNELBASE_HPP
