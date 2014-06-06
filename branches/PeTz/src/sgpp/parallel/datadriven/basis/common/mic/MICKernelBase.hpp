/* ****************************************************************************
* Copyright (C) 2013 Technische Universitaet Muenchen                         *
* This file is part of the SG++ project. For conditions of distribution and   *
* use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
**************************************************************************** */
// @author Alexander Heinecke (Alexander.Heinecke@mytum.de)
// @author Roman Karlstetter (karlstetter@mytum.de)

#ifndef MICKERNELBASE_HPP
#define MICKERNELBASE_HPP

#include "../KernelBase.hpp"

#ifdef USEMIC

#ifndef MIC_UNROLLING_WIDTH
#define MIC_UNROLLING_WIDTH 96
#endif


#if (MIC_UNROLLING_WIDTH != 24) && \
  (MIC_UNROLLING_WIDTH != 32) && \
  (MIC_UNROLLING_WIDTH != 48) && \
  (MIC_UNROLLING_WIDTH != 64) && \
  (MIC_UNROLLING_WIDTH != 72) && \
  (MIC_UNROLLING_WIDTH != 80) && \
  (MIC_UNROLLING_WIDTH != 96)
#error MIC_UNROLLING_WIDTH has to be one of 24, 32, 48, 64, 72, 80, 96.
#endif

#ifndef MIC_UNROLLING_WIDTH_SP
#define MIC_UNROLLING_WIDTH_SP 192
#endif

#if (MIC_UNROLLING_WIDTH_SP != 48) && \
  (MIC_UNROLLING_WIDTH_SP != 64) && \
  (MIC_UNROLLING_WIDTH_SP != 96) && \
  (MIC_UNROLLING_WIDTH_SP != 128) && \
  (MIC_UNROLLING_WIDTH_SP != 144) && \
  (MIC_UNROLLING_WIDTH_SP != 160) && \
  (MIC_UNROLLING_WIDTH_SP != 192)
#error MIC_UNROLLING_WIDTH_SP has to be one of 48, 64, 96, 128, 144, 160, 192.
#endif

#ifdef __INTEL_OFFLOAD
#pragma offload_attribute(push, target(mic))
#endif

namespace sg {
  namespace parallel {

    class MICKernelBase {
      public:
        static inline size_t getChunkGridPoints() {
          return 12;
        }
        static inline size_t getChunkDataPoints() {
          return MIC_UNROLLING_WIDTH;
        }
    };

    class SPMICKernelBase {
      public:
        static inline size_t getChunkGridPoints() {
          return 24;
        }
        static inline size_t getChunkDataPoints() {
          return MIC_UNROLLING_WIDTH_SP; // must be divisible by 96
        }
    };
  }
}

#ifdef __INTEL_OFFLOAD
#pragma offload_attribute(pop)
#endif

#endif // USEMIC
#endif // MICKERNELBASE_HPP
