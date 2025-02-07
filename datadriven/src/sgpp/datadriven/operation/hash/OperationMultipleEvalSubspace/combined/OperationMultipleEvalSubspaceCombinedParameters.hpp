// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#pragma once

/*
 * Don't remove the "ifndef", they are required to overwrite the parameters though compiler's "-D"
 *
 */

// number of data elements processed in parallel
// should be divisible by vector size
// corresponds to data chunk size
#ifndef X86COMBINED_PARALLEL_DATA_POINTS
// #define X86COMBINED_PARALLEL_DATA_POINTS 4
#define X86COMBINED_PARALLEL_DATA_POINTS 256
// #define X86COMBINED_PARALLEL_DATA_POINTS 192
#endif

#ifndef X86COMBINED_ENABLE_SUBSPACE_SKIPPING
#define X86COMBINED_ENABLE_SUBSPACE_SKIPPING 1
#endif

#ifndef X86COMBINED_UNROLL
#define X86COMBINED_UNROLL 0
// implies X86COMBINED_VEC_PADDING == 8
#endif

#ifndef X86COMBINED_VEC_PADDING
#define X86COMBINED_VEC_PADDING 4
// #define X86COMBINED_VEC_PADDING 8
// #define X86COMBINED_VEC_PADDING 24
#endif

#ifndef X86COMBINED_STREAMING_THRESHOLD
// good value: #define X86COMBINED_STREAMING_THRESHOLD 128
#define X86COMBINED_STREAMING_THRESHOLD 128
#endif

#ifndef X86COMBINED_LIST_RATIO
// good value: #define X86COMBINED_LIST_RATIO 0.1
#define X86COMBINED_LIST_RATIO 0.2
#endif

#ifndef X86COMBINED_ENABLE_PARTIAL_RESULT_REUSAGE
#define X86COMBINED_ENABLE_PARTIAL_RESULT_REUSAGE 1
#endif

// only set from the outside
// #define X86COMBINED_WRITE_STATS "stats.out"
