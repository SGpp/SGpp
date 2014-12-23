#pragma once

/*
 * Don't remove the "ifndef", they are required to overwrite the parameters though compiler's "-D"
 *
 */

//number of data elements processed in parallel
//should be divisible by vector size
//corresponds to data chunk size
#ifndef X86PERFECT_PARALLEL_DATA_POINTS
//#define X86PERFECT_PARALLEL_DATA_POINTS 4
#define X86PERFECT_PARALLEL_DATA_POINTS 256
//#define X86PERFECT_PARALLEL_DATA_POINTS 192
#endif

#ifndef X86PERFECT_ENABLE_SUBSPACE_SKIPPING
#define X86PERFECT_ENABLE_SUBSPACE_SKIPPING 1
#endif

#define X86PERFECT_VEC_PADDING 4
//#define X86PERFECT_VEC_PADDING 8
//#define X86PERFECT_VEC_PADDING 24

#ifndef X86PERFECT_STREAMING_THRESHOLD
// good value: #define X86PERFECT_STREAMING_THRESHOLD 128
#define X86PERFECT_STREAMING_THRESHOLD 128
#endif

#ifndef X86PERFECT_LIST_RATIO
// good value: #define X86PERFECT_LIST_RATIO 0.1
#define X86PERFECT_LIST_RATIO 0.2
#endif

//only set from the outside
//#define X86PERFECT_WRITE_STATS "stats.out"

