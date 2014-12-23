#pragma once

//number of data elements processed in parallel
//should be divisible by vector size
//corresponds to data chunk size
#ifndef X86VECTORIZED_PARALLEL_DATA_POINTS
#define X86VECTORIZED_PARALLEL_DATA_POINTS 128
//#define X86VECTORIZED_PARALLEL_DATA_POINTS 192
#endif

#ifndef X86VECTORIZED_ENABLE_SUBSPACE_SKIPPING
#define X86VECTORIZED_ENABLE_SUBSPACE_SKIPPING 1
#endif

#define X86VECTORIZED_VEC_PADDING 4
//#define X86VECTORIZED_VEC_PADDING 8
//#define X86VECTORIZED_VEC_PADDING 24
