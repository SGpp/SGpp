#pragma once

//number of data elements processed in parallel
//should be divisible by vector size
//corresponds to data chunk size
#ifndef X86BASELINE_PARALLEL_DATA_POINTS
#define X86BASELINE_PARALLEL_DATA_POINTS 128
#endif

#ifndef X86BASELINE_ENABLE_SUBSPACE_SKIPPING
#define X86BASELINE_ENABLE_SUBSPACE_SKIPPING 1
#endif
