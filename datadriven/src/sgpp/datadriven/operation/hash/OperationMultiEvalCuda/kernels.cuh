#include <stdint.h>
#include "basicCuda.hpp"

/// Initializes all constants used in devvice code
void initCudaConstants (const uint32_t DIM);

/// Kernel to compute index in Morton order curve
__global__ void gpu_zindex (double* pos, uint64_t* index);

/// Kernel for standard evaluation
__global__ void gpu_eval (double *res, gridnode_t* node, double* a, double* pos,
    uint32_t* limit, uint32_t subcnt, uint32_t* subs);

/// Kernel for transposed evaluation with additional FMA
__global__ void gpu_transevel (double* a, gridnode_t* node, double* pos, double* y,
    limit_t* limit, double* b, double c, double M_);

/// Kernel for preprocessing limitations of used subspaces for each evaluation point
__global__ void gpu_preprocess (uint32_t* limit, gridnode_t* node, double* pos,
    uint32_t maxlevel);

/// Kernel for transposed evaluation
__global__ void gpu_transevel (double* a, gridnode_t* node, double* pos, double* y,
    limit_t* limit);

/// Kernel for preprocessing of straming boundary limitations
__global__ void gpu_zbound (uint64_t* idx_p, gridnode_t* node, limit_t* limit, uint32_t M);
