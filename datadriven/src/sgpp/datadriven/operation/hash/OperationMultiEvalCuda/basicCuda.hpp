#ifndef BASICCUDA_HPP
#define BASICCUDA_HPP

#include <stdint.h>

/// Magic numbers for computing polynomial basis functions
const double idxtable[] = {1.0,2.0,-2.0,-1.0};
/// Magic numbers for computing polynomial basis functions
const int32_t idxtablei[]  = {1,2,-2,-1};

/// Struct containing all important information about a grid node in one dimension
typedef struct {
  uint32_t index; /// The node index
  uint32_t level; /// The node level
  uint32_t child[2];  /// Indices of child nodes
  uint32_t grad;  /// Maximal polynomial degree
  double x; /// Absolute position in hypercube
  double level2;  /// 2^level
} gridnode_t;

/// Struct for index limitation of evaluation points
typedef struct {
  uint32_t lower; /// Lower boundary
  uint32_t upper; /// Upper boundary
} limit_t;

/// Block size for kernel call on the GPU
#define CUDA_BLOCKSIZE 128

/// Index for children of leaf nodes
#define GRID_END 0xFFFFFFFF
/// Index of root node
#define GRID_START 0x00000000

#endif // BASICCUDA_HPP
