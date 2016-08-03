#ifndef MULTIEVALKERNEL_HPP
#define MULTIEVALKERNEL_HPP

#include <stdint.h>

/// Kernel wrapper for streaming limitation computation
void streamboundCuda (double* pos, gridnode_t* node, limit_t* limit, uint32_t M, uint32_t N);

/// Kernel wrapper for general preprocessing steps
void preprocessCuda (gridnode_t* node, double* pos, uint32_t* limit, 
  uint32_t maxlevel, uint32_t M, uint32_t DIM);

/// Kernel wrapper for standard evaluation
void evalCuda (double* res, double *a, gridnode_t* node, double* pos, uint32_t M,
  uint32_t maxlevel, uint32_t* limit, uint32_t subcnt, uint32_t* subs);
  
/// Kernel wrapper for transposed evaluation with additional FMA
void transposedCuda (double* a, gridnode_t* node, double* pos, double* y, 
  limit_t* limit, double* b, double c, uint32_t M, uint32_t N);
  
/// Kernel wrapper for transposed evaluation
void transposedCuda (double* a, gridnode_t* node, double* pos, double* y, 
  limit_t* limit, uint32_t M, uint32_t N);

#endif // MULTIEVALKERNEL_HPP
