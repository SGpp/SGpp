#ifndef MORTONORDERKERNEL_HPP
#define MORTONORDERKERNEL_HPP

#include <stdint.h>

/// Kernel wrapper for the Morton order generator
void zorder (double* pos, size_t* perm, uint32_t m, uint32_t DIM);

#endif // MORTONORDERKERNEL_HPP
