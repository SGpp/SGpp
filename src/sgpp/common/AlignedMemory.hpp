/******************************************************************************
* Copyright (C) 2010 Technische Universitaet Muenchen                         *
* This file is part of the SG++ project. For conditions of distribution and   *
* use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
******************************************************************************/
// @author Alexander Heinecke (Alexander.Heinecke@mytum.de)

#ifndef ALIGNEDMEMORY_HPP
#define ALIGNEDMEMORY_HPP

#include <new>
#include <exception>
#ifdef KNF
#include <stdlib.h>
#include <sys/types.h>
#include <sys/mman.h>
#undef aligned_malloc
#undef aligned_free
#define MAP_FAULTIN     0x80000000 /* fault in backing pages immediately */
#define MAP_PAGEORDER(n)      (((n) + 1 < 0xf ? (n) + 1 : 0xf) << 24)
#define MAP_PAGESIZE_4K       MAP_PAGEORDER(0)
#define MAP_PAGESIZE_64K      MAP_PAGEORDER(4)
#define MAP_PAGESIZE_2M       MAP_PAGEORDER(9)
#define aligned_malloc(size,alignment) \
  mmap(NULL, size, PROT_READ | PROT_WRITE, \
       MAP_PAGESIZE_4K | MAP_ANON | MAP_FAULTIN, -1, 0);
#define aligned_free(addr) munmap(addr, 0);
#else
#include <malloc.h>
#undef aligned_malloc
#undef aligned_free
#define aligned_malloc(size, alignment) memalign(alignment, size)
#define aligned_free(addr) free(addr)
#endif

/// define number of bytes that should used as alignment
#define SGPPMEMALIGNMENT 64

/**
 * Overrides normal new
 *
 * @param size size of object
 */
inline void* operator new (size_t size) throw (std::bad_alloc);

/**
 * Overrides normal new[]
 *
 * @param size size of object
 */
inline void* operator new[] (size_t size) throw (std::bad_alloc);

/**
 * Overrides normal delete
 *
 * @param p pointer to data to free
 */
inline void operator delete (void *p) throw();

/**
 * Overrides normal delete[]
 *
 * @param p pointer to data to free
 */
inline void operator delete[] (void *p) throw();

#endif /* ALIGNEDMEMORY_HPP */

