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
#include <lmmintrin.h>
#undef aligned_malloc
#undef aligned_free
#define aligned_malloc(size, alignment) _mm_malloc(size, alignment)
#define aligned_free(addr) _mm_free(addr)
#else
#ifdef WINDOWS
#include <pmmintrin.h>
#undef aligned_malloc
#undef aligned_free
#define aligned_malloc(size, alignment) _mm_malloc(size, alignment)
#define aligned_free(addr) _mm_free(addr)
#else
#include <malloc.h>
#undef aligned_malloc
#undef aligned_free
#define aligned_malloc(size, alignment) memalign(alignment, size)
#define aligned_free(addr) free(addr)
#endif
#endif

/// define number of bytes that should used as alignment
#define SGPPMEMALIGNMENT 64

/**
 * Overrides normal new
 *
 * @param size size of object
 */
void* operator new (size_t size) throw (std::bad_alloc);

/**
 * Overrides normal new[]
 *
 * @param size size of object
 */
void* operator new[] (size_t size) throw (std::bad_alloc);

/**
 * Overrides normal delete
 *
 * @param p pointer to data to free
 */
void operator delete (void *p) throw();

/**
 * Overrides normal delete[]
 *
 * @param p pointer to data to free
 */
void operator delete[] (void *p) throw();

#endif /* ALIGNEDMEMORY_HPP */

