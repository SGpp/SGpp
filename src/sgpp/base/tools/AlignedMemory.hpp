/* ****************************************************************************
* Copyright (C) 2010 Technische Universitaet Muenchen                         *
* This file is part of the SG++ project. For conditions of distribution and   *
* use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
**************************************************************************** */
// @author Alexander Heinecke (Alexander.Heinecke@mytum.de), David Pfander (David.Pfander@ipvs.uni-stuttgart.de)

#ifndef ALIGNEDMEMORY_HPP
#define ALIGNEDMEMORY_HPP

#include <new>
#include <exception>

// define number of bytes that should used as alignment
#define SGPPMEMALIGNMENT 64

#ifdef _WIN32
#include <pmmintrin.h>
#undef aligned_malloc
#undef aligned_free
#define aligned_malloc(p, size, alignment) p = _mm_malloc(size, alignment)
#define aligned_free(addr) _mm_free(addr)
#else //apple or linux or unknown
#include <stdlib.h>
#undef aligned_malloc
#undef aligned_free
#define POSIX_MEMALIGN
#define aligned_malloc(p, size, alignment) int success = posix_memalign(&p, SGPPMEMALIGNMENT, size);  
#define aligned_free(addr) free(addr)
#endif

/**
 * Overrides normal new
 *
 * @param size size of object
 */
#if defined(__GXX_EXPERIMENTAL_CXX0X__) && \
		(__GNUC__ == 4 && ((__GNUC_MINOR__ == 7) || (__GNUC_MINOR__ == 8)))
// g++ with C++11 enabled (at least 4.7 and 4.8) seem to have a different exception
// specifier for "new"
void* operator new (size_t size) _GLIBCXX_THROW (std::bad_alloc);
#else
void* operator new (size_t size) throw (std::bad_alloc);
#endif

/**
 * Overrides normal new[]
 *
 * @param size size of object
 */
#if defined(__GXX_EXPERIMENTAL_CXX0X__) && \
		(__GNUC__ == 4 && ((__GNUC_MINOR__ == 7) || (__GNUC_MINOR__ == 8)))
// g++ with C++11 enabled (at least 4.7 and 4.8) seem to have a different exception
// specifier for "new"
void* operator new[] (size_t size) _GLIBCXX_THROW (std::bad_alloc);
#else
void* operator new[] (size_t size) throw (std::bad_alloc);
#endif

/**
 * Overrides normal delete
 *
 * @param p pointer to data to free
 */
void operator delete (void* p) throw();

/**
 * Overrides normal delete[]
 *
 * @param p pointer to data to free
 */
void operator delete[] (void* p) throw();

#endif /* ALIGNEDMEMORY_HPP */

