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
#define memalign _mm_malloc
#define freealign _mm_free
#else
#include <malloc.h>
#define freealign free
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

