// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

/*

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
//void* operator new (size_t size)
//// to ensure compatibility wit C++11
//#if __cplusplus < 201103L
//throw (std::bad_alloc)
//#endif
//;
//
///**
// * Overrides normal new[]
// *
// * @param size size of object
// */
//void* operator new[] (size_t size)
//// to ensure compatibility wit C++11
//#if __cplusplus < 201103L
//throw (std::bad_alloc)
//#endif
//;
///**
// * Overrides normal delete
// *
// * @param p pointer to data to free
// */
//void operator delete (void* p) throw();
//
///**
// * Overrides normal delete[]
// *
// * @param p pointer to data to free
// */
//void operator delete[] (void* p) throw();

/*
#endif
*/
