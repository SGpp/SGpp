// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

// In this file, the operators "new" and "delete" (and their array
// counterparts are overriden).
// This file does not require a header as the overloaded functions are
// implicitly defined.
// The "overload" happens at link time of the library.

// TODO(valentjn): On MinGW, using aligned memory with _mm_malloc
// leads to crashes (e.g., in the Boost tests). posix_memalign isn't defined
// on MinGW. Somebody should enable aligned memory for MinGW...
#ifndef __MINGW64__

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
#else  // apple or linux or unknown
#include <stdlib.h>
#undef aligned_malloc
#undef aligned_free
#define POSIX_MEMALIGN
#define aligned_malloc(p, size, alignment) \
  int success = posix_memalign(&p, SGPPMEMALIGNMENT, size);
#define aligned_free(addr) free(addr)
#endif /* _WIN32 */

void* operator new(size_t size) {
  void* p;

  // Workaround for apples non-standard implementation of posix_memalign().
  // If a size of 0 is requested, posix_memalign() should return a pointer which
  // free() can be called successfully with. Unfortunately this is not possible
  // on OS X.
#ifdef __APPLE__

  if (size == 0) {
    size = 1;
  }

#endif

  // p = aligned_malloc(size, SGPPMEMALIGNMENT);
  aligned_malloc(p, size, SGPPMEMALIGNMENT);

#ifdef POSIX_MEMALIGN

  if (success != 0) {
    throw std::bad_alloc();
  }

#else /* POSIX_MEMALIGN */

  if (p == 0) {
    throw std::bad_alloc();
  }

#endif /* POSIX_MEMALIGN */

  return p;
}

void* operator new[](size_t size)
// to ensure compatibility wit C++11
#if __cplusplus < 201103L
throw(std::bad_alloc)
#endif
{
  void* p;

  // Workaround for apples non-standard implementation of posix_memalign().
  // If a size of 0 is requested, posix_memalign() should return a pointer which
  // free() can be called successfully with. Unfortunately this is not possible
  // on OS X.
  // long unsigned int oldSize = size;
#ifdef __APPLE__

  if (size == 0) {
    size = 64;
  }

#endif /* __APPLE__ */

  // p = aligned_malloc(size, SGPPMEMALIGNMENT);
  aligned_malloc(p, size, SGPPMEMALIGNMENT);

#ifdef POSIX_MEMALIGN

  if (success != 0) {
    throw std::bad_alloc();
  }

#else /* POSIX_MEMALIGN */

  if (p == 0) {
    throw std::bad_alloc();
  }

#endif /* POSIX_MEMALIGN */

  return p;
}

void operator delete(void* p) noexcept {
  aligned_free(p);
}

#if __cplusplus >= 201402L
void operator delete(void* p, size_t sz) noexcept {
  aligned_free(p);
}
#endif

void operator delete[](void* p) noexcept {
  aligned_free(p);
}

#if __cplusplus >= 201402L
void operator delete[](void* p, size_t sz) noexcept {
  aligned_free(p);
}
#endif

#endif /* __MINGW64__ */
