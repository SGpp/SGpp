/******************************************************************************
* Copyright (C) 2010 Technische Universitaet Muenchen                         *
* This file is part of the SG++ project. For conditions of distribution and   *
* use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
******************************************************************************/
// @author Alexander Heinecke (Alexander.Heinecke@mytum.de)

#include "base/tools/AlignedMemory.hpp"

#if defined(__GXX_EXPERIMENTAL_CXX0X__) && \
		(__GNUC__ == 4 && ((__GNUC_MINOR__ == 7) || (__GNUC_MINOR__ == 8)))
// g++ with C++11 enabled (at least 4.7 and 4.8) seem to have a different exception
// specifier for "new"
void* operator new (size_t size) _GLIBCXX_THROW (std::bad_alloc) {
#else
void* operator new (size_t size) throw (std::bad_alloc) {
#endif
  void* p;

  //Workaround for apples non-standard implementation of posix_memalign().
  //If a size of 0 is requested, posix_memalign() should return a pointer which
  //free() can be called successfully with. Unfortunately this is not possible 
  //on OS X.
  #ifdef __APPLE__
  if (size == 0) {
    size = 1;
  }
  #endif
  
  //p = aligned_malloc(size, SGPPMEMALIGNMENT);
  aligned_malloc(p, size, SGPPMEMALIGNMENT);

  #ifdef POSIX_MEMALIGN
  if (success != 0) {
    throw std::bad_alloc();
  }
  #else
  if (p == 0) {
    throw std::bad_alloc();
  }
  #endif

  return p;
}

#if defined(__GXX_EXPERIMENTAL_CXX0X__) && \
		(__GNUC__ == 4 && ((__GNUC_MINOR__ == 7) || (__GNUC_MINOR__ == 8)))
// g++ with C++11 enabled (at least 4.7 and 4.8) seem to have a different exception
// specifier for "new"
void* operator new[] (size_t size) _GLIBCXX_THROW (std::bad_alloc) {
#else
void* operator new[] (size_t size) throw (std::bad_alloc) {
#endif
  void* p;

  //Workaround for apples non-standard implementation of posix_memalign().
  //If a size of 0 is requested, posix_memalign() should return a pointer which
  //free() can be called successfully with. Unfortunately this is not possible 
  //on OS X.
  //long unsigned int oldSize = size;
  #ifdef __APPLE__
  if (size == 0) {
    size = 64;
  }
  #endif

  //p = aligned_malloc(size, SGPPMEMALIGNMENT);
  aligned_malloc(p, size, SGPPMEMALIGNMENT);

  #ifdef POSIX_MEMALIGN
  if (success != 0) {
    throw std::bad_alloc();
  }
  #else
  if (p == 0) {
    throw std::bad_alloc();
  }
  #endif

  return p;
}

void operator delete (void* p) throw () {
  aligned_free(p);
}

void operator delete[] (void* p) throw() {
  aligned_free(p);
}
