/******************************************************************************
* Copyright (C) 2010 Technische Universitaet Muenchen                         *
* This file is part of the SG++ project. For conditions of distribution and   *
* use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
******************************************************************************/
// @author Alexander Heinecke (Alexander.Heinecke@mytum.de)

#include "base/tools/AlignedMemory.hpp"

void* operator new (size_t size) throw (std::bad_alloc) {
  void* p;
  p = aligned_malloc(size, SGPPMEMALIGNMENT);

  if (p == 0) {
    throw std::bad_alloc();
  }

  return p;
}

void* operator new[] (size_t size) throw (std::bad_alloc) {
  void* p;
  p = aligned_malloc(size, SGPPMEMALIGNMENT);

  if (p == 0) {
    throw std::bad_alloc();
  }

  return p;
}

void operator delete (void* p) throw () {
  aligned_free(p);
}

void operator delete[] (void* p) throw() {
  aligned_free(p);
}
