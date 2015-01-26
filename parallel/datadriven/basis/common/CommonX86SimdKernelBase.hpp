/* ****************************************************************************
* Copyright (C) 2013 Technische Universitaet Muenchen                         *
* This file is part of the SG++ project. For conditions of distribution and   *
* use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
**************************************************************************** */
// @author Roman Karlstetter (karlstetter@mytum.de)

#ifndef KERNELMACROS_HPP
#define KERNELMACROS_HPP

#include "parallel/datadriven/basis/common/KernelBase.hpp"

#if defined(__SSE3__) || defined(__AVX__)
#include <immintrin.h>
#endif
#if defined(__FMA4__)
#include <x86intrin.h>
#endif

#ifdef __USEAVX128__
#undef __AVX__
#endif


#endif // KERNELMACROS_HPP
