// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at 
// sgpp.sparsegrids.org

#ifndef HASH_MAP_CONFIG
#define HASH_MAP_CONFIG

/*
// backward compatible: able to use the  standard gnu hashmap of linux (SGI/STLPort)
#ifndef USETRONE
#include <ext/hash_map>
namespace std {
  using namespace __gnu_cxx;
}
#endif

// if available you can use the upcoming standard: unordered_map
#ifdef USETRONE

// do some defines for icc, avoiding, icc only with gcc 4.3 or lower
// errors:
// See:
// http://software.intel.com/en-us/forums/intel-c-compiler/topic/65041/
#ifndef _WIN32
#define __aligned__   ignored
#include <tr1/unordered_map>
#undef __aligned__
#else
//#define __aligned__   ignored
#include <unordered_map>
//#undef __aligned__
#endif

#endif
*/

// forward declaration of hash function and hash comparison function
#include <sgpp/globaldef.hpp>


namespace SGPP {
  namespace base {
    template<class key>
    struct hash { };

    template<class key>
    struct eqIndex { };
  }
}

#endif /* HASH_MAP_CONFIG */
