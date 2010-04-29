/******************************************************************************
* Copyright (C) 2009 Technische Universitaet Muenchen                         *
* This file is part of the SG++ project. For conditions of distribution and   *
* use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
******************************************************************************/
// @author Jörg Blank (blankj@in.tum.de), Alexander Heinecke (Alexander.Heinecke@mytum.de)

#ifndef HASH_MAP_CONFIG
#define HASH_MAP_CONFIG

// use the MS hashmap (TR1 is not supported, yet) on Larrabee
#ifdef LARRABEENATIVE
#include <ext/hash_map>
namespace sg {
	template<class key>
	class LRBSGHasher;
}
#endif

// backward compatible: able to use the  standard gnu hashmap of linux (SGI/STLPort)
#ifndef LARRABEENATIVE
#ifndef USETRONE
#include <ext/hash_map>
namespace std {
    using namespace __gnu_cxx;
}
#endif
#endif

// if available you can use the upcoming standard: unordered_map
#ifdef USETRONE

// do some defines for icc, avoiding
// errors:
// See:
// http://software.intel.com/en-us/forums/intel-c-compiler/topic/65041/
#ifndef WINDOWS
#ifndef AIX_XLC
#define __aligned__   ignored
#include <tr1/unordered_map>
#undef __aligned__
#endif
#endif

#ifdef WINDOWS
#define __aligned__   ignored
#include <unordered_map>
#undef __aligned__
#endif

#ifdef AIX_XLC
#define __IBMCPP_TR1__ 1
#include <unordered_map>
#endif

#endif

// forward declaration of hash function and hash comparison function
#ifndef LARRABEENATIVE
namespace sg {
	template<class key>
	struct hash { };

	template<class key>
	struct eqIndex { };
}
#endif

#endif /* HASH_MAP_CONFIG */

