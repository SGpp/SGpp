/*****************************************************************************/
/* This file is part of sgpp, a program package making use of spatially      */
/* adaptive sparse grids to solve numerical problems                         */
/*                                                                           */
/* Copyright (C) 2008 Jörg Blank (blankj@in.tum.de)                          */
/* Copyright (C) 2009-2010 Alexander Heinecke (Alexander.Heinecke@mytum.de)  */
/*                                                                           */
/* sgpp is free software; you can redistribute it and/or modify              */
/* it under the terms of the GNU Lesser General Public License as published  */
/* by the Free Software Foundation; either version 3 of the License, or      */
/* (at your option) any later version.                                       */
/*                                                                           */
/* sgpp is distributed in the hope that it will be useful,                   */
/* but WITHOUT ANY WARRANTY; without even the implied warranty of            */
/* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the             */
/* GNU Lesser General Public License for more details.                       */
/*                                                                           */
/* You should have received a copy of the GNU Lesser General Public License  */
/* along with sgpp; if not, write to the Free Software                       */
/* Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA */
/* or see <http://www.gnu.org/licenses/>.                                    */
/*****************************************************************************/

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

// if available you can use the upcoming standard: unorded_map
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

