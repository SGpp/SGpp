/*
This file is part of sgpp, a program package making use of spatially adaptive sparse grids to solve numerical problems

Copyright (C) 2007  Jörg Blank (blankj@in.tum.de)

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU Lesser General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU Lesser General Public License for more details.

You should have received a copy of the GNU Lesser General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/

#ifndef HASH_MAP_CONFIG
#define HASH_MAP_CONFIG

#ifndef WINDOWS
#include <ext/hash_map>
namespace std {

    using namespace __gnu_cxx;

}
#endif
#ifdef WINDOWS
#include <hash_map>
#endif

namespace sg {
	
	template<class key>
	struct hash { };
	
	template<class key>
	struct eqIndex { };
	
}

#endif
