/*
This file is part of sg++, a program package making use of spatially adaptive sparse grids to solve numerical problems

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


#ifndef SGPP_HPP_
#define SGPP_HPP_

// Optimizations currently useless
//#define SGPP_OPTIMIZE

#include "GridIndex.hpp"
#include "GridStorage.hpp"
#include "algorithms.hpp"
#include "generation.hpp"

#include "basis/basis.hpp"

#include "basis/linear/operation/OperationLaplaceLinear.hpp"
#include "basis/modlinear/operation/OperationLaplaceModLinear.hpp"

#include "data/DataVector.h"

#include "GridFactory.hpp"
#include "Operations.hpp"

namespace sg
{

typedef linear_base<unsigned int, unsigned int> SLinearBase;
typedef modified_linear_base<unsigned int, unsigned int> SModLinearBase;
typedef poly_base<unsigned int, unsigned int> SPolyBase;
typedef modified_poly_base<unsigned int, unsigned int> SModPolyBase;

typedef AlgorithmB<SLinearBase> SGridOperationB;
typedef AlgorithmB<SModLinearBase> SGridModOperationB;

}


/*
typedef sg::GridIndex<unsigned int, unsigned int> SGridIndex;
typedef sg::GridStorage<SGridIndex> SGridStorage;
*/

#endif /*SGPP_HPP_*/
