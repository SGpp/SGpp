/*
This file is part of sg++, a program package making use of spatially adaptive sparse grids to solve numerical problems

Copyright (C) 2008 Joerg Blank (blankj@in.tum.de)

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

#include "grid/generation/StandardGridGenerator.hpp"
#include "grid/GridStorage.hpp"

#include "sgpp.hpp"

namespace sg
{

StandardGridGenerator::StandardGridGenerator(GridStorage* storage) : storage(storage)
{
}

StandardGridGenerator::~StandardGridGenerator()
{
}

void StandardGridGenerator::regular(size_t level)
{
	HashGenerator gen;
	gen.regular(this->storage, level);
}

void StandardGridGenerator::refine(RefinementFunctor* func)
{
	HashRefinement refine;
	refine.free_refine(this->storage, func);
}

}
