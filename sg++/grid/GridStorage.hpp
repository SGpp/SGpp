/*****************************************************************************/
/* This file is part of sg++, a program package making use of spatially      */
/* adaptive sparse grids to solve numerical problems                         */
/*                                                                           */
/* Copyright (C) 2008 Jörg Blank (blankj@in.tum.de)                          */
/* Copyright (C) 2009 Alexander Heinecke (Alexander.Heinecke@mytum.de)       */
/*                                                                           */
/* sg++ is free software; you can redistribute it and/or modify              */
/* it under the terms of the GNU General Public License as published by      */
/* the Free Software Foundation; either version 3 of the License, or         */
/* (at your option) any later version.                                       */
/*                                                                           */
/* sg++ is distributed in the hope that it will be useful,                   */
/* but WITHOUT ANY WARRANTY; without even the implied warranty of            */
/* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the             */
/* GNU General Public License for more details.                              */
/*                                                                           */
/* You should have received a copy of the GNU General Public License         */
/* along with sg++; if not, write to the Free Software                       */
/* Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA */
/* or see <http://www.gnu.org/licenses/>.                                    */
/*****************************************************************************/

#ifndef GRIDSTORAGE_HPP
#define GRIDSTORAGE_HPP

#include "grid/storage/hashmap/HashGridIndex.hpp"
#include "grid/storage/hashmap/HashGridStorage.hpp"
#include "grid/storage/hashmap/HashGridIterator.hpp"

#define SERIALIZATION_VERSION 1

namespace sg {

/**
 * Main typedef for GridIndex
 */
typedef HashGridIndex<unsigned int, unsigned int> GridIndex;
/**
 * Main typedef for GridStorage
 */
typedef HashGridStorage<GridIndex> GridStorage;

}

#endif /* GRIDSTORAGE_HPP */
