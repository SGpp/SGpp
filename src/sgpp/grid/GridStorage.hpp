/******************************************************************************
* Copyright (C) 2009 Technische Universitaet Muenchen                         *
* This file is part of the SG++ project. For conditions of distribution and   *
* use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
******************************************************************************/
// @author Jörg Blank (blankj@in.tum.de), Alexander Heinecke (Alexander.Heinecke@mytum.de)

#ifndef GRIDSTORAGE_HPP
#define GRIDSTORAGE_HPP

#include "grid/storage/hashmap/HashGridIndex.hpp"
#include "grid/storage/hashmap/HashGridStorage.hpp"
#include "grid/storage/hashmap/HashGridIterator.hpp"

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
