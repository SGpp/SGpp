/* ****************************************************************************
 * Copyright (C) 2009 Technische Universitaet Muenchen                         *
 * This file is part of the SG++ project. For conditions of distribution and   *
 * use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
 **************************************************************************** */
// @author Jörg Blank (blankj@in.tum.de), Alexander Heinecke (Alexander.Heinecke@mytum.de)
#ifndef GRIDSTORAGE_HPP
#define GRIDSTORAGE_HPP

#include <sgpp/base/grid/storage/hashmap/HashGridIndex.hpp>
#include <sgpp/base/grid/storage/hashmap/HashGridStorage.hpp>
#include <sgpp/base/grid/storage/hashmap/HashGridIterator.hpp>

namespace sg {
namespace base {

typedef unsigned int level_t;
typedef unsigned int index_t;

/**
 * Main typedef for GridIndex
 */
//typedef HashGridIndex<unsigned int, unsigned int> GridIndex;
typedef HashGridIndex<level_t, index_t> GridIndex;

/**
 * Main typedef for GridStorage
 */
typedef HashGridStorage<GridIndex> GridStorage;

}
}

#endif /* GRIDSTORAGE_HPP */
