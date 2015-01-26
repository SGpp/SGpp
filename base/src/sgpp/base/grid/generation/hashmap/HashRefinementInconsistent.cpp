/******************************************************************************
 * Copyright (C) 2014 Technische Universitaet Muenchen                         *
 * This file is part of the SG++ project. For conditions of distribution and   *
 * use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
 ******************************************************************************/
// @author Valeriy Khakhutskyy (khakhutv@in.tum.de)

#include <limits>

#include "base/grid/generation/hashmap/HashRefinementInconsistent.hpp"
#include "base/exception/generation_exception.hpp"


namespace sg {
namespace base {


void HashRefinementInconsistent::createGridpoint(GridStorage* storage, index_type& index)
{
    storage->insert(index);
}

}
}

