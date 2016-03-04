// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <sgpp/base/grid/generation/hashmap/HashRefinementInconsistent.hpp>
#include <sgpp/base/exception/generation_exception.hpp>

#include <sgpp/globaldef.hpp>

#include <limits>


namespace sgpp {
namespace base {


void HashRefinementInconsistent::createGridpoint(GridStorage& storage,
    index_type& index) {
  storage.insert(index);
}

}  // namespace base
}  // namespace sgpp
