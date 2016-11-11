// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#ifndef COMBIGRID_SRC_SGPP_COMBIGRID_STORAGE_TREE_TREESTORAGECONTEXT_HPP_
#define COMBIGRID_SRC_SGPP_COMBIGRID_STORAGE_TREE_TREESTORAGECONTEXT_HPP_

#include <sgpp/combigrid/definitions.hpp>

namespace sgpp {
namespace combigrid {

/**
 * Context for the TreeStorage class. Contains context information about the number of dimensions
 * and the function to generate default entries for the storage.
 * The context is referenced by all nodes.
 */
template <typename T>
class TreeStorageContext {
 public:
  typedef std::function<T(MultiIndex const &)> function_type;
  TreeStorageContext(size_t numDimensions, function_type func)
      : numDimensions(numDimensions), func(func) {}

  size_t numDimensions;
  function_type func;
};

} /* namespace combigrid */
} /* namespace sgpp */

#endif /* COMBIGRID_SRC_SGPP_COMBIGRID_STORAGE_TREE_TREESTORAGECONTEXT_HPP_ */
