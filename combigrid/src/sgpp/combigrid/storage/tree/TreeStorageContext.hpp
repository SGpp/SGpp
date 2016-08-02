/*
 * TreeStorageContext.hpp
 *
 *  Created on: 12.12.2015
 *      Author: david
 */

#ifndef COMBIGRID_SRC_SGPP_COMBIGRID_STORAGE_TREE_TREESTORAGECONTEXT_HPP_
#define COMBIGRID_SRC_SGPP_COMBIGRID_STORAGE_TREE_TREESTORAGECONTEXT_HPP_

#include "../../definitions.hpp"

#include <mutex>

namespace SGPP {
namespace combigrid {

/**
 * Context for the TreeStorage class. Contains context information about the number of dimensions and the function to generate default entries for the storage.
 * The context is referenced by all nodes.
 */
template<typename T> class TreeStorageContext {
public:
	typedef std::function<T(MultiIndex const &)> function_type;
	TreeStorageContext(size_t numDimensions, function_type func) : numDimensions(numDimensions), func(func) {

	}

	size_t numDimensions;
	function_type func;
};

} /* namespace combigrid */
} /* namespace SGPP */

#endif /* COMBIGRID_SRC_SGPP_COMBIGRID_STORAGE_TREE_TREESTORAGECONTEXT_HPP_ */
