/*
 * AbstractTreeStorageNode.hpp
 *
 *  Created on: 12.12.2015
 *      Author: david
 */

#ifndef COMBIGRID_SRC_SGPP_COMBIGRID_STORAGE_TREE_ABSTRACTTREESTORAGENODE_HPP_
#define COMBIGRID_SRC_SGPP_COMBIGRID_STORAGE_TREE_ABSTRACTTREESTORAGENODE_HPP_

#include <cstddef>
#include <memory>
#include "../../definitions.hpp"

namespace sgpp{
namespace combigrid {

enum class StorageStatus : uint8_t {
	NOT_STORED = 0,
	REQUESTED, // A thread has been started to compute the value, but it is not yet available
	STORED
};

/**
 * Interface class for TreeStorage-nodes
 */
template<typename T> class AbstractTreeStorageNode {
public:
	virtual ~AbstractTreeStorageNode(){}

	/**
	 * @return Returns the number of children
	 */
	virtual size_t numChildren() const = 0;

	/**
	 * @return Returns true iff the node is a leaf
	 */
	virtual bool isLeaf() const = 0;

	/**
	 * @param depth Level of the node in the tree, starting from zero.
	 * @return Returns a reference to the storage entry for the given multi-index.
	 * If the entry does not already exist, it is created and initialized with the value resulting from calling the default-value-function.
	 */
	virtual T &get(MultiIndex const &index, size_t depth = 0) = 0;

	/**
	 * Stores the given value at the given multi-index. If the storage entry does not exist, it is created. The default-value-function is not called.
	 * @param depth Level of the node in the tree, starting from zero
	 */
	virtual void set(MultiIndex const &index, T const &value, size_t depth = 0) = 0;
	virtual bool containsIndex(MultiIndex const &index, size_t depth = 0) const = 0;
};

} /* namespace combigrid */
} /* namespace sgpp*/

#endif /* COMBIGRID_SRC_SGPP_COMBIGRID_STORAGE_TREE_ABSTRACTTREESTORAGENODE_HPP_ */
