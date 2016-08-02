/*
 * InternalTreeStorageNode.hpp
 *
 *  Created on: 12.12.2015
 *      Author: david
 */

#ifndef COMBIGRID_SRC_SGPP_COMBIGRID_STORAGE_TREE_INTERNALTREESTORAGENODE_HPP_
#define COMBIGRID_SRC_SGPP_COMBIGRID_STORAGE_TREE_INTERNALTREESTORAGENODE_HPP_

#include "AbstractTreeStorageNode.hpp"
#include "LowestTreeStorageNode.hpp"
#include <vector>

namespace SGPP {
namespace combigrid {

/**
 * Each inner node of the storage contains a vector of other nodes.
 * For a d-dimensional multi-index, to find a storage entry, the components of the multi-index are iterated over.
 * Each component defines the index of the next child to be traversed.
 * Because traversal would be slower if each leaf only contained one number,
 * the nodes for the last dimension contain a vector of numbers (or Ts) instead of having children which themselves have one number.
 * These nodes are represented by LowestTreeStorageNode objects, all other nodes are InternalTreeStorageNode objects.
 */
template<typename T> class InternalTreeStorageNode: public AbstractTreeStorageNode<T> {
	MultiIndex restrictIndex(MultiIndex const &index, size_t depth) {
		MultiIndex result(depth + 1);
		for (size_t i = 0; i <= depth; ++i) {
			result[i] = index[i];
		}
		return result;
	}

public:
	std::vector<std::unique_ptr<AbstractTreeStorageNode<T>>> children;
	TreeStorageContext<T> &context;

	InternalTreeStorageNode(TreeStorageContext<T> &context, size_t remainingDimensions, MultiIndex const &index, size_t depth, size_t zeroDepth, bool createFirstChild = true)
		: children(), context(context) {
			if(!createFirstChild) {
				return;
			}

			if(remainingDimensions > 1) {
				children.emplace_back(new InternalTreeStorageNode<T>(context, remainingDimensions - 1, index, depth+1, zeroDepth, createFirstChild));
			} else {
				children.emplace_back(new LowestTreeStorageNode<T>(context, index, depth+1, zeroDepth, createFirstChild));
			}
		}

	virtual ~InternalTreeStorageNode() {}

	virtual size_t numChildren() const {
		return children.size();
	}

	virtual bool isLeaf() const {
		return false;
	}

	virtual T &get(MultiIndex const &index, size_t depth = 0) {
		size_t currentIndex = index[depth];
		size_t remainingDimensions = context.numDimensions - depth - 1;

		while (currentIndex >= children.size()) {
			MultiIndex restrictedIndex = restrictIndex(index, depth);
			restrictedIndex[depth] = children.size();
			bool createFirstChild = true;
			if(remainingDimensions >= 2) {
				children.emplace_back(new InternalTreeStorageNode<T>(context, remainingDimensions-1, index, depth+1, depth, createFirstChild));
			} else {
				children.emplace_back(new LowestTreeStorageNode<T>(context, index, depth+1, depth, createFirstChild));
			}
		}
		return children[currentIndex]->get(index, depth+1);
	}

	virtual void set(MultiIndex const &index, T const &value, size_t depth = 0) {
		size_t currentIndex = index[depth];
		size_t remainingDimensions = context.numDimensions - depth - 1;

		while (currentIndex >= children.size()) {
			MultiIndex restrictedIndex = restrictIndex(index, depth);
			restrictedIndex[depth] = children.size();
			bool createFirstChild = (currentIndex != children.size()); // if this is the child to set, then it will be created in the recursion or by setting the value
			if(remainingDimensions >= 2) {
				children.emplace_back(new InternalTreeStorageNode<T>(context, remainingDimensions-1, index, depth+1, depth, createFirstChild));
			} else {
				children.emplace_back(new LowestTreeStorageNode<T>(context, index, depth, depth+1, createFirstChild));
			}
		}

		children[currentIndex]->set(index, value, depth+1);
	}

	virtual bool containsIndex(MultiIndex const &index, size_t depth = 0) const {
		return index[depth] < children.size() && children[index[depth]]->containsIndex(index, depth+1);
	}
};

}
/* namespace combigrid */
} /* namespace SGPP */

#endif /* COMBIGRID_SRC_SGPP_COMBIGRID_STORAGE_TREE_INTERNALTREESTORAGENODE_HPP_ */
