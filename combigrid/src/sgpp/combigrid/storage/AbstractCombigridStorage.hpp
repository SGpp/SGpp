/*
 * AbstractCombigridStorage.hpp
 *
 *  Created on: 29.12.2015
 *      Author: david
 */

#ifndef COMBIGRID_SRC_SGPP_COMBIGRID_STORAGE_ABSTRACTCOMBIGRIDSTORAGE_HPP_
#define COMBIGRID_SRC_SGPP_COMBIGRID_STORAGE_ABSTRACTCOMBIGRIDSTORAGE_HPP_

#include <sgpp/combigrid/definitions.hpp>
#include <sgpp/combigrid/common/MultiIndexIterator.hpp>
#include "AbstractMultiStorageIterator.hpp"
#include <memory>

namespace SGPP {
namespace combigrid {

class AbstractCombigridStorage {
public:
	virtual ~AbstractCombigridStorage();

	/**
	 * @param level Level to traverse
	 * @param iterator Iterator that defines which values should be iterated over
	 * @param orderingConfiguration Defines for each dimension whether the grid points in that dimension should be traversed in sorted order.
	 * @return Returns an iterator that iterates along the multi-indices for a given level. These multi-indices are given by the parameter iterator.
	 * When moveToNext() is called on the returned iterator, it also calls moveToNext() on the underlying MultiIndexIterator.
	 * If the values are not already stored, they are created during iteration.
	 */
	virtual std::shared_ptr<AbstractMultiStorageIterator<SGPP::float_t> > getGuidedIterator(MultiIndex const &level, MultiIndexIterator &iterator,
			std::vector<bool> orderingConfiguration) = 0;

	/**
	 * @return Returns the number of stored values (cumulated over all levels).
	 */
	virtual size_t getNumEntries() = 0;

	virtual std::string serialize() = 0;
	virtual void deserialize(std::string const &str) = 0;
	virtual void set(MultiIndex const &level, MultiIndex const &index, float_t value) = 0;
};

} /* namespace combigrid */
} /* namespace SGPP */

#endif /* COMBIGRID_SRC_SGPP_COMBIGRID_STORAGE_ABSTRACTCOMBIGRIDSTORAGE_HPP_ */
