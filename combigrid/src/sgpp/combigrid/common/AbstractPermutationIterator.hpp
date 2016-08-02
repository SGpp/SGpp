/*
 * AbstractPermutationIterator.hpp
 *
 *  Created on: 12.12.2015
 *      Author: david
 */

#ifndef COMBIGRID_SRC_SGPP_COMBIGRID_COMMON_ABSTRACTPERMUTATIONITERATOR_HPP_
#define COMBIGRID_SRC_SGPP_COMBIGRID_COMMON_ABSTRACTPERMUTATIONITERATOR_HPP_

#include <sgpp/globaldef.hpp>
#include <cstddef>
#include <memory>

namespace SGPP {
namespace combigrid {

class AbstractPermutationIterator {
public:
	virtual ~AbstractPermutationIterator();

	/**
	 * Sets the iterator back to the start
	 */
	virtual void reset() = 0;

	virtual size_t value() = 0;

	virtual void moveToNext() = 0;

	virtual std::shared_ptr<AbstractPermutationIterator> clone() = 0;
};

} /* namespace combigrid */
} /* namespace SGPP */

#endif /* COMBIGRID_SRC_SGPP_COMBIGRID_COMMON_ABSTRACTPERMUTATIONITERATOR_HPP_ */
