/*
 * AbstractSerializationStrategy.hpp
 *
 *  Created on: 26.01.2016
 *      Author: david
 */

#ifndef COMBIGRID_SRC_SGPP_COMBIGRID_SERIALIZATION_ABSTRACTSERIALIZATIONSTRATEGY_HPP_
#define COMBIGRID_SRC_SGPP_COMBIGRID_SERIALIZATION_ABSTRACTSERIALIZATIONSTRATEGY_HPP_

#include <string>
#include <sgpp/globaldef.hpp>

namespace sgpp{
namespace combigrid {

template<typename T> class AbstractSerializationStrategy {
public:
	virtual ~AbstractSerializationStrategy(){}

	virtual std::string serialize(T const &value) = 0;
	virtual T deserialize(std::string const &input) = 0;
};

} /* namespace combigrid */
} /* namespace sgpp*/

#endif /* COMBIGRID_SRC_SGPP_COMBIGRID_SERIALIZATION_ABSTRACTSERIALIZATIONSTRATEGY_HPP_ */
