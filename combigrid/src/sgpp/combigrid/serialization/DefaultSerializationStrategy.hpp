/*
 * DefaultSerializationStrategy.hpp
 *
 *  Created on: 27.02.2016
 *      Author: david
 */

#ifndef COMBIGRID_SRC_SGPP_COMBIGRID_SERIALIZATION_DEFAULTSERIALIZATIONSTRATEGY_HPP_
#define COMBIGRID_SRC_SGPP_COMBIGRID_SERIALIZATION_DEFAULTSERIALIZATIONSTRATEGY_HPP_

#include "AbstractSerializationStrategy.hpp"
#include <sstream>

namespace SGPP {
namespace combigrid {

template<typename T> class DefaultSerializationStrategy: public AbstractSerializationStrategy<T> {
public:
	virtual ~DefaultSerializationStrategy(){}

	virtual std::string serialize(T const &value) {
		std::ostringstream stream;
		stream << value;
		return stream.str();
	}

	virtual T deserialize(std::string const &input) {
		std::istringstream stream(input);
		T value;
		stream >> value;
		return value;
	}
};

} /* namespace combigrid */
} /* namespace SGPP */

#endif /* COMBIGRID_SRC_SGPP_COMBIGRID_SERIALIZATION_DEFAULTSERIALIZATIONSTRATEGY_HPP_ */
