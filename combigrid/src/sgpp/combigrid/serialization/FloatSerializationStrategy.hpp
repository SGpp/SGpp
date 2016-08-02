/*
 * FloatSerializationStrategy.hpp
 *
 *  Created on: 26.01.2016
 *      Author: david
 */

#ifndef COMBIGRID_SRC_SGPP_COMBIGRID_SERIALIZATION_FLOATSERIALIZATIONSTRATEGY_HPP_
#define COMBIGRID_SRC_SGPP_COMBIGRID_SERIALIZATION_FLOATSERIALIZATIONSTRATEGY_HPP_

#include "AbstractSerializationStrategy.hpp"
#include <cmath>
#include <sstream>
#include <limits>
#include <sgpp/globaldef.hpp>

namespace SGPP {
namespace combigrid {

template<typename T> class FloatSerializationStrategy : public AbstractSerializationStrategy<T> {
	static const int64_t SHIFT_AMOUNT = 58;

	static const char normalPrefix = 'd';
	static const char infPrefix = 'u';
	static const char minusInfPrefix = 'l';
	static const char nanPrefix = 'n';

public:
	virtual ~FloatSerializationStrategy(){}

	// TODO: handle nan and inf
	virtual std::string serialize(T const &value) {
		if(std::isnan(value)) {
			return std::string(1, nanPrefix);
		} else if(std::isinf(value)) {
			if(value < 0) {
				return std::string(1, minusInfPrefix);
			} else {
				return std::string(1, infPrefix);
			}
		}

		int exponent = 0;

		T normalized = frexp(value, &exponent);

		int64_t shifted = static_cast<int64_t>(normalized * static_cast<T>(static_cast<int64_t>(1) << SHIFT_AMOUNT));

		std::ostringstream str;
		str << shifted << "e" << exponent;

		return normalPrefix + str.str();
	}

	virtual T deserialize(std::string const &input) {
		if(input[0] == infPrefix) {
			return std::numeric_limits<T>::infinity();
		} else if(input[0] == minusInfPrefix) {
			return -std::numeric_limits<T>::infinity();
		} else if(input[0] == nanPrefix) {
			return std::numeric_limits<T>::quiet_NaN();
		}

		std::istringstream str(input.substr(1));

		int exponent;
		int64_t shifted;

		str >> shifted;
		str.get();
		str >> exponent;

		int64_t one = 1;

		if(exponent >= 0) {
			return (static_cast<T>(shifted) / static_cast<T>(one << SHIFT_AMOUNT)) * static_cast<T>(one << exponent);
		} else {
			return (static_cast<T>(shifted) / static_cast<T>(one << SHIFT_AMOUNT)) / static_cast<T>(one << (-exponent));
		}
	}
};

} /* namespace combigrid */
} /* namespace SGPP */

#endif /* COMBIGRID_SRC_SGPP_COMBIGRID_SERIALIZATION_FLOATSERIALIZATIONSTRATEGY_HPP_ */
