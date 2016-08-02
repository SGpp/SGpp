/*
 * SingleFunction.cpp
 *
 *  Created on: 27.02.2016
 *      Author: david
 */

#include "SingleFunction.hpp"

namespace SGPP {
namespace combigrid {

SingleFunction::SingleFunction(SGPP::float_t (*ptr)(SGPP::float_t)) : func(ptr) {
}

float_t SingleFunction::operator()(float_t param) {
	return func(param);
}

float_t SingleFunction::call(float_t param) {
	return func(param);
}

} /* namespace combigrid */
} /* namespace SGPP */
