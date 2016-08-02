/*
 * MultiFunction.cpp
 *
 *  Created on: 13.01.2016
 *      Author: david
 */

#include "MultiFunction.hpp"

namespace SGPP {
namespace combigrid {

MultiFunction::MultiFunction(SGPP::float_t (*ptr)(const base::DataVector&)) : func(ptr) {
}

SGPP::float_t MultiFunction::operator ()(const base::DataVector& vec) {
	return func(vec);
}

SGPP::float_t MultiFunction::call(const base::DataVector& vec) {
	return func(vec);
}

} /* namespace combigrid */
} /* namespace SGPP */
