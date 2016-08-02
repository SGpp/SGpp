/*
 * MultiFunction.cpp
 *
 *  Created on: 13.01.2016
 *      Author: david
 */

#include "MultiFunction.hpp"

namespace sgpp{
namespace combigrid {

MultiFunction::MultiFunction(double (*ptr)(const base::DataVector&)) : func(ptr) {
}

double MultiFunction::operator ()(const base::DataVector& vec) {
	return func(vec);
}

double MultiFunction::call(const base::DataVector& vec) {
	return func(vec);
}

} /* namespace combigrid */
} /* namespace sgpp*/
