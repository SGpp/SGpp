/*
 * SingleFunction.cpp
 *
 *  Created on: 27.02.2016
 *      Author: david
 */

#include "SingleFunction.hpp"

namespace sgpp{
namespace combigrid {

SingleFunction::SingleFunction(double (*ptr)(double)) : func(ptr) {
}

double SingleFunction::operator()(double param) {
	return func(param);
}

double SingleFunction::call(double param) {
	return func(param);
}

} /* namespace combigrid */
} /* namespace sgpp*/
