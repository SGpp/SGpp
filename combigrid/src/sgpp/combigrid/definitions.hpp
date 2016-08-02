/*
 * definitions.hpp
 *
 *  Created on: 11.12.2015
 *      Author: david
 */

#ifndef COMBIGRID_SRC_SGPP_COMBIGRID_DEFINITIONS_HPP_
#define COMBIGRID_SRC_SGPP_COMBIGRID_DEFINITIONS_HPP_

#include <sgpp/globaldef.hpp>
#include <vector>
#include <cstddef>
#include <functional>

#define CGLOG(str)
// #include <iostream>
// #define CGLOG(str) std::cout << str << "\n"

namespace SGPP {
namespace combigrid {

typedef std::vector<size_t> MultiIndex;

template<typename In, typename Out> std::function<Out(In)> constantFunction(Out fixedValue = Out()) {
	return [=](In value){return fixedValue;};
}

template<typename Out> std::function<Out(MultiIndex const &)> multiIndexToDefaultValue(Out fixedValue = Out()) {
	return constantFunction<MultiIndex const &, Out>(fixedValue);
}

}
}
#endif /* COMBIGRID_SRC_SGPP_COMBIGRID_DEFINITIONS_HPP_ */
