// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#ifndef COMBIGRID_SRC_SGPP_COMBIGRID_DEFINITIONS_HPP_
#define COMBIGRID_SRC_SGPP_COMBIGRID_DEFINITIONS_HPP_

#include <sgpp/globaldef.hpp>
#include <vector>
#include <cstddef>
#include <functional>

#define CGLOG(str)
// #include <iostream>
// #define CGLOG(str) std::cout << str << "\n"

namespace sgpp {
namespace combigrid {

typedef std::vector<size_t> MultiIndex;

template <typename In, typename Out>
std::function<Out(In)> constantFunction(Out fixedValue = Out()) {
  return [=](In value) { return fixedValue; };
}

template <typename Out>
std::function<Out(MultiIndex const &)> multiIndexToDefaultValue(Out fixedValue = Out()) {
  return constantFunction<MultiIndex const &, Out>(fixedValue);
}
}  // namespace combigrid
}  // namespace sgpp
#endif /* COMBIGRID_SRC_SGPP_COMBIGRID_DEFINITIONS_HPP_ */
