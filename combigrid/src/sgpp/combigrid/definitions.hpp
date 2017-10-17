// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#ifndef COMBIGRID_SRC_SGPP_COMBIGRID_DEFINITIONS_HPP_
#define COMBIGRID_SRC_SGPP_COMBIGRID_DEFINITIONS_HPP_

#include <sgpp/globaldef.hpp>

#include <cstddef>

#include <functional>
#include <vector>

// defines for turning logging capabilities on or off
// #define CT_DEBUG

#ifndef CT_DEBUG
#define CGLOG(str)
#define CGLOG_SURROUND(cmd) cmd
#endif

#ifdef CT_DEBUG
#include <iostream>
#define CGLOG(str) std::cout << str << "\n"
#define CGLOG_SURROUND(cmd)                                                 \
  std::cout << #cmd << " before: " << __FILE__ << ", " << __LINE__ << "\n"; \
  cmd;                                                                      \
  std::cout << #cmd << " after: " << __FILE__ << ", " << __LINE__ << "\n"
#endif

namespace sgpp {
namespace combigrid {

typedef std::vector<size_t> MultiIndex;

/**
 * Returns a constant function with the given value. If no value is specified, the
 * default-constructed value is taken. The template-parameter In corresponding to the input
 * parameter type of the returned constant function has to be specified. If no fixed value is given,
 * the output type also has to be specified.
 */
template <typename In, typename Out>
std::function<Out(In)> constantFunction(Out fixedValue = Out()) {
  return [fixedValue](In value) { return fixedValue; };
}

/**
 * Returns a constant function with the given value and input type (MultiIndex const &). If no value
 * is specified, the default-constructed value is taken. If no fixed value is given, the output type
 * also has to be specified.
 */
template <typename Out>
std::function<Out(MultiIndex const &)> multiIndexToDefaultValue(Out fixedValue = Out()) {
  return constantFunction<MultiIndex const &, Out>(fixedValue);
}
}  // namespace combigrid
}  // namespace sgpp
#endif /* COMBIGRID_SRC_SGPP_COMBIGRID_DEFINITIONS_HPP_ */
