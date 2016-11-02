// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#ifndef COMBIGRID_SRC_SGPP_COMBIGRID_SERIALIZATION_DEFAULTSERIALIZATIONSTRATEGY_HPP_
#define COMBIGRID_SRC_SGPP_COMBIGRID_SERIALIZATION_DEFAULTSERIALIZATIONSTRATEGY_HPP_

#include <sgpp/combigrid/serialization/AbstractSerializationStrategy.hpp>

#include <sstream>
#include <string>
#include <vector>

namespace sgpp {
namespace combigrid {

/**
 * This serialization strategy uses the operators << and >> on streams that are overloaded for many
 * standard types. For floating-point types, FloatSerializationStrategy provides an exact
 * serialization method.
 */
template <typename T>
class DefaultSerializationStrategy : public AbstractSerializationStrategy<T> {
 public:
  virtual ~DefaultSerializationStrategy() {}

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
} /* namespace sgpp*/

#endif /* COMBIGRID_SRC_SGPP_COMBIGRID_SERIALIZATION_DEFAULTSERIALIZATIONSTRATEGY_HPP_ */
