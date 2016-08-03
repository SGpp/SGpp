// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#ifndef COMBIGRID_SRC_SGPP_COMBIGRID_SERIALIZATION_ABSTRACTSERIALIZATIONSTRATEGY_HPP_
#define COMBIGRID_SRC_SGPP_COMBIGRID_SERIALIZATION_ABSTRACTSERIALIZATIONSTRATEGY_HPP_

#include <sgpp/globaldef.hpp>

#include <string>

namespace sgpp {
namespace combigrid {

template <typename T>
class AbstractSerializationStrategy {
 public:
  virtual ~AbstractSerializationStrategy() {}

  virtual std::string serialize(T const &value) = 0;
  virtual T deserialize(std::string const &input) = 0;
};

} /* namespace combigrid */
} /* namespace sgpp*/

#endif /* COMBIGRID_SRC_SGPP_COMBIGRID_SERIALIZATION_ABSTRACTSERIALIZATIONSTRATEGY_HPP_ */
