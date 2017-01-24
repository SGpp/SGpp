// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <sgpp/datadriven/tools/Dataset.hpp>

#include <sgpp/globaldef.hpp>

namespace sgpp {
namespace datadriven {

Dataset::Dataset() : numberInstances(0), dimension(0), targets(0), data(0, 0) {}

Dataset::Dataset(size_t numberInstances, size_t dimension)
    : numberInstances(numberInstances),
      dimension(dimension),
      targets(numberInstances),
      data(numberInstances, dimension) {}

size_t Dataset::getNumberInstances() const { return numberInstances; }

size_t Dataset::getDimension() const { return dimension; }

sgpp::base::DataVector& Dataset::getTargets() {
  return const_cast<sgpp::base::DataVector&>(static_cast<Dataset&>(*this).getTargets());
}

sgpp::base::DataMatrix& Dataset::getData() {
  return const_cast<sgpp::base::DataMatrix&>(static_cast<Dataset&>(*this).getData());
}

const sgpp::base::DataVector& Dataset::getTargets() const { return targets; }

const sgpp::base::DataMatrix& Dataset::getData() const { return data; }

}  // namespace datadriven
}  // namespace sgpp
