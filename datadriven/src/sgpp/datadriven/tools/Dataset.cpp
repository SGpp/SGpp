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

sgpp::base::DataVector& Dataset::getTargets() { return targets; }

sgpp::base::DataMatrix& Dataset::getData() { return data; }

void Dataset::setData(const sgpp::base::DataMatrix& data) { this->data = data; }

void Dataset::setTargets(const sgpp::base::DataVector& targets) { this->targets = targets; }

}  // namespace datadriven
}  // namespace sgpp
