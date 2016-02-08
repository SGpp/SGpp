// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <sgpp/datadriven/tools/Dataset.hpp>

#include <sgpp/globaldef.hpp>

namespace SGPP {
namespace datadriven {

Dataset::Dataset() :
  numberInstances(0), dimension(0), classes(0), trainingData(0, 0) {
}

Dataset::Dataset(size_t numberInstances, size_t dimension) :
  numberInstances(numberInstances),
  dimension(dimension),
  classes(numberInstances),
  trainingData(numberInstances, dimension) {
}

size_t Dataset::getNumberInstances() const {
  return numberInstances;
}

size_t Dataset::getDimension() const {
  return dimension;
}

SGPP::base::DataVector& Dataset::getClasses() {
  return classes;
}

SGPP::base::DataMatrix& Dataset::getTrainingData() {
  return trainingData;
}

}  // namespace datadriven
}  // namespace SGPP

