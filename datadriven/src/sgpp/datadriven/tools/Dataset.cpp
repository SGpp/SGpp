/******************************************************************************
* Copyright (C) 2015 Technische Universitaet Muenchen                         *
* This file is part of the SG++ project. For conditions of distribution and   *
* use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
******************************************************************************/
// @author Julian Valentin (Julian.Valentin@ipvs.uni-stuttgart.de)

#include <sgpp/datadriven/tools/Dataset.hpp>

#include <sgpp/globaldef.hpp>

namespace sg {
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

    sg::base::DataVector* Dataset::getClasses() {
      return &classes;
    }

    sg::base::DataMatrix* Dataset::getTrainingData() {
      return &trainingData;
    }

  }
}
