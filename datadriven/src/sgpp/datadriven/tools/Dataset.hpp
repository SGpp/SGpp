// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#ifndef DATASET_HPP
#define DATASET_HPP

#include <sgpp/base/datatypes/DataVector.hpp>
#include <sgpp/base/datatypes/DataMatrix.hpp>

#include <cstddef>

#include <sgpp/globaldef.hpp>

namespace SGPP {
  namespace datadriven {

    class Dataset {
      public:
        /**
         * Constructs an empty dataset (zero size).
         */
        Dataset();

        /**
         * Constructs an empty dataset with given size.
         *
         * @param numberInstances number of instances in the dataset
         * @param dimension number of dimensions in the dataset
         */
        Dataset(size_t numberInstances, size_t dimension);

        /**
         * @return number of instances in the dataset
         */
        size_t getNumberInstances() const;

        /**
         * @return number of dimensions in the dataset
         */
        size_t getDimension() const;

        /**
         * @return classes data of the dataset
         */
        SGPP::base::DataVector& getClasses();

        /**
         * @return training data of the dataset
         */
        SGPP::base::DataMatrix& getTrainingData();

      protected:
        size_t numberInstances;
        size_t dimension;
        SGPP::base::DataVector classes;
        SGPP::base::DataMatrix trainingData;
    };

  }
}

#endif
