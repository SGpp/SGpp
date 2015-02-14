// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at 
// sgpp.sparsegrids.org

#ifndef DATASETGENERATOR_HPP
#define DATASETGENERATOR_HPP

#include <sgpp/base/datatypes/DataVector.hpp>
#include <sgpp/base/datatypes/DataMatrix.hpp>

#include <sgpp/globaldef.hpp>


namespace SGPP {
  namespace datadriven {
    class DatasetGenerator {
      public:
        virtual ~DatasetGenerator();
        virtual float_t uniform(float_t a, float_t b);
        virtual float_t normal(float_t mean, float_t stddev);
        virtual void createData(size_t offset, size_t size, SGPP::base::DataMatrix& trainingData, SGPP::base::DataVector& classes) = 0;
        virtual size_t getDims() = 0;
    };

    class Friedman1Generator : public DatasetGenerator {
      public:
        virtual void createData(size_t offset, size_t size, SGPP::base::DataMatrix& trData, SGPP::base::DataVector& classes);
        virtual size_t getDims();
    };

    class Friedman2Generator : public DatasetGenerator {
      public:
        virtual void createData(size_t offset, size_t size, SGPP::base::DataMatrix& trData, SGPP::base::DataVector& classes);
        virtual size_t getDims();
    };

    class Friedman3Generator : public DatasetGenerator {
      public:
        virtual void createData(size_t offset, size_t size, SGPP::base::DataMatrix& trData, SGPP::base::DataVector& classes);
        virtual size_t getDims();
    };

  }
}

#endif // DATASETGENERATOR_HPP