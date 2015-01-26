/* ****************************************************************************
* Copyright (C) 2013 Technische Universitaet Muenchen                         *
* This file is part of the SG++ project. For conditions of distribution and   *
* use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
**************************************************************************** */
// @author Roman Karlstetter (karlstet@in.tum.de)

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
        virtual double uniform(double a, double b);
        virtual double normal(double mean, double stddev);
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
