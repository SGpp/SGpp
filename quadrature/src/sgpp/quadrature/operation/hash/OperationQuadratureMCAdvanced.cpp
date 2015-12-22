// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <sgpp/quadrature/operation/hash/OperationQuadratureMCAdvanced.hpp>
#include <sgpp/base/operation/BaseOpFactory.hpp>
#include <sgpp/base/datatypes/DataMatrix.hpp>
#include <sgpp/base/datatypes/DataVector.hpp>
#include <cmath>
#include <iostream>
#include <sgpp/globaldef.hpp>
#include <sgpp/quadrature/Random.hpp>
#include <sgpp/quadrature/sampling/HaltonSampleGenerator.hpp>
#include <sgpp/quadrature/sampling/LatinHypercubeSampleGenerator.hpp>
#include <sgpp/quadrature/sampling/NaiveSampleGenerator.hpp>
#include <sgpp/quadrature/sampling/StratifiedSampleGenerator.hpp>

namespace SGPP {
  namespace quadrature {

    OperationQuadratureMCAdvanced::OperationQuadratureMCAdvanced(
      SGPP::base::Grid& grid, size_t numberOfSamples, std::uint64_t seed) :
      grid(&grid), numberOfSamples(numberOfSamples), seed(seed) {
      dimensions = grid.getStorage()->dim();
      myGenerator = new SGPP::quadrature::NaiveSampleGenerator(dimensions, seed);
    }

    OperationQuadratureMCAdvanced::OperationQuadratureMCAdvanced(size_t dimensions,
        size_t numberOfSamples, std::uint64_t seed) :
      grid(NULL), numberOfSamples(numberOfSamples), dimensions(dimensions), seed(seed) {
      myGenerator = new SGPP::quadrature::NaiveSampleGenerator(dimensions, seed);
    }

    OperationQuadratureMCAdvanced::~OperationQuadratureMCAdvanced() {
      if (myGenerator != NULL) {
        delete myGenerator;
      }
    }

    void OperationQuadratureMCAdvanced::useNaiveMonteCarlo() {
      if (myGenerator != NULL) {
        delete myGenerator;
      }
      myGenerator = new SGPP::quadrature::NaiveSampleGenerator(dimensions, seed);
    }

    void OperationQuadratureMCAdvanced::useStratifiedMonteCarlo(
      std::vector<size_t>& strataPerDimension) {
      if (myGenerator != NULL) {
        delete myGenerator;
      }
      myGenerator = new SGPP::quadrature::StratifiedSampleGenerator(
        strataPerDimension, seed);
    }

    void OperationQuadratureMCAdvanced::useLatinHypercubeMonteCarlo() {
      if (myGenerator != NULL) {
        delete myGenerator;
      }
      myGenerator = new SGPP::quadrature::LatinHypercubeSampleGenerator(
        dimensions, numberOfSamples, seed);
    }

    void OperationQuadratureMCAdvanced::useQuasiMonteCarloWithHaltonSequences() {
      if (myGenerator != NULL) {
        delete myGenerator;
      }
      myGenerator = new SGPP::quadrature::HaltonSampleGenerator(dimensions);
    }

    float_t OperationQuadratureMCAdvanced::doQuadrature(
      SGPP::base::DataVector& alpha) {

      SGPP::base::DataMatrix dm(numberOfSamples, dimensions);

      myGenerator->getSamples(dm);

      SGPP::base::OperationMultipleEval* opEval =
        SGPP::op_factory::createOperationMultipleEval(*grid, dm);
      SGPP::base::DataVector res = SGPP::base::DataVector(numberOfSamples);
      opEval->mult(alpha, res);
      return res.sum() / static_cast<float_t>(numberOfSamples);
    }

    float_t OperationQuadratureMCAdvanced::doQuadratureFunc(FUNC func,
        void* clientdata) {
      //float_t* p = new float_t[dimensions];

      SGPP::base::DataMatrix dm(numberOfSamples, dimensions);
      myGenerator->getSamples(dm);

      // create number of paths (uniformly drawn from [0,1]^d)
      float_t res = 0;

      for (size_t i = 0; i < numberOfSamples; i++) {
        SGPP::base::DataVector dv(dimensions);
        dm.getRow(i, dv);
        res += func(*reinterpret_cast<int*>(&dimensions), dv.getPointer(),
                    clientdata);
      }

      //delete p;
      return res / static_cast<float_t>(numberOfSamples);
    }

    float_t OperationQuadratureMCAdvanced::doQuadratureL2Error(FUNC func,
        void* clientdata, SGPP::base::DataVector& alpha) {
      float_t x;
      float_t* p = new float_t[dimensions];

      SGPP::base::DataVector point(dimensions);
      SGPP::base::OperationEval* opEval = SGPP::op_factory::createOperationEval(
                                            *grid);
      // create number of paths (uniformly drawn from [0,1]^d)
      float_t res = 0;

      for (size_t i = 0; i < numberOfSamples; i++) {
        for (size_t d = 0; d < dimensions; d++) {
          x = Random::random_double();
          p[d] = x;
          point[d] = x;
        }

        res += pow(
                 func(*reinterpret_cast<int*>(&dimensions), p, clientdata)
                 - opEval->eval(alpha, point), 2);
      }

      delete p;
      return sqrt(res / static_cast<float_t>(numberOfSamples));
    }

    size_t OperationQuadratureMCAdvanced::getDimensions() {
      return dimensions;
    }

  }

}
