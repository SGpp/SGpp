/* ****************************************************************************
* Copyright (C) 2014 Universitaet Stuttgart                                   *
* This file is part of the SG++ project. For conditions of distribution and   *
* use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
**************************************************************************** */
// @author Andreas Doerr, Marcel Schneider, Matthias Moegerle

#include "OperationQuadratureMCAdvanced.hpp"

#include <sgpp/base/operation/BaseOpFactory.hpp>
#include <sgpp/base/datatypes/DataMatrix.hpp>
#include <sgpp/base/datatypes/DataVector.hpp>
#include <sgpp/mcm/SampleGenerator.hpp>
#include <sgpp/mcm/NaiveSampleGenerator.hpp>
#include <sgpp/mcm/SobolSampleGenerator.hpp>
#include <sgpp/mcm/ScrambledSobolSampleGenerator.hpp>
#include <sgpp/mcm/LatinHypercubeSampleGenerator.hpp>
#include <sgpp/mcm/StratifiedSampleGenerator.hpp>
#include <sgpp/mcm/SSobolSampleGenerator.hpp>

#include <cmath>
#include <iostream>
#include <sgpp/mcm/Random.hpp>

#include <sgpp/globaldef.hpp>


namespace SGPP {
  namespace mcm {


    OperationQuadratureMCAdvanced::OperationQuadratureMCAdvanced(SGPP::base::Grid& grid, int numberOfSamples) : grid(&grid), numberOfSamples(numberOfSamples) {
      this->dimensions = grid.getStorage()->dim();
      myGenerator = new SGPP::mcm::NaiveSampleGenerator(this->dimensions);
    }
    
    OperationQuadratureMCAdvanced::OperationQuadratureMCAdvanced(size_t dimensions, int numberOfSamples) : numberOfSamples(numberOfSamples) {
      this->dimensions = dimensions;
      myGenerator = new SGPP::mcm::NaiveSampleGenerator(this->dimensions);
      grid = NULL;
    }
    
    void OperationQuadratureMCAdvanced::useNaiveMonteCarlo(){
      myGenerator = new SGPP::mcm::NaiveSampleGenerator(dimensions);
    }

    void OperationQuadratureMCAdvanced::useQuasiMonteCarlo(){
      myGenerator = new SGPP::mcm::SobolSampleGenerator(dimensions, 0);
    }
    void OperationQuadratureMCAdvanced::useQuasiMonteCarloScrambled(){
      myGenerator = new SGPP::mcm::ScrambledSobolSampleGenerator(dimensions, 0);
    }

    void OperationQuadratureMCAdvanced::useStratifiedMonteCarlo(long long int* strataPerDimension){
      myGenerator = new SGPP::mcm::StratifiedSampleGenerator(dimensions, strataPerDimension);
    }

    void OperationQuadratureMCAdvanced::useLatinHypercubeMonteCarlo(){
      myGenerator = new SGPP::mcm::LatinHypercubeSampleGenerator(dimensions, numberOfSamples);
    }
    
    void OperationQuadratureMCAdvanced::useSSobol(int scrambling){
      myGenerator = new SGPP::mcm::SSobolSampleGenerator(dimensions, (int) numberOfSamples, scrambling);
    }
    
    double OperationQuadratureMCAdvanced::doQuadrature(SGPP::base::DataVector& alpha) {
      
      SGPP::base::DataMatrix dm(numberOfSamples, dimensions);
      
      myGenerator->getSamples(dm);
      
      SGPP::base::OperationMultipleEval* opEval = SGPP::op_factory::createOperationMultipleEval(*grid, &dm);
      SGPP::base::DataVector res = SGPP::base::DataVector(numberOfSamples);
      opEval->mult(alpha, res);
      return res.sum() / static_cast<double>(numberOfSamples);
    }

    double OperationQuadratureMCAdvanced::doQuadratureFunc(FUNC func, void* clientdata) {
      //double* p = new double[dimensions];

      SGPP::base::DataMatrix dm(numberOfSamples, dimensions);
      myGenerator->getSamples(dm);

      // create number of paths (uniformly drawn from [0,1]^d)
      double res = 0;

      for (size_t i = 0; i < numberOfSamples; i++) {
        SGPP::base::DataVector dv(dimensions);
        dm.getRow(i,dv);
        res += func(*reinterpret_cast<int*>(&dimensions), dv.getPointer(), clientdata);
      }

      //delete p;
      return res / static_cast<double>(numberOfSamples);
    }

    double OperationQuadratureMCAdvanced::doQuadratureL2Error(FUNC func, void* clientdata, SGPP::base::DataVector& alpha) {
      double x;
      double* p = new double[dimensions];

      SGPP::base::DataVector point(dimensions);
      SGPP::base::OperationEval* opEval = SGPP::op_factory::createOperationEval(*grid);
      // create number of paths (uniformly drawn from [0,1]^d)
      double res = 0;

      for (size_t i = 0; i < numberOfSamples; i++) {
        for (size_t d = 0; d < dimensions; d++) {
          x = Random::random_double();
          p[d] = x;
          point[d] = x;
        }

        res += pow(func(*reinterpret_cast<int*>(&dimensions), p, clientdata) - opEval->eval(alpha, point), 2);
      }

      delete p;
      return sqrt(res / static_cast<double>(numberOfSamples));
    }

    size_t OperationQuadratureMCAdvanced::getDimensions() {
      return dimensions;
    }

  }
  
  
}
