/* ****************************************************************************
* Copyright (C) 2014 Universitaet Stuttgart                                   *
* This file is part of the SG++ project. For conditions of distribution and   *
* use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
**************************************************************************** */
// @author Andreas Doerr, Marcel Schneider, Matthias Moegerle

#include "OperationQuadratureMCAdvanced.hpp"

#include "base/operation/BaseOpFactory.hpp"
#include "base/datatypes/DataMatrix.hpp"
#include "base/datatypes/DataVector.hpp"
#include "mcm/SampleGenerator.hpp"
#include "mcm/NaiveSampleGenerator.hpp"
#include "mcm/SobolSampleGenerator.hpp"
#include "mcm/ScrambledSobolSampleGenerator.hpp"
#include "mcm/LatinHypercubeSampleGenerator.hpp"
#include "mcm/StratifiedSampleGenerator.hpp"
#include "mcm/SSobolSampleGenerator.hpp"

#include <cmath>
#include <iostream>
#include "mcm/Random.hpp"

namespace sg {
  namespace mcm {


    OperationQuadratureMCAdvanced::OperationQuadratureMCAdvanced(sg::base::Grid& grid, int numberOfSamples) : grid(&grid), numberOfSamples(numberOfSamples) {
      this->dimensions = grid.getStorage()->dim();
      myGenerator = new sg::mcm::NaiveSampleGenerator(this->dimensions);
    }
    
    OperationQuadratureMCAdvanced::OperationQuadratureMCAdvanced(size_t dimensions, int numberOfSamples) : numberOfSamples(numberOfSamples) {
      this->dimensions = dimensions;
      myGenerator = new sg::mcm::NaiveSampleGenerator(this->dimensions);
      grid = NULL;
    }
    
    void OperationQuadratureMCAdvanced::useNaiveMonteCarlo(){
      myGenerator = new sg::mcm::NaiveSampleGenerator(dimensions);
    }

    void OperationQuadratureMCAdvanced::useQuasiMonteCarlo(){
      myGenerator = new sg::mcm::SobolSampleGenerator(dimensions, 0);
    }
    void OperationQuadratureMCAdvanced::useQuasiMonteCarloScrambled(){
      myGenerator = new sg::mcm::ScrambledSobolSampleGenerator(dimensions, 0);
    }

    void OperationQuadratureMCAdvanced::useStratifiedMonteCarlo(long long int* strataPerDimension){
      myGenerator = new sg::mcm::StratifiedSampleGenerator(dimensions, strataPerDimension);
    }

    void OperationQuadratureMCAdvanced::useLatinHypercubeMonteCarlo(){
      myGenerator = new sg::mcm::LatinHypercubeSampleGenerator(dimensions, numberOfSamples);
    }
    
    void OperationQuadratureMCAdvanced::useSSobol(int scrambling){
      myGenerator = new sg::mcm::SSobolSampleGenerator(dimensions, (int) numberOfSamples, scrambling);
    }
    
    double OperationQuadratureMCAdvanced::doQuadrature(sg::base::DataVector& alpha) {
      
      sg::base::DataMatrix dm(numberOfSamples, dimensions);
      
      myGenerator->getSamples(dm);
      
      sg::base::OperationMultipleEval* opEval = sg::op_factory::createOperationMultipleEval(*grid, &dm);
      sg::base::DataVector res = sg::base::DataVector(numberOfSamples);
      opEval->mult(alpha, res);
      return res.sum() / static_cast<double>(numberOfSamples);
    }

    double OperationQuadratureMCAdvanced::doQuadratureFunc(FUNC func, void* clientdata) {
      //double* p = new double[dimensions];

      sg::base::DataMatrix dm(numberOfSamples, dimensions);
      myGenerator->getSamples(dm);

      // create number of paths (uniformly drawn from [0,1]^d)
      double res = 0;

      for (size_t i = 0; i < numberOfSamples; i++) {
        sg::base::DataVector dv(dimensions);
        dm.getRow(i,dv);
        res += func(*reinterpret_cast<int*>(&dimensions), dv.getPointer(), clientdata);
      }

      //delete p;
      return res / static_cast<double>(numberOfSamples);
    }

    double OperationQuadratureMCAdvanced::doQuadratureL2Error(FUNC func, void* clientdata, sg::base::DataVector& alpha) {
      double x;
      double* p = new double[dimensions];

      sg::base::DataVector point(dimensions);
      sg::base::OperationEval* opEval = sg::op_factory::createOperationEval(*grid);
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
