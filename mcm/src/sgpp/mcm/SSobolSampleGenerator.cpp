// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at 
// sgpp.sparsegrids.org

#include "SSobolSampleGenerator.hpp"

using namespace SGPP::base;

#include <sgpp/globaldef.hpp>


namespace SGPP {
  namespace mcm {
      int SSobolSampleGenerator::isOk(){
          return this->ok;
      } 

      void SSobolSampleGenerator::getSample(SGPP::base::DataVector& dv) {
          gen.next(dv.getPointer());
      }
    }
  }