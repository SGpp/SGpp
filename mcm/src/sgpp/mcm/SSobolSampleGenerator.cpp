/* ****************************************************************************
* Copyright (C) 2014 Universitaet Stuttgart                                   *
* This file is part of the SG++ project. For conditions of distribution and   *
* use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
**************************************************************************** */
// @author Andreas Doerr, Marcel Schneider, Matthias Moegerle

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
