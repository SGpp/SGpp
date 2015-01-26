/* ****************************************************************************
* Copyright (C) 2011 Technische Universitaet Muenchen                         *
* This file is part of the SG++ project. For conditions of distribution and   *
* use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
**************************************************************************** */
// @author Stefanie Schraufstetter (schraufs@in.tum.de)

#include <sgpp/finance/tools/VariableDiscountFactor.hpp>

namespace sg {
  namespace finance {

    VariableDiscountFactor::VariableDiscountFactor(sg::base::GridStorage* storage, int dim_r): myBoundingBox(storage->getBoundingBox()), storage(storage), dim_r(dim_r) {
    }

    VariableDiscountFactor::~VariableDiscountFactor() {
    }

    void VariableDiscountFactor::getDiscountFactor(sg::base::DataVector& factor, double T) {
      double tmp;

      for (size_t i = 0; i < storage->size(); i++) {
        std::string coords = (*storage)[i]->getCoordsStringBB(*this->myBoundingBox);
        std::stringstream coordsStream(coords);
        double dblFuncValues[2];

        for (size_t j = 0; j < 2; j++) {
          coordsStream >> tmp;
          dblFuncValues[j] = tmp;
        }

        //std::cout<<dblFuncValues[1]<<std::endl;
        //factor.set(i, exp((-1.0)*dblFuncValues[1]*T));
        factor.set(i, exp((-1.0)*dblFuncValues[this->dim_r]*T));
      }
    }

  }
}
