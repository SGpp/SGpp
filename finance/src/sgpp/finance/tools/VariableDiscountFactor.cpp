// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <sgpp/finance/tools/VariableDiscountFactor.hpp>

#include <sgpp/globaldef.hpp>


namespace SGPP {
  namespace finance {

    VariableDiscountFactor::VariableDiscountFactor(SGPP::base::GridStorage* storage, int dim_r): myBoundingBox(storage->getBoundingBox()), storage(storage), dim_r(dim_r) {
    }

    VariableDiscountFactor::~VariableDiscountFactor() {
    }

    void VariableDiscountFactor::getDiscountFactor(SGPP::base::DataVector& factor, float_t T) {
      float_t tmp;

      for (size_t i = 0; i < storage->size(); i++) {
        std::string coords = (*storage)[i]->getCoordsStringBB(*this->myBoundingBox);
        std::stringstream coordsStream(coords);
        float_t dblFuncValues[2];

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