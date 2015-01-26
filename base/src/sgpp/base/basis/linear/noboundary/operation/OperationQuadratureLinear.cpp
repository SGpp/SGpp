/******************************************************************************
* Copyright (C) 2009 Technische Universitaet Muenchen                         *
* This file is part of the SG++ project. For conditions of distribution and   *
* use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
******************************************************************************/
// @author Dirk Pflueger (pflueged@in.tum.de)


#include <sgpp/base/basis/linear/noboundary/operation/OperationQuadratureLinear.hpp>


namespace sg {
  namespace base {

    double OperationQuadratureLinear::doQuadrature(DataVector& alpha) {
      double res = 0;
      GridStorage::index_type index;
      GridStorage::grid_map_iterator end_iter = storage->end();

      for (GridStorage::grid_map_iterator iter = storage->begin(); iter != end_iter; iter++) {
        //    index = *(iter->first);
        //    std::cout << iter->second << " " << iter->first->getLevelSum() << " " << pow(2.0, -static_cast<double>(iter->first->getLevelSum())) << std::endl;
        res += pow(2.0, -static_cast<double>(iter->first->getLevelSum())) * alpha.get(iter->second);
      }

      return res;
    }

  }
}
