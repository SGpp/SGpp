/******************************************************************************
* Copyright (C) 2009 Technische Universitaet Muenchen                         *
* This file is part of the SG++ project. For conditions of distribution and   *
* use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
******************************************************************************/
// @author Dirk Pflueger (pflueged@in.tum.de)
// @author Benjamin


#include "base/basis/linear/noboundary/operation/OperationFirstMomentLinear.hpp"


namespace sg {
  namespace base {

    double OperationFirstMomentLinear::doQuadrature(DataVector& alpha) {
      double res = 0;
      double tmpres = 1;
      GridStorage::index_type index;
      GridStorage::grid_map_iterator end_iter = storage->end();

      for (GridStorage::grid_map_iterator iter = storage->begin(); iter != end_iter; iter++) {
        //    index = *(iter->first);
        //    std::cout << iter->second << " " << iter->first->getLevelSum() << " " << pow(2.0, -static_cast<double>(iter->first->getLevelSum())) << std::endl;
        tmpres = 1.;

        for (size_t dim = 0; dim < storage->dim(); dim++)
          tmpres *= iter->first->getIndex(dim) * pow(4.0, -static_cast<double>(iter->first->getLevel(dim)));

        res += alpha.get(iter->second) * tmpres;
      }

      return res;
    }

  }
}
