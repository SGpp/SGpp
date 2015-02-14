// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at 
// sgpp.sparsegrids.org

#include <sgpp/base/operation/hash/OperationQuadratureLinearBoundary.hpp>


#include <sgpp/globaldef.hpp>


namespace SGPP {
  namespace base {

    float_t OperationQuadratureLinearBoundary::doQuadrature(DataVector& alpha) {
      float_t res = 0;
      float_t tmp;
      int nr_boundaries = 0; //nr. of boundaries a point touches
      int cur_ind, cur_lev;
      GridStorage::index_type index;
      GridStorage::grid_map_iterator end_iter = storage->end();

      for (GridStorage::grid_map_iterator iter = storage->begin(); iter != end_iter; iter++) {
        tmp = pow(2.0, -static_cast<float_t>(iter->first->getLevelSum())) * alpha.get(iter->second);

        if (!iter->first->isInnerPoint()) {
          nr_boundaries = 0;

          for (size_t d = 0; d < iter->first->dim(); d++) {
            cur_ind = iter->first->getIndex(d);
            cur_lev = iter->first->getLevel(d);

            if ((cur_ind == 0) || (pow(2.0, cur_lev) == cur_ind))
              nr_boundaries++;
          }

          tmp *= (pow(2.0, -nr_boundaries));
        }

        res += tmp;
      }

      return res;
    }

  }
}