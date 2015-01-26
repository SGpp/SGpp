/******************************************************************************
* Copyright (C) 2009 Technische Universitaet Muenchen                         *
* This file is part of the SG++ project. For conditions of distribution and   *
* use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
******************************************************************************/
// @author Dirk Pflueger (pflueged@in.tum.de)
// @author Benjamin (pehersto@in.tum.de)


#include <sgpp/base/basis/linear/boundary/operation/OperationQuadratureLinearBoundary.hpp>


namespace sg {
  namespace base {

    double OperationQuadratureLinearBoundary::doQuadrature(DataVector& alpha) {
      double res = 0;
      double tmp;
      int nr_boundaries = 0; //nr. of boundaries a point touches
      int cur_ind, cur_lev;
      GridStorage::index_type index;
      GridStorage::grid_map_iterator end_iter = storage->end();

      for (GridStorage::grid_map_iterator iter = storage->begin(); iter != end_iter; iter++) {
        tmp = pow(2.0, -static_cast<double>(iter->first->getLevelSum())) * alpha.get(iter->second);

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
