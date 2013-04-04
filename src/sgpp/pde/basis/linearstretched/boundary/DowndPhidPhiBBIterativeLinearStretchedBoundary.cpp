/******************************************************************************
* Copyright (C) 2011 Technische Universitaet Muenchen                         *
* This file is part of the SG++ project. For conditions of distribution and   *
* use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
******************************************************************************/
// @author Alexander Heinecke (Alexander.Heinecke@mytum.de, Sarpkan Selcuk (Sarpkan.Selcuk@mytum.de)

#include "pde/basis/linearstretched/boundary/DowndPhidPhiBBIterativeLinearStretchedBoundary.hpp"
#include "base/grid/common/Stretching.hpp"

namespace sg {
  namespace pde {

    DowndPhidPhiBBIterativeLinearStretchedBoundary::DowndPhidPhiBBIterativeLinearStretchedBoundary(sg::base::GridStorage* storage) : storage(storage) {
    }

    DowndPhidPhiBBIterativeLinearStretchedBoundary::~DowndPhidPhiBBIterativeLinearStretchedBoundary() {
    }

    void DowndPhidPhiBBIterativeLinearStretchedBoundary::operator()(sg::base::DataVector& alpha, sg::base::DataVector& result, size_t dim) {

      // Bounding Box handling
      sg::base::Stretching* stretching = this->storage->getStretching();
      double q = stretching->getIntervalWidth(dim);
      double Qqout = 1.0 / q;

      // init the coefficients of the ansatz functions with boundary
      result.setAll(0.0);

      if (q != 1.0) {
        // traverse all basis function by sequence number
        for (size_t i = 0; i < storage->size(); i++) {
          sg::base::GridStorage::index_type::level_type level;
          sg::base::GridStorage::index_type::index_type index;
          (*storage)[i]->get(dim, level, index);

          if (level == 0) {
            if (index == 0) {
              if (!stretching->hasDirichletBoundaryLeft(dim)) {
                //only affects the diagonal of the stiffness matrix
                result[i] += Qqout * alpha[i];

                // down
                if (index == 0) {
                  sg::base::GridIndex index_one = (*storage)[i];
                  index_one.set(dim, 0, 1);

                  if (!stretching->hasDirichletBoundaryRight(dim)) {
                    result[(*storage)[&index_one]] += ((-1.0 * Qqout) * alpha[i]);
                  }
                }
              }
            }

            if (index == 1) {
              if (!stretching->hasDirichletBoundaryRight(dim)) {
                //only affects the diagonal of the stiffness matrix
                result[i] += Qqout * alpha[i];
              }
            }
          }
          //only affects the diagonal of the stiffness matrix
          else {
            result[i] = alpha[i] * (Qqout * pow(2.0, static_cast<int>(level + 1)));
          }
        }
      } else {
        // traverse all basis function by sequence number
        for (size_t i = 0; i < storage->size(); i++) {
          sg::base::GridStorage::index_type::level_type level;
          sg::base::GridStorage::index_type::index_type index;
          (*storage)[i]->get(dim, level, index);

          if (level == 0) {
            if (index == 0) {
              if (!stretching->hasDirichletBoundaryLeft(dim)) {
                //only affects the diagonal of the stiffness matrix
                result[i] += alpha[i];

                // down
                if (index == 0) {
                  sg::base::GridIndex index_one = (*storage)[i];
                  index_one.set(dim, 0, 1);

                  if (!stretching->hasDirichletBoundaryRight(dim)) {
                    result[(*storage)[&index_one]] += ((-1.0) * alpha[i]);
                  }
                }
              }
            }

            if (index == 1) {
              if (!stretching->hasDirichletBoundaryRight(dim)) {
                //only affects the diagonal of the stiffness matrix
                result[i] += alpha[i];
              }
            }
          }
          //only affects the diagonal of the stiffness matrix
          else {
            result[i] = alpha[i] * pow(2.0, static_cast<int>(level + 1));
          }
        }
      }
    }

  }
}
