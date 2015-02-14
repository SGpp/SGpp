// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at 
// sgpp.sparsegrids.org

#include <sgpp/pde/operation/hash/DowndPhidPhiBBIterativeLinearBoundary.hpp>
#include <sgpp/base/grid/common/BoundingBox.hpp>

#include <sgpp/globaldef.hpp>


namespace SGPP {
  namespace pde {

    DowndPhidPhiBBIterativeLinearBoundary::DowndPhidPhiBBIterativeLinearBoundary(SGPP::base::GridStorage* storage) : storage(storage) {
    }

    DowndPhidPhiBBIterativeLinearBoundary::~DowndPhidPhiBBIterativeLinearBoundary() {
    }

    void DowndPhidPhiBBIterativeLinearBoundary::operator()(SGPP::base::DataVector& alpha, SGPP::base::DataVector& result, size_t dim) {
      // Bounding Box handling
      SGPP::base::BoundingBox* boundingBox = this->storage->getBoundingBox();
      float_t q = boundingBox->getIntervalWidth(dim);
      float_t Qqout = 1.0 / q;

      // init the coefficients of the ansatz functions with boundary
      result.setAll(0.0);

      if (q != 1.0) {
        // traverse all basis function by sequence number
        for (size_t i = 0; i < storage->size(); i++) {
          SGPP::base::GridStorage::index_type::level_type level;
          SGPP::base::GridStorage::index_type::index_type index;
          (*storage)[i]->get(dim, level, index);

          if (level == 0) {
            if (index == 0) {
              if (!boundingBox->hasDirichletBoundaryLeft(dim)) {
                //only affects the diagonal of the stiffness matrix
                result[i] += Qqout * alpha[i];

                // down
                if (index == 0) {
                  SGPP::base::GridIndex index_one = (*storage)[i];
                  index_one.set(dim, 0, 1);

                  if (!boundingBox->hasDirichletBoundaryRight(dim)) {
                    result[(*storage)[&index_one]] += ((-1.0 * Qqout) * alpha[i]);
                  }
                }
              }
            }

            if (index == 1) {
              if (!boundingBox->hasDirichletBoundaryRight(dim)) {
                //only affects the diagonal of the stiffness matrix
                result[i] += Qqout * alpha[i];
              }
            }
          }
          //only affects the diagonal of the stiffness matrix
          else {
            result[i] = alpha[i] * (Qqout * static_cast<float_t>(1 << (level + 1)));
          }
        }
      } else {
        // traverse all basis function by sequence number
        for (size_t i = 0; i < storage->size(); i++) {
          SGPP::base::GridStorage::index_type::level_type level;
          SGPP::base::GridStorage::index_type::index_type index;
          (*storage)[i]->get(dim, level, index);

          if (level == 0) {
            if (index == 0) {
              if (!boundingBox->hasDirichletBoundaryLeft(dim)) {
                //only affects the diagonal of the stiffness matrix
                result[i] += alpha[i];

                // down
                if (index == 0) {
                  SGPP::base::GridIndex index_one = (*storage)[i];
                  index_one.set(dim, 0, 1);

                  if (!boundingBox->hasDirichletBoundaryRight(dim)) {
                    result[(*storage)[&index_one]] += ((-1.0) * alpha[i]);
                  }
                }
              }
            }

            if (index == 1) {
              if (!boundingBox->hasDirichletBoundaryRight(dim)) {
                //only affects the diagonal of the stiffness matrix
                result[i] += alpha[i];
              }
            }
          }
          //only affects the diagonal of the stiffness matrix
          else {
            result[i] = alpha[i] * static_cast<float_t>(1 << (level + 1));
          }
        }
      }
    }

  }
}