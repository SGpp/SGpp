// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <sgpp/pde/operation/hash/UpdPhidPhiBBIterativeLinearBoundary.hpp>
#include <sgpp/base/grid/common/BoundingBox.hpp>

#include <sgpp/globaldef.hpp>


namespace SGPP {
  namespace pde {

    UpdPhidPhiBBIterativeLinearBoundary::UpdPhidPhiBBIterativeLinearBoundary(SGPP::base::GridStorage* storage) : storage(storage) {
    }

    UpdPhidPhiBBIterativeLinearBoundary::~UpdPhidPhiBBIterativeLinearBoundary() {
    }

    void UpdPhidPhiBBIterativeLinearBoundary::operator()(SGPP::base::DataVector& alpha, SGPP::base::DataVector& result, size_t dim) {
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
            // up
            if (index == 1) {
              SGPP::base::GridIndex index_zero = *(*storage)[i];
              index_zero.set(dim, 0, 0);

              if (!boundingBox->hasDirichletBoundaryLeft(dim)) {
                result[(*storage)[&index_zero]] += ((-1.0 * Qqout) * alpha[i]);
              }
            }
          }
        }
      } else {
        // traverse all basis function by sequence number
        for (size_t i = 0; i < storage->size(); i++) {
          SGPP::base::GridStorage::index_type::level_type level;
          SGPP::base::GridStorage::index_type::index_type index;
          (*storage)[i]->get(dim, level, index);

          if (level == 0) {
            // up
            if (index == 1) {
              SGPP::base::GridIndex index_zero = *(*storage)[i];
              index_zero.set(dim, 0, 0);

              if (!boundingBox->hasDirichletBoundaryLeft(dim)) {
                result[(*storage)[&index_zero]] += ((-1.0) * alpha[i]);
              }
            }
          }
        }
      }
    }

  }
}
