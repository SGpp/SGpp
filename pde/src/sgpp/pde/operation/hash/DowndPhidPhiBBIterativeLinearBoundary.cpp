// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <sgpp/pde/operation/hash/DowndPhidPhiBBIterativeLinearBoundary.hpp>
#include <sgpp/base/grid/common/BoundingBox.hpp>

#include <sgpp/globaldef.hpp>

namespace sgpp {
namespace pde {

DowndPhidPhiBBIterativeLinearBoundary::DowndPhidPhiBBIterativeLinearBoundary(
    sgpp::base::GridStorage* storage)
    : storage(storage) {}

DowndPhidPhiBBIterativeLinearBoundary::~DowndPhidPhiBBIterativeLinearBoundary() {}

void DowndPhidPhiBBIterativeLinearBoundary::operator()(sgpp::base::DataVector& alpha,
                                                       sgpp::base::DataVector& result, size_t dim) {
  // Bounding Box handling
  sgpp::base::BoundingBox* boundingBox = this->storage->getBoundingBox();
  double q = boundingBox->getIntervalWidth(dim);
  double Qqout = 1.0 / q;

  // init the coefficients of the ansatz functions with boundary
  result.setAll(0.0);

  if (q != 1.0) {
    // traverse all basis function by sequence number
    for (size_t i = 0; i < storage->getSize(); i++) {
      sgpp::base::level_t level;
      sgpp::base::index_t index;
      (*storage)[i].get(dim, level, index);

      if (level == 0) {
        if (index == 0) {
          if (!boundingBox->hasDirichletBoundaryLeft(dim)) {
            // only affects the diagonal of the stiffness matrix
            result[i] += Qqout * alpha[i];

            // down
            if (index == 0) {
              sgpp::base::GridPoint index_one = (*storage)[i];
              index_one.set(dim, 0, 1);

              if (!boundingBox->hasDirichletBoundaryRight(dim)) {
                result[storage->getSequenceNumber(index_one)] += ((-1.0 * Qqout) * alpha[i]);
              }
            }
          }
        }

        if (index == 1) {
          if (!boundingBox->hasDirichletBoundaryRight(dim)) {
            // only affects the diagonal of the stiffness matrix
            result[i] += Qqout * alpha[i];
          }
        }
      } else {
        // only affects the diagonal of the stiffness matrix
        result[i] = alpha[i] * (Qqout * static_cast<double>(1 << (level + 1)));
      }
    }
  } else {
    // traverse all basis function by sequence number
    for (size_t i = 0; i < storage->getSize(); i++) {
      sgpp::base::level_t level;
      sgpp::base::index_t index;
      (*storage)[i].get(dim, level, index);

      if (level == 0) {
        if (index == 0) {
          if (!boundingBox->hasDirichletBoundaryLeft(dim)) {
            // only affects the diagonal of the stiffness matrix
            result[i] += alpha[i];

            // down
            if (index == 0) {
              sgpp::base::GridPoint index_one = (*storage)[i];
              index_one.set(dim, 0, 1);

              if (!boundingBox->hasDirichletBoundaryRight(dim)) {
                result[storage->getSequenceNumber(index_one)] += ((-1.0) * alpha[i]);
              }
            }
          }
        }

        if (index == 1) {
          if (!boundingBox->hasDirichletBoundaryRight(dim)) {
            // only affects the diagonal of the stiffness matrix
            result[i] += alpha[i];
          }
        }
      } else {
        // only affects the diagonal of the stiffness matrix
        result[i] = alpha[i] * static_cast<double>(1 << (level + 1));
      }
    }
  }
}
}  // namespace pde
}  // namespace sgpp
