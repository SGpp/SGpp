// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <sgpp/pde/basis/linearstretched/boundary/DowndPhidPhiBBIterativeLinearStretchedBoundary.hpp>
#include <sgpp/base/grid/common/Stretching.hpp>

#include <sgpp/globaldef.hpp>

namespace sgpp {
namespace pde {

DowndPhidPhiBBIterativeLinearStretchedBoundary::DowndPhidPhiBBIterativeLinearStretchedBoundary(
    sgpp::base::GridStorage* storage)
    : storage(storage) {}

DowndPhidPhiBBIterativeLinearStretchedBoundary::~DowndPhidPhiBBIterativeLinearStretchedBoundary() {}

void DowndPhidPhiBBIterativeLinearStretchedBoundary::operator()(sgpp::base::DataVector& alpha,
                                                                sgpp::base::DataVector& result,
                                                                size_t dim) {
  // Bounding Box handling
  sgpp::base::Stretching* stretching = this->storage->getStretching();
  double q = stretching->getIntervalWidth(dim);
  double Qqout = 1.0 / q;

  // init the coefficients of the ansatz functions with boundary
  result.setAll(0.0);

  if (q != 1.0) {
    // traverse all basis function by sequence number
    for (size_t i = 0; i < storage->getSize(); i++) {
      sgpp::base::GridStorage::index_type::level_type level;
      sgpp::base::GridStorage::index_type::index_type index;
      (*storage)[i]->get(dim, level, index);

      if (level == 0) {
        if (index == 0) {
          if (!stretching->hasDirichletBoundaryLeft(dim)) {
            // only affects the diagonal of the stiffness matrix
            result[i] += Qqout * alpha[i];

            // down
            if (index == 0) {
              sgpp::base::GridIndex index_one = *(*storage)[i];
              index_one.set(dim, 0, 1);

              if (!stretching->hasDirichletBoundaryRight(dim)) {
                result[storage->getSequenceNumber(&index_one)] += ((-1.0 * Qqout) * alpha[i]);
              }
            }
          }
        }

        if (index == 1) {
          if (!stretching->hasDirichletBoundaryRight(dim)) {
            // only affects the diagonal of the stiffness matrix
            result[i] += Qqout * alpha[i];
          }
        }
      } else {
        // only affects the diagonal of the stiffness matrix
        result[i] = alpha[i] * (Qqout * pow(2.0, static_cast<int>(level + 1)));
      }
    }
  } else {
    // traverse all basis function by sequence number
    for (size_t i = 0; i < storage->getSize(); i++) {
      sgpp::base::GridStorage::index_type::level_type level;
      sgpp::base::GridStorage::index_type::index_type index;
      (*storage)[i]->get(dim, level, index);

      if (level == 0) {
        if (index == 0) {
          if (!stretching->hasDirichletBoundaryLeft(dim)) {
            // only affects the diagonal of the stiffness matrix
            result[i] += alpha[i];

            // down
            if (index == 0) {
              sgpp::base::GridIndex index_one = *(*storage)[i];
              index_one.set(dim, 0, 1);

              if (!stretching->hasDirichletBoundaryRight(dim)) {
                result[storage->getSequenceNumber(&index_one)] += ((-1.0) * alpha[i]);
              }
            }
          }
        }

        if (index == 1) {
          if (!stretching->hasDirichletBoundaryRight(dim)) {
            // only affects the diagonal of the stiffness matrix
            result[i] += alpha[i];
          }
        }
      } else {
        // only affects the diagonal of the stiffness matrix
        result[i] = alpha[i] * pow(2.0, static_cast<int>(level + 1));
      }
    }
  }
}
}  // namespace pde
}  // namespace sgpp
