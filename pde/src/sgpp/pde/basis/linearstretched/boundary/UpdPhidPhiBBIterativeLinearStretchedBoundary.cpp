// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <sgpp/pde/basis/linearstretched/boundary/UpdPhidPhiBBIterativeLinearStretchedBoundary.hpp>
#include <sgpp/base/grid/common/Stretching.hpp>

#include <sgpp/globaldef.hpp>

namespace SGPP {
namespace pde {

UpdPhidPhiBBIterativeLinearStretchedBoundary::UpdPhidPhiBBIterativeLinearStretchedBoundary(
    SGPP::base::GridStorage* storage)
    : storage(storage) {}

UpdPhidPhiBBIterativeLinearStretchedBoundary::~UpdPhidPhiBBIterativeLinearStretchedBoundary() {}

void UpdPhidPhiBBIterativeLinearStretchedBoundary::operator()(SGPP::base::DataVector& alpha,
                                                              SGPP::base::DataVector& result,
                                                              size_t dim) {
  // Bounding Box handling
  SGPP::base::Stretching* stretching = this->storage->getStretching();
  float_t q = stretching->getIntervalWidth(dim);
  float_t Qqout = 1.0 / q;

  // init the coefficients of the ansatz functions with boundary
  result.setAll(0.0);

  if (q != 1.0) {
    // traverse all basis function by sequence number
    for (size_t i = 0; i < storage->getSize(); i++) {
      SGPP::base::GridStorage::index_type::level_type level;
      SGPP::base::GridStorage::index_type::index_type index;
      (*storage)[i]->get(dim, level, index);

      if (level == 0) {
        // up
        if (index == 1) {
          SGPP::base::GridIndex index_zero = *(*storage)[i];
          index_zero.set(dim, 0, 0);

          if (!stretching->hasDirichletBoundaryLeft(dim)) {
            result[storage->seq(&index_zero)] += ((-1.0 * Qqout) * alpha[i]);
          }
        }
      }
    }
  } else {
    // traverse all basis function by sequence number
    for (size_t i = 0; i < storage->getSize(); i++) {
      SGPP::base::GridStorage::index_type::level_type level;
      SGPP::base::GridStorage::index_type::index_type index;
      (*storage)[i]->get(dim, level, index);

      if (level == 0) {
        // up
        if (index == 1) {
          SGPP::base::GridIndex index_zero = *(*storage)[i];
          index_zero.set(dim, 0, 0);

          if (!stretching->hasDirichletBoundaryLeft(dim)) {
            result[storage->seq(&index_zero)] += ((-1.0) * alpha[i]);
          }
        }
      }
    }
  }
}
}  // namespace pde
}  // namespace SGPP
