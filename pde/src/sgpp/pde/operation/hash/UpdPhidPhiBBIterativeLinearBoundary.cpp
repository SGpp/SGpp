// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <sgpp/pde/operation/hash/UpdPhidPhiBBIterativeLinearBoundary.hpp>
#include <sgpp/base/grid/common/BoundingBox.hpp>

#include <sgpp/globaldef.hpp>

namespace sgpp {
namespace pde {

UpdPhidPhiBBIterativeLinearBoundary::UpdPhidPhiBBIterativeLinearBoundary(
    sgpp::base::GridStorage* storage)
    : storage(storage) {}

UpdPhidPhiBBIterativeLinearBoundary::~UpdPhidPhiBBIterativeLinearBoundary() {}

void UpdPhidPhiBBIterativeLinearBoundary::operator()(sgpp::base::DataVector& alpha,
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
      sgpp::base::GridStorage::index_type::level_type level;
      sgpp::base::GridStorage::index_type::index_type index;
      (*storage)[i]->get(dim, level, index);

      if (level == 0) {
        // up
        if (index == 1) {
          sgpp::base::GridIndex index_zero = *(*storage)[i];
          index_zero.set(dim, 0, 0);

          if (!boundingBox->hasDirichletBoundaryLeft(dim)) {
            result[storage->getSequenceNumber(&index_zero)] += ((-1.0 * Qqout) * alpha[i]);
          }
        }
      }
    }
  } else {
    // traverse all basis function by sequence number
    for (size_t i = 0; i < storage->getSize(); i++) {
      sgpp::base::GridStorage::index_type::level_type level;
      sgpp::base::GridStorage::index_type::index_type index;
      (*storage)[i]->get(dim, level, index);

      if (level == 0) {
        // up
        if (index == 1) {
          sgpp::base::GridIndex index_zero = *(*storage)[i];
          index_zero.set(dim, 0, 0);

          if (!boundingBox->hasDirichletBoundaryLeft(dim)) {
            result[storage->getSequenceNumber(&index_zero)] += ((-1.0) * alpha[i]);
          }
        }
      }
    }
  }
}
}  // namespace pde
}  // namespace sgpp
