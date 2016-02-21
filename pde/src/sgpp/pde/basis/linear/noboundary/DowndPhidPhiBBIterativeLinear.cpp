// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <sgpp/pde/basis/linear/noboundary/DowndPhidPhiBBIterativeLinear.hpp>
#include <sgpp/base/grid/common/BoundingBox.hpp>

#include <sgpp/globaldef.hpp>

namespace SGPP {
namespace pde {

DowndPhidPhiBBIterativeLinear::DowndPhidPhiBBIterativeLinear(SGPP::base::GridStorage* storage)
    : storage(storage) {}

DowndPhidPhiBBIterativeLinear::~DowndPhidPhiBBIterativeLinear() {}

void DowndPhidPhiBBIterativeLinear::operator()(SGPP::base::DataVector& alpha,
                                               SGPP::base::DataVector& result, size_t dim) {
  // Bounding Box handling
  SGPP::base::BoundingBox* boundingBox = this->storage->getBoundingBox();
  float_t q = boundingBox->getIntervalWidth(dim);
  float_t Qqout = 1.0 / q;

  // init the coefficients of the ansatz functions with boundary
  result.setAll(0.0);

  if (q != 1.0) {
    // traverse all basis function by sequence number
    for (size_t i = 0; i < storage->getSize(); i++) {
      SGPP::base::GridStorage::index_type::level_type level;
      SGPP::base::GridStorage::index_type::index_type index;
      (*storage)[i]->get(dim, level, index);
      // only affects the diagonal of the stiffness matrix
      result[i] = alpha[i] * (Qqout * (static_cast<float_t>(1 << (level + 1))));
    }
  } else {
    // traverse all basis function by sequence number
    for (size_t i = 0; i < storage->getSize(); i++) {
      SGPP::base::GridStorage::index_type::level_type level;
      SGPP::base::GridStorage::index_type::index_type index;
      (*storage)[i]->get(dim, level, index);
      // only affects the diagonal of the stiffness matrix
      result[i] = alpha[i] * static_cast<float_t>(1 << (level + 1));
    }
  }
}
}  // namespace pde
}  // namespace SGPP
