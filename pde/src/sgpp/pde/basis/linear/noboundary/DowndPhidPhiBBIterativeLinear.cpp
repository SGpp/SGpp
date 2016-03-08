// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <sgpp/pde/basis/linear/noboundary/DowndPhidPhiBBIterativeLinear.hpp>
#include <sgpp/base/grid/common/BoundingBox.hpp>

#include <sgpp/globaldef.hpp>

namespace sgpp {
namespace pde {

DowndPhidPhiBBIterativeLinear::DowndPhidPhiBBIterativeLinear(sgpp::base::GridStorage* storage)
    : storage(storage) {}

DowndPhidPhiBBIterativeLinear::~DowndPhidPhiBBIterativeLinear() {}

void DowndPhidPhiBBIterativeLinear::operator()(sgpp::base::DataVector& alpha,
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
      // only affects the diagonal of the stiffness matrix
      result[i] = alpha[i] * (Qqout * (static_cast<double>(1 << (level + 1))));
    }
  } else {
    // traverse all basis function by sequence number
    for (size_t i = 0; i < storage->getSize(); i++) {
      sgpp::base::GridStorage::index_type::level_type level;
      sgpp::base::GridStorage::index_type::index_type index;
      (*storage)[i]->get(dim, level, index);
      // only affects the diagonal of the stiffness matrix
      result[i] = alpha[i] * static_cast<double>(1 << (level + 1));
    }
  }
}
}  // namespace pde
}  // namespace sgpp
