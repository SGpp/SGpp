// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <sgpp/pde/algorithm/UpDownOneOpDimWithShadow.hpp>
#include <sgpp/base/tools/SGppStopwatch.hpp>

#include <sgpp/globaldef.hpp>

namespace sgpp {
namespace pde {

UpDownOneOpDimWithShadow::UpDownOneOpDimWithShadow(sgpp::base::GridStorage* storage,
                                                   sgpp::base::GridStorage* shadowStorage) {
  this->storage = storage;
  this->shadowStorage = shadowStorage;
}

UpDownOneOpDimWithShadow::~UpDownOneOpDimWithShadow() {}

void UpDownOneOpDimWithShadow::mult(sgpp::base::DataVector& alpha, sgpp::base::DataVector& result) {
  expandGrid();

  // Create new Datavectors for the grid including shadow points.
  sgpp::base::DataVector alpha_temp(storage->getSize());
  sgpp::base::DataVector result_temp(storage->getSize());

  alpha_temp.setAll(0.0);
  result_temp.setAll(0.0);

  // Copy the alpha vector ... the old points remain the same, the new points are zero.
  for (size_t i = 0; i < alpha.getSize(); i++) {
    alpha_temp[i] = alpha[i];
  }

  sgpp::base::DataVector beta(result_temp.getSize());
  beta.setAll(0.0);

  for (size_t i = 0; i < storage->getDimension(); i++) {
    this->updown(alpha_temp, beta, storage->getDimension() - 1, i);

    result_temp.add(beta);
  }

  // Remove shadow grid points from the grid
  shrinkGrid();

  // Copy the result in the actual result vector. The values for the shadow
  // grid points are just needed for data transport, thus they can be omitted
  for (size_t i = 0; i < storage->getSize(); i++) {
    result[i] = result_temp[i];
  }
}

void UpDownOneOpDimWithShadow::expandGrid() {
  for (size_t i = 0; i < shadowStorage->getSize(); i++) {
    storage->insert(*shadowStorage->getGridIndex(i));
  }
}

void UpDownOneOpDimWithShadow::shrinkGrid() {
  for (size_t i = 0; i < shadowStorage->getSize(); i++) {
    storage->deleteLast();
  }
}

void UpDownOneOpDimWithShadow::updown(sgpp::base::DataVector& alpha, sgpp::base::DataVector& result,
                                      size_t dim, size_t op_dim) {
  if (dim == op_dim) {
    specialOP(alpha, result, dim, op_dim);
  } else {
    // Unidirectional scheme
    if (dim > 0) {
      // Reordering ups and downs
      sgpp::base::DataVector temp(alpha.getSize());
      temp.setAll(0.0);
      up(alpha, temp, dim);
      updown(temp, result, dim - 1, op_dim);

      // Same from the other direction:
      sgpp::base::DataVector result_temp(alpha.getSize());
      result_temp.setAll(0.0);
      updown(alpha, temp, dim - 1, op_dim);
      down(temp, result_temp, dim);

      result.add(result_temp);
    } else {
      // Terminates dimension recursion
      up(alpha, result, dim);

      sgpp::base::DataVector temp(alpha.getSize());
      temp.setAll(0.0);
      down(alpha, temp, dim);

      result.add(temp);
    }
  }
}

void UpDownOneOpDimWithShadow::specialOP(sgpp::base::DataVector& alpha,
                                         sgpp::base::DataVector& result, size_t dim,
                                         size_t op_dim) {
  // Unidirectional scheme
  if (dim > 0) {
    // Reordering ups and downs
    sgpp::base::DataVector temp(alpha.getSize());
    temp.setAll(0.0);
    upOpDim(alpha, temp, dim);
    updown(temp, result, dim - 1, op_dim);

    // Same from the other direction:
    sgpp::base::DataVector result_temp(alpha.getSize());
    result_temp.setAll(0.0);
    updown(alpha, temp, dim - 1, op_dim);
    downOpDim(temp, result_temp, dim);

    result.add(result_temp);
  } else {
    // Terminates dimension recursion
    upOpDim(alpha, result, dim);

    sgpp::base::DataVector temp(alpha.getSize());
    temp.setAll(0.0);
    downOpDim(alpha, temp, dim);

    result.add(temp);
  }
}
}  // namespace pde
}  // namespace sgpp
