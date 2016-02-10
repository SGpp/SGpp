// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <sgpp/pde/algorithm/UpDownOneOpDimWithShadow.hpp>
#include <sgpp/base/tools/SGppStopwatch.hpp>

#include <sgpp/globaldef.hpp>


namespace SGPP {
namespace pde {

UpDownOneOpDimWithShadow::UpDownOneOpDimWithShadow(SGPP::base::GridStorage*
    storage,
    SGPP::base::GridStorage* shadowStorage) {
  this->storage = storage;
  this->shadowStorage = shadowStorage;
}

UpDownOneOpDimWithShadow::~UpDownOneOpDimWithShadow() {
}

void UpDownOneOpDimWithShadow::mult(SGPP::base::DataVector& alpha,
                                    SGPP::base::DataVector& result) {

  expandGrid();

  //Create new Datavectors for the grid including shadow points.
  SGPP::base::DataVector alpha_temp(storage->size());
  SGPP::base::DataVector result_temp(storage->size());

  alpha_temp.setAll(0.0);
  result_temp.setAll(0.0);

  //Copy the alpha vector ... the old points remain the same, the new points are zero.
  for (size_t i = 0; i < alpha.getSize(); i++) {
    alpha_temp[i] = alpha[i];
  }

  SGPP::base::DataVector beta(result_temp.getSize());
  beta.setAll(0.0);

  for (size_t i = 0; i < storage->dim(); i++) {
    this->updown(alpha_temp, beta, storage->dim() - 1, i);

    result_temp.add(beta);

  }

  //Remove shadow grid points from the grid
  shrinkGrid();

  //Copy the result in the actual result vector. The values for the shadow
  //grid points are just needed for data transport, thus they can be omitted
  for (size_t i = 0; i < storage->size(); i++) {
    result[i] = result_temp[i];
  }
}

void UpDownOneOpDimWithShadow::expandGrid() {
  for (size_t i = 0; i < shadowStorage->size(); i++) {
    storage->insert(*shadowStorage->get(i));
  }
}

void UpDownOneOpDimWithShadow::shrinkGrid() {
  for (size_t i = 0; i < shadowStorage->size(); i++) {
    storage->deleteLast();
  }
}

void UpDownOneOpDimWithShadow::updown(SGPP::base::DataVector& alpha,
                                      SGPP::base::DataVector& result,
                                      size_t dim, size_t op_dim) {
  if (dim == op_dim) {
    specialOP(alpha, result, dim, op_dim);
  } else {
    //Unidirectional scheme
    if (dim > 0) {
      // Reordering ups and downs
      SGPP::base::DataVector temp(alpha.getSize());
      temp.setAll(0.0);
      up(alpha, temp, dim);
      updown(temp, result, dim - 1, op_dim);

      // Same from the other direction:
      SGPP::base::DataVector result_temp(alpha.getSize());
      result_temp.setAll(0.0);
      updown(alpha, temp, dim - 1, op_dim);
      down(temp, result_temp, dim);

      result.add(result_temp);
    } else {
      // Terminates dimension recursion
      up(alpha, result, dim);

      SGPP::base::DataVector temp(alpha.getSize());
      temp.setAll(0.0);
      down(alpha, temp, dim);

      result.add(temp);
    }

  }
}

void UpDownOneOpDimWithShadow::specialOP(SGPP::base::DataVector& alpha,
    SGPP::base::DataVector& result,
    size_t dim, size_t op_dim) {
  //Unidirectional scheme
  if (dim > 0) {
    // Reordering ups and downs
    SGPP::base::DataVector temp(alpha.getSize());
    temp.setAll(0.0);
    upOpDim(alpha, temp, dim);
    updown(temp, result, dim - 1, op_dim);

    // Same from the other direction:
    SGPP::base::DataVector result_temp(alpha.getSize());
    result_temp.setAll(0.0);
    updown(alpha, temp, dim - 1, op_dim);
    downOpDim(temp, result_temp, dim);

    result.add(result_temp);
  } else {
    // Terminates dimension recursion
    upOpDim(alpha, result, dim);

    SGPP::base::DataVector temp(alpha.getSize());
    temp.setAll(0.0);
    downOpDim(alpha, temp, dim);

    result.add(temp);
  }
}

}
}