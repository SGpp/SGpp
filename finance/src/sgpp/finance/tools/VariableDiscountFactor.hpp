// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#ifndef VARIABLEDISCOUNTFACTOR_HPP
#define VARIABLEDISCOUNTFACTOR_HPP

#include <sgpp/base/grid/common/BoundingBox.hpp>
#include <sgpp/base/grid/GridStorage.hpp>
#include <sgpp/base/datatypes/DataVector.hpp>

#include <sgpp/globaldef.hpp>

namespace sgpp {
namespace finance {

/**
 * This class provides a method for discounting when the risk-free rate corresponds to an axis of
 * the grid
 */
class VariableDiscountFactor {
 private:
  /// bounding box of the grid
  sgpp::base::BoundingBox* myBoundingBox;
  /// sgpp::base::Grid Storage object
  sgpp::base::GridStorage* storage;
  /// dimension of the risk-free rate (HW axis)
  int dim_r;

 public:
  /**
   * Std-Constructor
   *
   * @param storage the grid's storage object; needed to determine the bounding box and to iterate
   * of the entries in the coefficient vector
   * @param dim_r dimension of the risk-free rate (HW axis)
   */
  VariableDiscountFactor(sgpp::base::GridStorage* storage, int dim_r);

  /**
   * Std-Destructor
   */
  ~VariableDiscountFactor();
  /**
   * get a vector which contains all the factors needed to multiply with another vector
   *@param factor the vector that should be calculated to multiply with another vector
   *@param T timestepsize
   */
  void getDiscountFactor(sgpp::base::DataVector& factor, double T);
};
}  // namespace finance
}  // namespace sgpp

#endif /* VARIABLEDISCOUNTFACTOR_HPP */
