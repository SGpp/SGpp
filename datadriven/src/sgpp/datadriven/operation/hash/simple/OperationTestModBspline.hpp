// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#ifndef OPERATIONTESTMODBSPLINE_HPP
#define OPERATIONTESTMODBSPLINE_HPP

#include <sgpp/datadriven/operation/hash/simple/OperationTest.hpp>
#include <sgpp/base/grid/GridStorage.hpp>

#include <sgpp/base/operation/hash/common/basis/BsplineModifiedBasis.hpp>

#include <sgpp/globaldef.hpp>

namespace sgpp {
namespace datadriven {

/**
 * This class implements OperationTest for a grids with modified bspline basis functions with a
 * certain degree
 *
 */
class OperationTestModBspline : public OperationTest {
 public:
  /**
   * Constructor
   *
   * @param storage the grid's base::GridStorage object
   * @param degree the bspline's degree
   */
  explicit OperationTestModBspline(base::GridStorage* storage, size_t degree)
      : storage(storage), base(degree) {}

  /**
   * Destructor
   */
  virtual ~OperationTestModBspline() {}

  virtual double test(base::DataVector& alpha, base::DataMatrix& data, base::DataVector& classes);
  virtual double testMSE(base::DataVector& alpha, base::DataMatrix& data,
                          base::DataVector& refValues);
  virtual double testWithCharacteristicNumber(sgpp::base::DataVector& alpha,
                                               sgpp::base::DataMatrix& data,
                                               sgpp::base::DataVector& classes,
                                               sgpp::base::DataVector& charaNumbers);
  virtual void calculateROCcurve(sgpp::base::DataVector& alpha, sgpp::base::DataMatrix& data,
                                 sgpp::base::DataVector& classes,
                                 sgpp::base::DataVector& thresholds,
                                 sgpp::base::DataMatrix& ROC_curve);

 protected:
  /// Pointer to base::GridStorage object
  base::GridStorage* storage;
  /// Mod Bspline Basis object
  base::SBsplineModifiedBase base;
};

}  // namespace datadriven
}  // namespace sgpp

#endif /* OPERATIONTESTMODBSPLINE_HPP */
