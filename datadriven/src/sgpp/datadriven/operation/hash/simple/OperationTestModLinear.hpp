// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#ifndef OPERATIONTESTMODLINEAR_HPP
#define OPERATIONTESTMODLINEAR_HPP

#include <sgpp/datadriven/operation/hash/simple/OperationTest.hpp>
#include <sgpp/base/grid/GridStorage.hpp>

#include <sgpp/globaldef.hpp>

namespace sgpp {
namespace datadriven {

/**
 * This class implements sgpp::base::OperationEval for a grids with mod linear basis ansatzfunctions
 * with
 *
 */
class OperationTestModLinear : public OperationTest {
 public:
  /**
   * Constructor
   *
   * @param storage the grid's sgpp::base::GridStorage object
   */
  explicit OperationTestModLinear(sgpp::base::GridStorage* storage) : storage(storage) {}

  /**
   * Destructor
   */
  virtual ~OperationTestModLinear() {}

  virtual double test(sgpp::base::DataVector& alpha, sgpp::base::DataMatrix& data,
                       sgpp::base::DataVector& classes);
  virtual double testMSE(sgpp::base::DataVector& alpha, sgpp::base::DataMatrix& data,
                          sgpp::base::DataVector& refValues);
  virtual double testWithCharacteristicNumber(sgpp::base::DataVector& alpha,
                                               sgpp::base::DataMatrix& data,
                                               sgpp::base::DataVector& classes,
                                               sgpp::base::DataVector& charaNumbers);
  virtual void calculateROCcurve(sgpp::base::DataVector& alpha, sgpp::base::DataMatrix& data,
                                 sgpp::base::DataVector& classes,
                                 sgpp::base::DataVector& thresholds,
                                 sgpp::base::DataMatrix& ROC_curve);

 protected:
  /// Pointer to sgpp::base::GridStorage object
  sgpp::base::GridStorage* storage;
};

}  // namespace datadriven
}  // namespace sgpp

#endif /* OPERATIONTESTMODLINEAR_HPP */
