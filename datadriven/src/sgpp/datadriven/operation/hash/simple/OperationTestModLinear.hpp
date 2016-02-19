// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#ifndef OPERATIONTESTMODLINEAR_HPP
#define OPERATIONTESTMODLINEAR_HPP

#include <sgpp/datadriven/operation/hash/simple/OperationTest.hpp>
#include <sgpp/base/grid/GridStorage.hpp>

#include <sgpp/globaldef.hpp>

namespace SGPP {
namespace datadriven {

/**
 * This class implements SGPP::base::OperationEval for a grids with mod linear basis ansatzfunctions
 * with
 *
 */
class OperationTestModLinear : public OperationTest {
 public:
  /**
   * Constructor
   *
   * @param storage the grid's SGPP::base::GridStorage object
   */
  explicit OperationTestModLinear(SGPP::base::GridStorage* storage) : storage(storage) {}

  /**
   * Destructor
   */
  virtual ~OperationTestModLinear() {}

  virtual float_t test(SGPP::base::DataVector& alpha, SGPP::base::DataMatrix& data,
                       SGPP::base::DataVector& classes);
  virtual float_t testMSE(SGPP::base::DataVector& alpha, SGPP::base::DataMatrix& data,
                          SGPP::base::DataVector& refValues);
  virtual float_t testWithCharacteristicNumber(SGPP::base::DataVector& alpha,
                                               SGPP::base::DataMatrix& data,
                                               SGPP::base::DataVector& classes,
                                               SGPP::base::DataVector& charaNumbers);
  virtual void calculateROCcurve(SGPP::base::DataVector& alpha, SGPP::base::DataMatrix& data,
                                 SGPP::base::DataVector& classes,
                                 SGPP::base::DataVector& thresholds,
                                 SGPP::base::DataMatrix& ROC_curve);

 protected:
  /// Pointer to SGPP::base::GridStorage object
  SGPP::base::GridStorage* storage;
};

}  // namespace datadriven
}  // namespace SGPP

#endif /* OPERATIONTESTMODLINEAR_HPP */
