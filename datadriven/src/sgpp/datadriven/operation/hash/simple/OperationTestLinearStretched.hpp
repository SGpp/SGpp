// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#ifndef OPERATIONTESTLINEARSTRETCHED_HPP
#define OPERATIONTESTLINEARSTRETCHED_HPP

#include <sgpp/datadriven/operation/hash/simple/OperationTest.hpp>
#include <sgpp/base/grid/GridStorage.hpp>

#include <sgpp/globaldef.hpp>

namespace sgpp {
namespace datadriven {

/**
 * This class implements OperationTest for a grids with linearstretched basis ansatzfunctions
 * without boundaries
 */
class OperationTestLinearStretched : public OperationTest {
 public:
  /**
   * Constructor of OperationTestLinearStretched
   *
   * @param storage Pointer to the grid's gridstorage obejct
   */
  explicit OperationTestLinearStretched(sgpp::base::GridStorage* storage) : storage(storage) {}

  /**
   * Destructor
   */
  virtual ~OperationTestLinearStretched() {}

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
  /// Pointer to the grid's GridStorage object
  sgpp::base::GridStorage* storage;
};

}  // namespace datadriven
}  // namespace sgpp

#endif /* OPERATIONTESTLINEARSTRETCHED_HPP */
