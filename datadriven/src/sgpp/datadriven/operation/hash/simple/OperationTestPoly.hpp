// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#ifndef OPERATIONTESTPOLY_HPP
#define OPERATIONTESTPOLY_HPP

#include <sgpp/base/operation/hash/OperationEval.hpp>
#include <sgpp/base/grid/GridStorage.hpp>
#include <sgpp/base/operation/hash/common/basis/PolyBasis.hpp>
#include <sgpp/datadriven/operation/hash/simple/OperationTest.hpp>
#include <sgpp/base/datatypes/DataVector.hpp>
#include <sgpp/base/datatypes/DataMatrix.hpp>

#include <sgpp/globaldef.hpp>

namespace sgpp {
namespace datadriven {

/**
 * This class implements OperationTest for a grids with poly basis ansatzfunctions with
 *
 */
class OperationTestPoly : public OperationTest {
 public:
  /**
   * Constructor
   *
   * @param storage the grid's base::GridStorage object
   * @param degree the polynom's max. degree
   */
  explicit OperationTestPoly(base::GridStorage* storage, size_t degree)
      : storage(storage), base(degree) {}

  /**
   * Destructor
   */
  virtual ~OperationTestPoly() {}

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
  /// Poly Basis object
  base::SPolyBase base;
};
}  // namespace datadriven
}  // namespace sgpp
#endif /* OPERATIONTESTPOLY_HPP */
