// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#ifndef OPERATIONEVALMODBSPLINE_HPP
#define OPERATIONEVALMODBSPLINE_HPP

#include <sgpp/base/operation/hash/OperationEval.hpp>
#include <sgpp/base/grid/GridStorage.hpp>
#include <sgpp/base/operation/hash/common/basis/BsplineModifiedBasis.hpp>
#include <sgpp/base/datatypes/DataVector.hpp>

#include <sgpp/globaldef.hpp>


namespace SGPP {
namespace base {

/**
 * This class implements OperationEval for a grids with modified Bspline basis functions with a certain degree
 *
 */
class OperationEvalModBspline : public OperationEval {
 public:
  /**
   * Constructor
   *
   * @param storage the grid's GridStorage object
   * @param degree the polynom's max. degree
   */
  OperationEvalModBspline(GridStorage& storage, size_t degree) :
    storage(storage),
    base(degree) {}

  /**
   * Destructor
   */
  ~OperationEvalModBspline() override {}

  float_t eval(const DataVector& alpha,
               const DataVector& point) override;

 protected:
  /// Pointer to GridStorage object
  GridStorage& storage;
  /// Mod Bspline Basis object
  SBsplineModifiedBase base;
};

}  // namespace base
}  // namespace SGPP

#endif /* OPERATIONEVALMODBSPLINE_HPP */
