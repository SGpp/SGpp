// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#ifndef OPERATIONEVALMODPOLY_HPP
#define OPERATIONEVALMODPOLY_HPP

#include <sgpp/base/operation/hash/OperationEval.hpp>
#include <sgpp/base/grid/GridStorage.hpp>
#include <sgpp/base/operation/hash/common/basis/PolyModifiedBasis.hpp>
#include <sgpp/base/datatypes/DataVector.hpp>

#include <sgpp/globaldef.hpp>


namespace SGPP {
namespace base {

/**
 * This class implements OperationEval for a grids with mod poly basis ansatzfunctions with
 *
 */
class OperationEvalModPoly : public OperationEval {
 public:
  /**
   * Constructor
   *
   * @param storage the grid's GridStorage object
   * @param degree the polynom's max. degree
   */
  OperationEvalModPoly(GridStorage* storage, size_t degree) : storage(storage),
    base(degree) {}

  /**
   * Destructor
   */
  virtual ~OperationEvalModPoly() override {}

  virtual float_t eval(const DataVector& alpha,
                       const DataVector& point) override;

 protected:
  /// Pointer to GridStorage object
  GridStorage* storage;
  /// Mod Poly Basis object
  SGPP::base::SPolyModifiedBase base;
};

}  // namespace base
}  // namespace SGPP

#endif /* OPERATIONEVALMODPOLY_HPP */
