/* ****************************************************************************
* Copyright (C) 2009 Technische Universitaet Muenchen                         *
* This file is part of the SG++ project. For conditions of distribution and   *
* use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
**************************************************************************** */
// @author Jörg Blank (blankj@in.tum.de), Alexander Heinecke (Alexander.Heinecke@mytum.de)

#ifndef OPERATIONHIERARCHISATIONPOLYBOUNDARY_HPP
#define OPERATIONHIERARCHISATIONPOLYBOUNDARY_HPP

#include <sgpp/base/operation/hash/OperationHierarchisation.hpp>
#include <sgpp/base/grid/GridStorage.hpp>
#include <sgpp/base/operation/hash/common/basis/PolyBoundaryBasis.hpp>
#include <sgpp/base/datatypes/DataVector.hpp>

#include <sgpp/globaldef.hpp>


namespace SGPP {
namespace base {

/**
 * Hierarchisation on sparse grid, poly case
 */
class OperationHierarchisationPolyBoundary : public OperationHierarchisation {
 public:
  /**
   * Constructor
   *
   * @param storage the grid's GridStorage object
   * @param degree the polynom's max. degree
   */
  OperationHierarchisationPolyBoundary(GridStorage* storage,
                                       size_t degree) : storage(storage), base(degree) {}

  /**
   * Destructor
   */
  virtual ~OperationHierarchisationPolyBoundary() override {}

  /**
   * Implements the hierarchisation on a sprase grid with poly base functions
   *
   * @param node_values the functions values in the node base
   *
   */
  virtual void doHierarchisation(DataVector& node_values) override;

  /**
   * Implements the dehierarchisation on a sprase grid with poly base functions
   *
   * @param alpha the coefficients of the sparse grid's base functions
   *
   */
  virtual void doDehierarchisation(DataVector& alpha) override;

 protected:
  /// Pointer to GridStorage object
  GridStorage* storage;
  /// Poly Basis object
  SPolyBoundaryBase base;
};

}  // namespace base
}  // namespace SGPP

#endif /* OPERATIONHIERARCHISATIONPOLYBOUNDARY_HPP */
