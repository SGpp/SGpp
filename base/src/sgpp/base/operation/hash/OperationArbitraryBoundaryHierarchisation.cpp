// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <sgpp/globaldef.hpp>
#include <sgpp/base/operation/hash/OperationArbitraryBoundaryHierarchisation.hpp>
#include <sgpp/base/operation/hash/OperationHierarchisation.hpp>
#include <sgpp/base/operation/BaseOpFactory.hpp>
#include <sgpp/base/exception/factory_exception.hpp>
#include <sgpp/base/exception/algorithm_exception.hpp>

namespace sgpp {
namespace base {

OperationArbitraryBoundaryHierarchisation::OperationArbitraryBoundaryHierarchisation(Grid& grid)
    : grid(grid) {
  // init inner grid
  RegularGridConfiguration innerGridConfig;
  innerGridConfig.dim_ = grid.getDimension();
  innerGridConfig.type_ = grid.getZeroBoundaryType();
  innerGridConfig.maxDegree_ = grid.getBasis().getDegree();
  innerGridConfig.boundaryLevel_ = 0;

  innerGrid.reset(Grid::createGrid(innerGridConfig));
  HashGridStorage& innerGs = innerGrid->getStorage();
  // init outer grid
  RegularGridConfiguration boundaryGridConfig;
  boundaryGridConfig.dim_ = grid.getDimension();
  boundaryGridConfig.type_ = grid.getType();
  boundaryGridConfig.maxDegree_ = grid.getBasis().getDegree();
  boundaryGridConfig.boundaryLevel_ = 0;

  boundaryGrid.reset(Grid::createGrid(boundaryGridConfig));
  HashGridStorage& boundaryGs = boundaryGrid->getStorage();

  // split inner and boundary grid
  HashGridStorage& gs = grid.getStorage();
  for (size_t i = 0; i < gs.getSize(); i++) {
    HashGridPoint gp = gs.getPoint(i);
    if (gp.isInnerPoint()) {
      innerGs.insert(gp);
    } else {
      boundaryGs.insert(gp);
    }
  }

  // check if this operation is really needed
  if (innerGrid->getSize() == 0) {
    throw algorithm_exception(
        "there are no inner points available. Use the standard hierarchisation operation "
        "'base::op_factory::createOperationHierarchisation' instead");
  }

  if (boundaryGrid->getSize() == 0) {
    throw algorithm_exception(
        "there are no boundary points available. Use the standard hierarchisation operation "
        "'base::op_factory::createOperationHierarchisation' instead");
  }
}

OperationArbitraryBoundaryHierarchisation::~OperationArbitraryBoundaryHierarchisation() {}

OperationMultipleEval* OperationArbitraryBoundaryHierarchisation::createOperationMultipleEval(
    Grid& grid, DataMatrix& coordinates) {
  OperationMultipleEval* opEval;
  try {
    opEval = op_factory::createOperationMultipleEval(grid, coordinates);
  } catch (base::factory_exception&) {
    try {
      opEval = op_factory::createOperationMultipleEvalNaive(grid, coordinates);
    } catch (base::factory_exception&) {
      throw;
    }
  }

  return opEval;
}

void OperationArbitraryBoundaryHierarchisation::doHierarchisation(DataVector& nodal_values) {
  // collect nodal values for boundary grids and inner grids
  DataVector boundaryNodalValues(boundaryGrid->getSize());
  DataVector innerNodalValues(innerGrid->getSize());

  // split nodal values for inner and boundary grid
  HashGridStorage& gs = grid.getStorage();
  HashGridStorage& innerGs = innerGrid->getStorage();
  HashGridStorage& boundaryGs = boundaryGrid->getStorage();

  size_t j = 0;
  for (size_t i = 0; i < gs.getSize(); i++) {
    HashGridPoint gp = gs.getPoint(i);
    if (gp.isInnerPoint()) {
      j = innerGs.getSequenceNumber(gp);
      innerNodalValues[j] = nodal_values[i];
    } else {
      j = boundaryGs.getSequenceNumber(gp);
      boundaryNodalValues[j] = nodal_values[i];
    }
  }

  // hierarchize boundary grid with standard hierarchization
  std::unique_ptr<OperationHierarchisation> boundaryHierOp(
      op_factory::createOperationHierarchisation(*boundaryGrid));
  boundaryHierOp->doHierarchisation(boundaryNodalValues);

  // subtract the function values from the boundary grid from the inner nodal values
  DataMatrix coordinates;
  DataVector boundaryValues(innerGs.getSize());
  innerGs.getCoordinateArrays(coordinates);
  std::unique_ptr<OperationMultipleEval> opEval(
      createOperationMultipleEval(*boundaryGrid, coordinates));
  opEval->eval(boundaryNodalValues, boundaryValues);
  innerNodalValues.sub(boundaryValues);

  // hierarchize the inner grid with the updated nodal values
  std::unique_ptr<OperationHierarchisation> innerHierOp(
      op_factory::createOperationHierarchisation(*innerGrid));
  innerHierOp->doHierarchisation(innerNodalValues);

  // join coefficients
  for (size_t i = 0; i < gs.getSize(); i++) {
    HashGridPoint gp = gs.getPoint(i);
    if (gp.isInnerPoint()) {
      j = innerGs.getSequenceNumber(gp);
      nodal_values[i] = innerNodalValues[j];
    } else {
      j = boundaryGs.getSequenceNumber(gp);
      nodal_values[i] = boundaryNodalValues[j];
    }
  }
}

void OperationArbitraryBoundaryHierarchisation::doDehierarchisation(DataVector& alpha) {
  // split coefficients for inner and boundary grid
  HashGridStorage& gs = grid.getStorage();
  HashGridStorage& innerGs = innerGrid->getStorage();
  HashGridStorage& boundaryGs = boundaryGrid->getStorage();

  DataVector boundaryAlpha(boundaryGrid->getSize());
  DataVector innerAlpha(innerGrid->getSize());

  size_t j = 0;
  for (size_t i = 0; i < gs.getSize(); i++) {
    HashGridPoint gp = gs.getPoint(i);
    if (gp.isInnerPoint()) {
      j = innerGs.getSequenceNumber(gp);
      innerAlpha[j] = alpha[i];
    } else {
      j = boundaryGs.getSequenceNumber(gp);
      boundaryAlpha[j] = alpha[i];
    }
  }

  // evaluate the inner and the outer grid at the actual grid points
  DataMatrix coordinates;
  gs.getCoordinateArrays(coordinates);

  // evaluate the boundary grid
  DataVector boundaryValues(gs.getSize());
  std::unique_ptr<OperationMultipleEval> opEvalBoundary(
      createOperationMultipleEval(*boundaryGrid, coordinates));
  opEvalBoundary->eval(boundaryAlpha, boundaryValues);

  // evaluate the inner grid
  DataVector innerValues(gs.getSize());
  std::unique_ptr<OperationMultipleEval> opEvalInner(
      createOperationMultipleEval(*innerGrid, coordinates));
  opEvalInner->eval(innerAlpha, innerValues);

  for (size_t i = 0; i < gs.getSize(); i++) {
    alpha[i] = innerValues[i] + boundaryValues[i];
  }
}

}  // namespace base
}  // namespace sgpp
