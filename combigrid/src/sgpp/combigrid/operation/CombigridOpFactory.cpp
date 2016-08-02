// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// SGPPpp.sparsegrids.org

#include <sgpp/combigrid/operation/CombigridOpFactory.hpp>
#include <sgpp/base/exception/factory_exception.hpp>
#include <sgpp/base/operation/hash/OperationQuadrature.hpp>

#include <cstring>
#include <sgpp/globaldef.hpp>

namespace SGPP {

  namespace op_factory {

//    base::OperationQuadrature* createOperationQuadrature(base::Grid& grid) {
//
//      if (grid.getType() == base::GridType::Linear) {
//        return new base::OperationQuadratureLinear(grid.getStorage());
//      } else if (grid.getType() == base::GridType::LinearL0Boundary
//                 || grid.getType() == base::GridType::LinearBoundary) {
//        return new base::OperationQuadratureLinearBoundary(grid.getStorage());
//      } else if (grid.getType() == base::GridType::Poly) {
//        return new base::OperationQuadraturePoly(grid.getStorage(),
//               dynamic_cast<base::PolyGrid*>(&grid)->getDegree());
//      } else if (grid.getType() == base::GridType::PolyBoundary) {
//        return new base::OperationQuadraturePolyBoundary(grid.getStorage(),
//               dynamic_cast<base::PolyBoundaryGrid*>(&grid)->getDegree());
//      } else
//        throw base::factory_exception(
//          "OperationQuadrature is not implemented for this grid type.");
//    }
//
//    base::OperationEval* createOperationEval(base::Grid& grid) {
//
//      if (grid.getType() == base::GridType::Linear) {
//        return new base::OperationEvalLinear(grid.getStorage());
//      } else if (grid.getType() == base::GridType::LinearL0Boundary
//                 || grid.getType() == base::GridType::LinearBoundary
//                 || grid.getType() == base::GridType::LinearTruncatedBoundary
//                 || grid.getType() == base::GridType::SquareRoot) {
//        return new base::OperationEvalLinearBoundary(grid.getStorage());
//      } else if (grid.getType() == base::GridType::ModLinear) {
//        return new base::OperationEvalModLinear(grid.getStorage());
//      } else if (grid.getType() == base::GridType::Poly) {
//        return new base::OperationEvalPoly(grid.getStorage(),
//                                           dynamic_cast<base::PolyGrid*>(&grid)->getDegree());
//      } else if (grid.getType() == base::GridType::PolyBoundary) {
//        return new base::OperationEvalPolyBoundary(grid.getStorage(),
//               dynamic_cast<base::PolyBoundaryGrid*>(&grid)->getDegree());
//      } else if (grid.getType() == base::GridType::ModPoly) {
//        return new base::OperationEvalModPoly(grid.getStorage(),
//                                              dynamic_cast<base::ModPolyGrid*>(&grid)->getDegree());
//      } else if (grid.getType() == base::GridType::ModBspline) {
//        return new base::OperationEvalModBspline(grid.getStorage(),
//               dynamic_cast<base::ModBsplineGrid*>(&grid)->getDegree());
//      } else if (grid.getType() == base::GridType::ModWavelet) {
//        return new base::OperationEvalModWavelet(grid.getStorage());
//      } else if (grid.getType() == base::GridType::Prewavelet) {
//        return new base::OperationEvalPrewavelet(grid.getStorage());
//      } else if (grid.getType() == base::GridType::LinearStretched) {
//        return new base::OperationEvalLinearStretched(grid.getStorage());
//      } else if (grid.getType() == base::GridType::LinearStretchedBoundary) {
//        return new base::OperationEvalLinearStretchedBoundary(grid.getStorage());
//      } else if (grid.getType() == base::GridType::Periodic) {
//        return new base::OperationEvalPeriodic(grid.getStorage());
//      } else
//        throw base::factory_exception(
//          "OperationEval is not implemented for this grid type.");
//    }
  }
}
