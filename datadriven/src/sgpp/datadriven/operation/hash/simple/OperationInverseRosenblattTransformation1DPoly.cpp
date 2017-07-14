// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <sgpp/datadriven/operation/hash/simple/OperationInverseRosenblattTransformation1DPoly.hpp>
#include <sgpp/datadriven/operation/hash/simple/OperationRosenblattTransformation1DPoly.hpp>
#include <sgpp/base/exception/operation_exception.hpp>
#include <sgpp/base/exception/algorithm_exception.hpp>
#include <sgpp/base/operation/hash/OperationEval.hpp>
#include <sgpp/datadriven/DatadrivenOpFactory.hpp>
#include <sgpp/base/datatypes/DataVector.hpp>
#include <sgpp/base/grid/type/PolyGrid.hpp>
#include <sgpp/base/tools/GaussLegendreQuadRule1D.hpp>
#include <sgpp_optimization.hpp>
#include <sgpp_datadriven.hpp>

#include <sgpp/globaldef.hpp>
#include <map>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <utility>
#include <vector>
#include <algorithm>
#include <functional>

namespace sgpp {
namespace datadriven {
/**
 * WARNING: the grid must be a 1D grid!
 */
OperationInverseRosenblattTransformation1DPoly::OperationInverseRosenblattTransformation1DPoly(
  base::Grid* grid)
  : grid(grid) {}

OperationInverseRosenblattTransformation1DPoly::~OperationInverseRosenblattTransformation1DPoly() {}

double OperationInverseRosenblattTransformation1DPoly::doTransformation1D(base::DataVector* alpha1d,
                                                                          double coord1d) {

  std::cout << "INVERSE" << std::endl;
  OperationRosenblattTransformation1DPoly* opRosenblatt
    = static_cast<OperationRosenblattTransformation1DPoly*>
    (op_factory::createOperationRosenblattTransformation1D(*grid));

  std::function<double(const base::DataVector&)> optFunc =
    [coord1d, opRosenblatt, alpha1d](const base::DataVector& x) -> double {
    double F_x = opRosenblatt->doTransformation1D(alpha1d, x[0]);
    return  (F_x - coord1d) * (F_x - coord1d);};

  optimization::WrapperScalarFunction f(1, optFunc);
  optimization::optimizer::NelderMead nelderMead(f);
  nelderMead.optimize();
  return nelderMead.getOptimalPoint()[0];
}
}  // namespace datadriven
}  // namespace sgpp
