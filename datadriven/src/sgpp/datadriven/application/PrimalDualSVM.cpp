// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <sgpp/base/operation/BaseOpFactory.hpp>
#include <sgpp/base/operation/hash/OperationMultipleEval.hpp>
#include <sgpp/datadriven/application/PrimalDualSVM.hpp>

#include <sgpp/globaldef.hpp>

#include <cmath>

namespace sgpp {
namespace datadriven {

PrimalDualSVM::PrimalDualSVM(size_t dim, size_t dataDim, size_t budget,
                             bool useBias)
    : svs(sgpp::base::DataMatrix(0, dataDim)),
      alphas(sgpp::base::DataVector(0)),
      norms(sgpp::base::DataVector(0)),
      w(sgpp::base::DataVector(dim, 0.0)),
      w2(sgpp::base::DataVector(dim, 0.0)),
      budget(budget),
      useBias(useBias),
      bias(0.0) {
  svs.setInc(budget);
}

PrimalDualSVM::~PrimalDualSVM() {}

double PrimalDualSVM::predictRaw(sgpp::base::Grid& grid,
                                 sgpp::base::DataVector& x, size_t dataDim,
                                 bool trans) {
  sgpp::base::DataVector xTrans(grid.getSize());
  if (trans) {
    xTrans = x;
  } else {
    // SG-kernel evaluation
    sgpp::base::DataMatrix xMatrix(1, dataDim);
    xMatrix.setRow(0, x);
    sgpp::base::DataVector unitAlpha(1, 1.0);
    std::unique_ptr<base::OperationMultipleEval> multEval(
        op_factory::createOperationMultipleEval(grid, xMatrix));
    multEval->multTranspose(unitAlpha, xTrans);
  }
  double res = w.dotProduct(xTrans);
  if (useBias) {
    res += bias;
  }
  return res;
}

int PrimalDualSVM::predict(sgpp::base::Grid& grid, sgpp::base::DataVector& x,
                           size_t dataDim) {
  bool sign = std::signbit(this->predictRaw(grid, x, dataDim));
  if (sign) {
    return -1;
  } else {
    return 1;
  }
}

void PrimalDualSVM::multiply(double scalar) {
  if (scalar != 1.0) {
    w.mult(scalar);
    w2.mult(scalar);
    if (alphas.getSize() > 0) {
      alphas.mult(scalar);
    }
    if (useBias) {
      bias += scalar;
    }
  }
}

void PrimalDualSVM::add(sgpp::base::Grid& grid, sgpp::base::DataVector& x,
                        double alpha, size_t dataDim) {
  if (svs.getNrows() < budget) {
    sgpp::base::DataVector xTrans(grid.getSize());
    // SG-kernel evaluation
    sgpp::base::DataMatrix xMatrix(1, dataDim);
    xMatrix.setRow(0, x);
    sgpp::base::DataVector unitAlpha(1, 1.0);
    std::unique_ptr<base::OperationMultipleEval> multEval(
        op_factory::createOperationMultipleEval(grid, xMatrix));
    multEval->multTranspose(unitAlpha, xTrans);

    svs.appendRow(x);
    alphas.append(alpha);
    norms.append(xTrans.dotProduct(xTrans));

    sgpp::base::DataVector xTransCopy(xTrans);
    xTrans.mult(alpha);
    w.add(xTrans);
    xTransCopy.mult(std::abs(alpha));
    w2.add(xTransCopy);

    if (useBias) {
      bias += alpha;
    }
  }
}

/*void PrimalDualSVM::remove(size_t idx) {
  // ToDo:
}*/

}  // namespace datadriven
}  // namespace sgpp
