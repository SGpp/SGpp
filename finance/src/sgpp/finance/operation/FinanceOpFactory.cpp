// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <sgpp/finance/operation/FinanceOpFactory.hpp>

#include <sgpp/base/exception/factory_exception.hpp>

#include <sgpp/finance/operation/hash/OperationGammaLinear.hpp>
#include <sgpp/finance/operation/hash/OperationGammaLinearStretched.hpp>
#include <sgpp/finance/operation/hash/OperationGammaLogLinearStretchedBoundary.hpp>
#include <sgpp/finance/operation/hash/OperationGammaLinearBoundary.hpp>

#include <sgpp/finance/operation/hash/OperationGammaLogLinear.hpp>
#include <sgpp/finance/operation/hash/OperationGammaLogLinearStretched.hpp>
#include <sgpp/finance/operation/hash/OperationGammaLinearStretchedBoundary.hpp>
#include <sgpp/finance/operation/hash/OperationGammaLogLinearBoundary.hpp>

#include <sgpp/finance/operation/hash/OperationLBLinear.hpp>
#include <sgpp/finance/operation/hash/OperationLBLinearBoundary.hpp>

#include <sgpp/finance/operation/hash/OperationLELinear.hpp>
#include <sgpp/finance/operation/hash/OperationLELinearBoundary.hpp>

#include <sgpp/finance/operation/hash/OperationLDLinear.hpp>
#include <sgpp/finance/operation/hash/OperationLDLinearBoundary.hpp>

#include <sgpp/finance/operation/hash/OperationLFLinear.hpp>
#include <sgpp/finance/operation/hash/OperationLFLinearBoundary.hpp>

#include <sgpp/finance/operation/hash/OperationDeltaLinear.hpp>
#include <sgpp/finance/operation/hash/OperationDeltaLinearStretched.hpp>
#include <sgpp/finance/operation/hash/OperationDeltaLinearStretchedBoundary.hpp>
#include <sgpp/finance/operation/hash/OperationDeltaLinearBoundary.hpp>

#include <sgpp/finance/operation/hash/OperationDeltaLogLinear.hpp>
#include <sgpp/finance/operation/hash/OperationDeltaLogLinearStretched.hpp>
#include <sgpp/finance/operation/hash/OperationDeltaLogLinearStretchedBoundary.hpp>
#include <sgpp/finance/operation/hash/OperationDeltaLogLinearBoundary.hpp>

#include <sgpp/finance/operation/hash/OperationHestonBLinear.hpp>
#include <sgpp/finance/operation/hash/OperationHestonCLinear.hpp>
#include <sgpp/finance/operation/hash/OperationHestonDLinear.hpp>
#include <sgpp/finance/operation/hash/OperationHestonELinear.hpp>
#include <sgpp/finance/operation/hash/OperationHestonFLinear.hpp>
#include <sgpp/finance/operation/hash/OperationHestonGLinear.hpp>
#include <sgpp/finance/operation/hash/OperationHestonHLinear.hpp>
#include <sgpp/finance/operation/hash/OperationHestonKLinear.hpp>
#include <sgpp/finance/operation/hash/OperationHestonWLinear.hpp>
#include <sgpp/finance/operation/hash/OperationHestonXLinear.hpp>
#include <sgpp/finance/operation/hash/OperationHestonYLinear.hpp>
#include <sgpp/finance/operation/hash/OperationHestonZLinear.hpp>

#include <sgpp/finance/operation/hash/OperationHestonBLinearBoundary.hpp>
#include <sgpp/finance/operation/hash/OperationHestonCLinearBoundary.hpp>
#include <sgpp/finance/operation/hash/OperationHestonDLinearBoundary.hpp>
#include <sgpp/finance/operation/hash/OperationHestonELinearBoundary.hpp>
#include <sgpp/finance/operation/hash/OperationHestonFLinearBoundary.hpp>
#include <sgpp/finance/operation/hash/OperationHestonGLinearBoundary.hpp>
#include <sgpp/finance/operation/hash/OperationHestonHLinearBoundary.hpp>
#include <sgpp/finance/operation/hash/OperationHestonKLinearBoundary.hpp>
#include <sgpp/finance/operation/hash/OperationHestonWLinearBoundary.hpp>
#include <sgpp/finance/operation/hash/OperationHestonXLinearBoundary.hpp>
#include <sgpp/finance/operation/hash/OperationHestonYLinearBoundary.hpp>
#include <sgpp/finance/operation/hash/OperationHestonZLinearBoundary.hpp>

#include <sgpp/globaldef.hpp>

#include <cstring>

namespace sgpp {

namespace op_factory {

std::unique_ptr<base::OperationMatrix> createOperationGamma(
    base::Grid& grid, base::DataMatrix& coef) {
  if (grid.getType() == base::GridType::Linear) {
    return std::unique_ptr<base::OperationMatrix>(
        new finance::OperationGammaLinear(&grid.getStorage(), coef));
  } else if (grid.getType() == base::GridType::LinearStretched) {
    return std::unique_ptr<base::OperationMatrix>(
        new finance::OperationGammaLinearStretched(&grid.getStorage(), coef));
  } else if (grid.getType() == base::GridType::LinearStretchedBoundary) {
    return std::unique_ptr<base::OperationMatrix>(
        new finance::OperationGammaLinearStretchedBoundary(&grid.getStorage(), coef));
  } else if (grid.getType() == base::GridType::LinearL0Boundary ||
             grid.getType() == base::GridType::LinearBoundary) {
    return std::unique_ptr<base::OperationMatrix>(
        new finance::OperationGammaLinearBoundary(&grid.getStorage(), coef));
  } else {
    throw base::factory_exception("OperationGamma is not implemented for this grid type.");
  }
}

std::unique_ptr<base::OperationMatrix> createOperationGammaLog(
    base::Grid& grid, base::DataMatrix& coef) {
  if (grid.getType() == base::GridType::Linear) {
    return std::unique_ptr<base::OperationMatrix>(
        new finance::OperationGammaLogLinear(&grid.getStorage(), coef));
  } else if (grid.getType() == base::GridType::LinearStretched) {
    return std::unique_ptr<base::OperationMatrix>(
        new finance::OperationGammaLogLinearStretched(&grid.getStorage(), coef));
  } else if (grid.getType() == base::GridType::LinearStretchedBoundary) {
    return std::unique_ptr<base::OperationMatrix>(
        new finance::OperationGammaLogLinearStretchedBoundary(&grid.getStorage(), coef));
  } else if (grid.getType() == base::GridType::LinearL0Boundary ||
             grid.getType() == base::GridType::LinearBoundary) {
    return std::unique_ptr<base::OperationMatrix>(
        new finance::OperationGammaLogLinearBoundary(&grid.getStorage(), coef));
  } else {
    throw base::factory_exception("OperationGammaLog is not implemented for this grid type.");
  }
}

std::unique_ptr<base::OperationMatrix> createOperationLB(base::Grid& grid) {
  if (grid.getType() == base::GridType::Linear) {
    return std::unique_ptr<base::OperationMatrix>(
        new finance::OperationLBLinear(&grid.getStorage()));
  } else if (grid.getType() == base::GridType::LinearL0Boundary ||
             grid.getType() == base::GridType::LinearBoundary) {
    return std::unique_ptr<base::OperationMatrix>(
        new finance::OperationLBLinearBoundary(&grid.getStorage()));
  } else {
    throw base::factory_exception("OperationLB is not implemented for this grid type.");
  }
}

std::unique_ptr<base::OperationMatrix> createOperationLE(base::Grid& grid) {
  if (grid.getType() == base::GridType::Linear) {
    return std::unique_ptr<base::OperationMatrix>(
        new finance::OperationLELinear(&grid.getStorage()));
  } else if (grid.getType() == base::GridType::LinearL0Boundary ||
             grid.getType() == base::GridType::LinearBoundary) {
    return std::unique_ptr<base::OperationMatrix>(
        new finance::OperationLELinearBoundary(&grid.getStorage()));
  } else {
    throw base::factory_exception("OperationLE is not implemented for this grid type.");
  }
}

std::unique_ptr<base::OperationMatrix> createOperationLD(base::Grid& grid) {
  if (grid.getType() == base::GridType::Linear) {
    return std::unique_ptr<base::OperationMatrix>(
        new finance::OperationLDLinear(&grid.getStorage()));
  } else if (grid.getType() == base::GridType::LinearL0Boundary ||
             grid.getType() == base::GridType::LinearBoundary) {
    return std::unique_ptr<base::OperationMatrix>(
        new finance::OperationLDLinearBoundary(&grid.getStorage()));
  } else {
    throw base::factory_exception("OperationLD is not implemented for this grid type.");
  }
}

std::unique_ptr<base::OperationMatrix> createOperationLF(base::Grid& grid) {
  if (grid.getType() == base::GridType::Linear) {
    return std::unique_ptr<base::OperationMatrix>(
        new finance::OperationLFLinear(&grid.getStorage()));
  } else if (grid.getType() == base::GridType::LinearL0Boundary ||
             grid.getType() == base::GridType::LinearBoundary) {
    return std::unique_ptr<base::OperationMatrix>(
        new finance::OperationLFLinearBoundary(&grid.getStorage()));
  } else {
    throw base::factory_exception("OperationLF is not implemented for this grid type.");
  }
}

std::unique_ptr<base::OperationMatrix> createOperationDelta(
    base::Grid& grid, base::DataVector& coef) {
  if (grid.getType() == base::GridType::Linear) {
    return std::unique_ptr<base::OperationMatrix>(
        new finance::OperationDeltaLinear(&grid.getStorage(), coef));
  } else if (grid.getType() == base::GridType::LinearStretched) {
    return std::unique_ptr<base::OperationMatrix>(
        new finance::OperationDeltaLinearStretched(&grid.getStorage(), coef));
  } else if (grid.getType() == base::GridType::LinearStretchedBoundary) {
    return std::unique_ptr<base::OperationMatrix>(
        new finance::OperationDeltaLinearStretchedBoundary(&grid.getStorage(), coef));
  } else if (grid.getType() == base::GridType::LinearL0Boundary ||
             grid.getType() == base::GridType::LinearBoundary) {
    return std::unique_ptr<base::OperationMatrix>(
        new finance::OperationDeltaLinearBoundary(&grid.getStorage(), coef));
  } else {
    throw base::factory_exception("OperationDelta is not implemented for this grid type.");
  }
}

std::unique_ptr<base::OperationMatrix> createOperationDeltaLog(
    base::Grid& grid, base::DataVector& coef) {
  if (grid.getType() == base::GridType::Linear) {
    return std::unique_ptr<base::OperationMatrix>(
        new finance::OperationDeltaLogLinear(&grid.getStorage(), coef));
  } else if (grid.getType() == base::GridType::LinearStretched) {
    return std::unique_ptr<base::OperationMatrix>(
        new finance::OperationDeltaLogLinearStretched(&grid.getStorage(), coef));
  } else if (grid.getType() == base::GridType::LinearStretchedBoundary) {
    return std::unique_ptr<base::OperationMatrix>(
        new finance::OperationDeltaLogLinearStretchedBoundary(&grid.getStorage(), coef));
  } else if (grid.getType() == base::GridType::LinearL0Boundary ||
             grid.getType() == base::GridType::LinearBoundary) {
    return std::unique_ptr<base::OperationMatrix>(
        new finance::OperationDeltaLogLinearBoundary(&grid.getStorage(), coef));
  } else {
    throw base::factory_exception("OperationDeltaLog is not implemented for this grid type.");
  }
}

std::unique_ptr<base::OperationMatrix> createOperationHestonBLog(
    base::Grid& grid, base::DataMatrix& coef) {
  if (grid.getType() == base::GridType::Linear) {
    return std::unique_ptr<base::OperationMatrix>(
        new finance::OperationHestonBLinear(&grid.getStorage(), coef));
  } else if (grid.getType() == base::GridType::LinearL0Boundary ||
             grid.getType() == base::GridType::LinearBoundary) {
    return std::unique_ptr<base::OperationMatrix>(
        new finance::OperationHestonBLinearBoundary(&grid.getStorage(), coef));
  } else {
    throw base::factory_exception("OperationHestonB is not implemented for this grid type.");
  }
}

std::unique_ptr<base::OperationMatrix> createOperationHestonCLog(
    base::Grid& grid, base::DataMatrix& coef) {
  if (grid.getType() == base::GridType::Linear) {
    return std::unique_ptr<base::OperationMatrix>(
        new finance::OperationHestonCLinear(&grid.getStorage(), coef));
  } else if (grid.getType() == base::GridType::LinearL0Boundary ||
             grid.getType() == base::GridType::LinearBoundary) {
    return std::unique_ptr<base::OperationMatrix>(
        new finance::OperationHestonCLinearBoundary(&grid.getStorage(), coef));
  } else {
    throw base::factory_exception("OperationHestonC is not implemented for this grid type.");
  }
}

std::unique_ptr<base::OperationMatrix> createOperationHestonDLog(
    base::Grid& grid, base::DataVector& coef) {
  if (grid.getType() == base::GridType::Linear) {
    return std::unique_ptr<base::OperationMatrix>(
        new finance::OperationHestonDLinear(&grid.getStorage(), coef));
  } else if (grid.getType() == base::GridType::LinearL0Boundary ||
             grid.getType() == base::GridType::LinearBoundary) {
    return std::unique_ptr<base::OperationMatrix>(
        new finance::OperationHestonDLinearBoundary(&grid.getStorage(), coef));
  } else {
    throw base::factory_exception("OperationHestonD is not implemented for this grid type.");
  }
}

std::unique_ptr<base::OperationMatrix> createOperationHestonELog(
    base::Grid& grid, base::DataVector& coef) {
  if (grid.getType() == base::GridType::Linear) {
    return std::unique_ptr<base::OperationMatrix>(
        new finance::OperationHestonELinear(&grid.getStorage(), coef));
  } else if (grid.getType() == base::GridType::LinearL0Boundary ||
             grid.getType() == base::GridType::LinearBoundary) {
    return std::unique_ptr<base::OperationMatrix>(
        new finance::OperationHestonELinearBoundary(&grid.getStorage(), coef));
  } else {
    throw base::factory_exception("OperationHestonE is not implemented for this grid type.");
  }
}

std::unique_ptr<base::OperationMatrix> createOperationHestonFLog(
    base::Grid& grid, base::DataVector& coef) {
  if (grid.getType() == base::GridType::Linear) {
    return std::unique_ptr<base::OperationMatrix>(
        new finance::OperationHestonFLinear(&grid.getStorage(), coef));
  } else if (grid.getType() == base::GridType::LinearL0Boundary ||
             grid.getType() == base::GridType::LinearBoundary) {
    return std::unique_ptr<base::OperationMatrix>(
        new finance::OperationHestonFLinearBoundary(&grid.getStorage(), coef));
  } else {
    throw base::factory_exception("OperationHestonF is not implemented for this grid type.");
  }
}

std::unique_ptr<base::OperationMatrix> createOperationHestonGLog(
    base::Grid& grid, base::DataVector& coef) {
  if (grid.getType() == base::GridType::Linear) {
    return std::unique_ptr<base::OperationMatrix>(
        new finance::OperationHestonGLinear(&grid.getStorage(), coef));
  } else if (grid.getType() == base::GridType::LinearL0Boundary ||
             grid.getType() == base::GridType::LinearBoundary) {
    return std::unique_ptr<base::OperationMatrix>(
        new finance::OperationHestonGLinearBoundary(&grid.getStorage(), coef));
  } else {
    throw base::factory_exception("OperationHestonG is not implemented for this grid type.");
  }
}

std::unique_ptr<base::OperationMatrix> createOperationHestonHLog(
    base::Grid& grid, base::DataMatrix& coef) {
  if (grid.getType() == base::GridType::Linear) {
    return std::unique_ptr<base::OperationMatrix>(
        new finance::OperationHestonHLinear(&grid.getStorage(), coef));
  } else if (grid.getType() == base::GridType::LinearL0Boundary ||
             grid.getType() == base::GridType::LinearBoundary) {
    return std::unique_ptr<base::OperationMatrix>(
        new finance::OperationHestonHLinearBoundary(&grid.getStorage(), coef));
  } else {
    throw base::factory_exception("OperationHestonH is not implemented for this grid type.");
  }
}

std::unique_ptr<base::OperationMatrix> createOperationHestonKLog(
    base::Grid& grid, double***** coef) {
  if (grid.getType() == base::GridType::Linear) {
    return std::unique_ptr<base::OperationMatrix>(
        new finance::OperationHestonKLinear(&grid.getStorage(), coef));
  } else if (grid.getType() == base::GridType::LinearL0Boundary ||
             grid.getType() == base::GridType::LinearBoundary) {
    return std::unique_ptr<base::OperationMatrix>(
        new finance::OperationHestonKLinearBoundary(&grid.getStorage(), coef));
  } else {
    throw base::factory_exception("OperationHestonK is not implemented for this grid type.");
  }
}

std::unique_ptr<base::OperationMatrix> createOperationHestonX(
    base::Grid& grid, base::DataMatrix& coef) {
  if (grid.getType() == base::GridType::Linear) {
    return std::unique_ptr<base::OperationMatrix>(
        new finance::OperationHestonXLinear(&grid.getStorage(), coef));
  } else if (grid.getType() == base::GridType::LinearL0Boundary ||
             grid.getType() == base::GridType::LinearBoundary) {
    return std::unique_ptr<base::OperationMatrix>(
        new finance::OperationHestonXLinearBoundary(&grid.getStorage(), coef));
  } else {
    throw base::factory_exception("OperationHestonX is not implemented for this grid type.");
  }
}

std::unique_ptr<base::OperationMatrix> createOperationHestonY(
    base::Grid& grid, base::DataMatrix& coef) {
  if (grid.getType() == base::GridType::Linear) {
    return std::unique_ptr<base::OperationMatrix>(
        new finance::OperationHestonYLinear(&grid.getStorage(), coef));
  } else if (grid.getType() == base::GridType::LinearL0Boundary ||
             grid.getType() == base::GridType::LinearBoundary) {
    return std::unique_ptr<base::OperationMatrix>(
        new finance::OperationHestonYLinearBoundary(&grid.getStorage(), coef));
  } else {
    throw base::factory_exception("OperationHestonY is not implemented for this grid type.");
  }
}

std::unique_ptr<base::OperationMatrix> createOperationHestonW(
    base::Grid& grid, base::DataMatrix& coef) {
  if (grid.getType() == base::GridType::Linear) {
    return std::unique_ptr<base::OperationMatrix>(
        new finance::OperationHestonWLinear(&grid.getStorage(), coef));
  } else if (grid.getType() == base::GridType::LinearL0Boundary ||
             grid.getType() == base::GridType::LinearBoundary) {
    return std::unique_ptr<base::OperationMatrix>(
        new finance::OperationHestonWLinearBoundary(&grid.getStorage(), coef));
  } else {
    throw base::factory_exception("OperationHestonW is not implemented for this grid type.");
  }
}

std::unique_ptr<base::OperationMatrix> createOperationHestonZ(
    base::Grid& grid, base::DataVector& coef) {
  if (grid.getType() == base::GridType::Linear) {
    return std::unique_ptr<base::OperationMatrix>(
        new finance::OperationHestonZLinear(&grid.getStorage(), coef));
  } else if (grid.getType() == base::GridType::LinearL0Boundary ||
             grid.getType() == base::GridType::LinearBoundary) {
    return std::unique_ptr<base::OperationMatrix>(
        new finance::OperationHestonZLinearBoundary(&grid.getStorage(), coef));
  } else {
    throw base::factory_exception("OperationHestonZ is not implemented for this grid type.");
  }
}
}  // namespace op_factory
}  // namespace sgpp
