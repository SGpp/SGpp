// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <sgpp/datadriven/DatadrivenOpFactory.hpp>

#include <sgpp/base/exception/factory_exception.hpp>

#include <sgpp/base/grid/type/PolyGrid.hpp>
#include <sgpp/base/grid/type/ModPolyGrid.hpp>
#include <sgpp/base/grid/type/PrewaveletGrid.hpp>
#include <sgpp/base/grid/type/ModBsplineGrid.hpp>

#include <sgpp/datadriven/operation/hash/simple/OperationTestLinear.hpp>
#include <sgpp/datadriven/operation/hash/simple/OperationTestLinearBoundary.hpp>
#include <sgpp/datadriven/operation/hash/simple/OperationTestModBspline.hpp>
#include <sgpp/datadriven/operation/hash/simple/OperationTestModLinear.hpp>
#include <sgpp/datadriven/operation/hash/simple/OperationTestPoly.hpp>
#include <sgpp/datadriven/operation/hash/simple/OperationTestModPoly.hpp>
#include <sgpp/datadriven/operation/hash/simple/OperationTestModWavelet.hpp>
#include <sgpp/datadriven/operation/hash/simple/OperationTestPrewavelet.hpp>
#include <sgpp/datadriven/operation/hash/simple/OperationTestLinearStretchedBoundary.hpp>
#include <sgpp/datadriven/operation/hash/simple/OperationTestLinearStretched.hpp>
#include <sgpp/datadriven/operation/hash/simple/OperationDensityMarginalizeLinear.hpp>
#include <sgpp/datadriven/operation/hash/simple/OperationDensityMargTo1DLinear.hpp>
#include <sgpp/datadriven/operation/hash/simple/OperationDensitySampling1DLinear.hpp>
#include <sgpp/datadriven/operation/hash/simple/OperationDensitySamplingLinear.hpp>
#include <sgpp/datadriven/operation/hash/simple/OperationDensityRejectionSamplingLinear.hpp>
#include <sgpp/datadriven/operation/hash/simple/OperationDensityConditionalLinear.hpp>
#include <sgpp/datadriven/operation/hash/simple/OperationRosenblattTransformationLinear.hpp>
#include <sgpp/datadriven/operation/hash/simple/OperationRosenblattTransformation1DLinear.hpp>
#include <sgpp/datadriven/operation/hash/simple/OperationInverseRosenblattTransformationLinear.hpp>
#include <sgpp/datadriven/operation/hash/simple/OperationInverseRosenblattTransformation1DLinear.hpp>
#include <sgpp/datadriven/operation/hash/simple/OperationRegularizationDiagonalLinearBoundary.hpp>

#include <sgpp/datadriven/operation/hash/OperationMultiEvalModMaskStreaming/OperationMultiEvalModMaskStreaming.hpp>
#include <sgpp/datadriven/operation/hash/OperationMultiEvalStreaming/OperationMultiEvalStreaming.hpp>

#ifdef __AVX__
#include <sgpp/datadriven/operation/hash/OperationMultipleEvalSubspace/combined/OperationMultipleEvalSubspaceCombined.hpp>
#include <sgpp/datadriven/operation/hash/OperationMultipleEvalSubspace/simple/OperationMultipleEvalSubspaceSimple.hpp>
#endif

#ifdef USE_OCL
#include "operation/hash/OperationMultipleEvalStreamingOCLMultiPlatform/OperatorFactory.hpp"
#include "operation/hash/OperationMultipleEvalStreamingModOCLMaskMultiPlatform/OperatorFactory.hpp"
#include "operation/hash/OperationMultipleEvalStreamingModOCLFastMultiPlattform/OperatorFactory.hpp"
#include "operation/hash/OperationMultipleEvalStreamingModOCLOpt/OperatorFactory.hpp"
#include "operation/hash/OperationMultipleEvalAdaptiveOCL/AdaptiveOCLOperatorFactory.hpp"
#include "operation/hash/OperationMultipleEvalStreamingBSplineOCL/StreamingBSplineOCLOperatorFactory.hpp"
#endif

#include <sgpp/base/operation/BaseOpFactory.hpp>
#include <sgpp/globaldef.hpp>

#include <cstring>

namespace sgpp {
namespace op_factory {

std::unique_ptr<datadriven::OperationTest> createOperationTest(base::Grid& grid) {
  if (grid.getType() == base::GridType::Linear) {
    return std::unique_ptr<datadriven::OperationTest>(
        new datadriven::OperationTestLinear(&grid.getStorage()));
  } else if (grid.getType() == base::GridType::LinearL0Boundary ||
             grid.getType() == base::GridType::LinearBoundary) {
    return std::unique_ptr<datadriven::OperationTest>(
        new datadriven::OperationTestLinearBoundary(&grid.getStorage()));
  } else if (grid.getType() == base::GridType::ModBspline) {
    return std::unique_ptr<datadriven::OperationTest>(new datadriven::OperationTestModBspline(
        &grid.getStorage(), ((base::ModBsplineGrid*)&grid)->getDegree()));
  } else if (grid.getType() == base::GridType::ModLinear) {
    return std::unique_ptr<datadriven::OperationTest>(
        new datadriven::OperationTestModLinear(&grid.getStorage()));
  } else if (grid.getType() == base::GridType::Poly) {
    return std::unique_ptr<datadriven::OperationTest>(new datadriven::OperationTestPoly(
        &grid.getStorage(), ((base::PolyGrid*)&grid)->getDegree()));
  } else if (grid.getType() == base::GridType::ModPoly) {
    return std::unique_ptr<datadriven::OperationTest>(new datadriven::OperationTestModPoly(
        &grid.getStorage(), ((base::ModPolyGrid*)&grid)->getDegree()));
  } else if (grid.getType() == base::GridType::ModWavelet) {
    return std::unique_ptr<datadriven::OperationTest>(
        new datadriven::OperationTestModWavelet(&grid.getStorage()));
  } else if (grid.getType() == base::GridType::Prewavelet) {
    return std::unique_ptr<datadriven::OperationTest>(
        new datadriven::OperationTestPrewavelet(&grid.getStorage()));
  } else if (grid.getType() == base::GridType::LinearStretched) {
    return std::unique_ptr<datadriven::OperationTest>(
        new datadriven::OperationTestLinearStretched(&grid.getStorage()));
  } else if (grid.getType() == base::GridType::LinearStretchedBoundary) {
    return std::unique_ptr<datadriven::OperationTest>(
        new datadriven::OperationTestLinearStretchedBoundary(&grid.getStorage()));
  } else {
    throw base::factory_exception("OperationTest is not implemented for this grid type.");
  }
}

std::unique_ptr<base::OperationMatrix> createOperationRegularizationDiagonal(base::Grid& grid,
                                                                             int mode, double k) {
  if (grid.getType() == base::GridType::Linear ||
      grid.getType() == base::GridType::LinearL0Boundary ||
      grid.getType() == base::GridType::LinearBoundary ||
      grid.getType() == base::GridType::ModLinear) {
    return std::unique_ptr<base::OperationMatrix>(
        new datadriven::OperationRegularizationDiagonalLinearBoundary(&grid.getStorage(), mode, k));
  } else {
    throw base::factory_exception(
        "OperationRegularizationDiagonal is not implemented for this grid type.");
  }
}

std::unique_ptr<datadriven::OperationDensityMarginalize> createOperationDensityMarginalize(
    base::Grid& grid) {
  if (grid.getType() == base::GridType::Linear)
    return std::unique_ptr<datadriven::OperationDensityMarginalize>(
        new datadriven::OperationDensityMarginalizeLinear(&grid));
  else
    throw base::factory_exception(
        "OperationDensityMarginalize is not implemented for this grid type.");
}

std::unique_ptr<datadriven::OperationDensityMargTo1D> createOperationDensityMargTo1D(
    base::Grid& grid) {
  if (grid.getType() == base::GridType::Linear)
    return std::unique_ptr<datadriven::OperationDensityMargTo1D>(
        new datadriven::OperationDensityMargTo1DLinear(&grid));
  else
    throw base::factory_exception(
        "OperationDensityMargTo1D is not implemented for this grid type.");
}

std::unique_ptr<datadriven::OperationDensitySampling1D> createOperationDensitySampling1D(
    base::Grid& grid) {
  if (grid.getType() == base::GridType::Linear)
    return std::unique_ptr<datadriven::OperationDensitySampling1D>(
        new datadriven::OperationDensitySampling1DLinear(&grid));
  else
    throw base::factory_exception(
        "OperationDensitySampling1D is not implemented for this grid type.");
}

std::unique_ptr<datadriven::OperationDensitySampling> createOperationDensitySampling(
    base::Grid& grid) {
  if (grid.getType() == base::GridType::Linear)
    return std::unique_ptr<datadriven::OperationDensitySampling>(
        new datadriven::OperationDensitySamplingLinear(&grid));
  else
    throw base::factory_exception(
        "OperationDensitySampling is not implemented for this grid type.");
}

std::unique_ptr<datadriven::OperationDensityRejectionSampling>
createOperationDensityRejectionSampling(base::Grid& grid) {
  if (grid.getType() == base::GridType::Linear)
    return std::unique_ptr<datadriven::OperationDensityRejectionSampling>(
        new datadriven::OperationDensityRejectionSamplingLinear(&grid));
  else
    throw base::factory_exception(
        "OperationDensityRejectionSampling is not implemented for this grid type.");
}

std::unique_ptr<datadriven::OperationDensityConditional> createOperationDensityConditional(
    base::Grid& grid) {
  if (grid.getType() == base::GridType::Linear)
    return std::unique_ptr<datadriven::OperationDensityConditional>(
        new datadriven::OperationDensityConditionalLinear(&grid));
  else
    throw base::factory_exception(
        "OperationDensityConditional is not implemented for this grid type.");
}

std::unique_ptr<datadriven::OperationRosenblattTransformation>
createOperationRosenblattTransformation(base::Grid& grid) {
  if (grid.getType() == base::GridType::Linear)
    return std::unique_ptr<datadriven::OperationRosenblattTransformation>(
        new datadriven::OperationRosenblattTransformationLinear(&grid));
  else
    throw base::factory_exception(
        "OperationRosenblattTransformation is not implemented for this grid type.");
}

std::unique_ptr<datadriven::OperationTransformation1D> createOperationRosenblattTransformation1D(
    base::Grid& grid) {
  if (grid.getType() == base::GridType::Linear)
    return std::unique_ptr<datadriven::OperationTransformation1D>(
        new datadriven::OperationRosenblattTransformation1DLinear(&grid));
  else
    throw base::factory_exception(
        "OperationRosenblattTransformation1D is not implemented for this grid type.");
}

std::unique_ptr<datadriven::OperationInverseRosenblattTransformation>
createOperationInverseRosenblattTransformation(base::Grid& grid) {
  if (grid.getType() == base::GridType::Linear)
    return std::unique_ptr<datadriven::OperationInverseRosenblattTransformation>(
        new datadriven::OperationInverseRosenblattTransformationLinear(&grid));
  else
    throw base::factory_exception(
        "OperationInverseRosenblattTransformation is not implemented for this grid type.");
}

std::unique_ptr<datadriven::OperationTransformation1D>
createOperationInverseRosenblattTransformation1D(base::Grid& grid) {
  if (grid.getType() == base::GridType::Linear)
    return std::unique_ptr<datadriven::OperationTransformation1D>(
        new datadriven::OperationInverseRosenblattTransformation1DLinear(&grid));
  else
    throw base::factory_exception(
        "OperationInverseRosenblattTransformation1D is not implemented for this grid type.");
}

std::unique_ptr<datadriven::OperationInverseRosenblattTransformationKDE>
createOperationInverseRosenblattTransformationKDE(datadriven::KernelDensityEstimator& kde) {
  return std::unique_ptr<datadriven::OperationInverseRosenblattTransformationKDE>(
      new datadriven::OperationInverseRosenblattTransformationKDE(kde));
}

std::unique_ptr<datadriven::OperationRosenblattTransformationKDE>
createOperationRosenblattTransformationKDE(datadriven::KernelDensityEstimator& kde) {
  return std::unique_ptr<datadriven::OperationRosenblattTransformationKDE>(
      new datadriven::OperationRosenblattTransformationKDE(kde));
}

std::unique_ptr<datadriven::OperationDensityMarginalizeKDE> createOperationDensityMarginalizeKDE(
    datadriven::KernelDensityEstimator& kde) {
  return std::unique_ptr<datadriven::OperationDensityMarginalizeKDE>(
      new datadriven::OperationDensityMarginalizeKDE(kde));
}

std::unique_ptr<datadriven::OperationDensityConditionalKDE> createOperationDensityConditionalKDE(
    datadriven::KernelDensityEstimator& kde) {
  return std::unique_ptr<datadriven::OperationDensityConditionalKDE>(
      new datadriven::OperationDensityConditionalKDE(kde));
}

std::unique_ptr<base::OperationMultipleEval> createOperationMultipleEval(
    base::Grid& grid, base::DataMatrix& dataset,
    sgpp::datadriven::OperationMultipleEvalConfiguration& configuration) {
  if (configuration.getType() == sgpp::datadriven::OperationMultipleEvalType::DEFAULT) {
    return createOperationMultipleEval(grid, dataset);
  }

  if (grid.getType() == base::GridType::Linear) {
    if (configuration.getType() == datadriven::OperationMultipleEvalType::DEFAULT ||
        configuration.getType() == datadriven::OperationMultipleEvalType::STREAMING) {
      if (configuration.getSubType() == sgpp::datadriven::OperationMultipleEvalSubType::DEFAULT) {
        return std::unique_ptr<base::OperationMultipleEval>(
            new datadriven::OperationMultiEvalStreaming(grid, dataset));
      }
      if (configuration.getSubType() == sgpp::datadriven::OperationMultipleEvalSubType::OCLMP) {
#ifdef USE_OCL
        return std::unique_ptr<base::OperationMultipleEval>(
            datadriven::createStreamingOCLMultiPlatformConfigured(grid, dataset, configuration));
#else
        throw base::factory_exception(
            "Error creating function: the library wasn't compiled with OpenCL support");
#endif
      }
    } else if (configuration.getType() == datadriven::OperationMultipleEvalType::SUBSPACELINEAR) {
      if (configuration.getSubType() == sgpp::datadriven::OperationMultipleEvalSubType::DEFAULT ||
          configuration.getSubType() == sgpp::datadriven::OperationMultipleEvalSubType::COMBINED) {
#ifdef __AVX__
        return std::unique_ptr<base::OperationMultipleEval>(
            new datadriven::OperationMultipleEvalSubspaceCombined(grid, dataset));
#else
        throw base::factory_exception(
            "Error creating function: the library wasn't compiled with AVX");
#endif
      } else if (configuration.getSubType() ==
                 sgpp::datadriven::OperationMultipleEvalSubType::SIMPLE) {
#ifdef __AVX__
        return std::unique_ptr<base::OperationMultipleEval>(
            new datadriven::OperationMultipleEvalSubspaceSimple(grid, dataset));
#else
        throw base::factory_exception(
            "Error creating function: the library wasn't compiled with AVX");
#endif
      }
    } else if (configuration.getType() == datadriven::OperationMultipleEvalType::ADAPTIVE) {
      if (configuration.getSubType() == sgpp::datadriven::OperationMultipleEvalSubType::DEFAULT) {
#ifdef USE_OCL
        return std::unique_ptr<base::OperationMultipleEval>(
            datadriven::createAdaptiveOCLConfigured(grid, dataset, configuration));
#else
        throw base::factory_exception(
            "Error creating function: the library wasn't compiled with OpenCL support");
#endif
      }
    }
  } else if (grid.getType() == base::GridType::ModLinear) {
    if (configuration.getType() == datadriven::OperationMultipleEvalType::STREAMING) {
      if (configuration.getSubType() == sgpp::datadriven::OperationMultipleEvalSubType::DEFAULT) {
        return std::unique_ptr<base::OperationMultipleEval>(
            new datadriven::OperationMultiEvalModMaskStreaming(grid, dataset));
      }
      if (configuration.getSubType() == sgpp::datadriven::OperationMultipleEvalSubType::OCLFASTMP) {
#ifdef USE_OCL
        return std::unique_ptr<base::OperationMultipleEval>(
            datadriven::createStreamingModOCLFastMultiPlatformConfigured(grid, dataset,
                                                                         configuration));
#else
        throw base::factory_exception(
            "Error creating function: the library wasn't compiled with OpenCL support");
#endif
      } else if (configuration.getSubType() ==
                 sgpp::datadriven::OperationMultipleEvalSubType::OCLMASKMP) {
#ifdef USE_OCL
        return std::unique_ptr<base::OperationMultipleEval>(
            datadriven::createStreamingModOCLMaskMultiPlatformConfigured(grid, dataset,
                                                                         configuration));
#else
        throw base::factory_exception(
            "Error creating function: the library wasn't compiled with OpenCL support");
#endif
      } else if (configuration.getSubType() ==
                 sgpp::datadriven::OperationMultipleEvalSubType::OCLOPT) {
#ifdef USE_OCL
        return std::unique_ptr<base::OperationMultipleEval>(
            datadriven::createStreamingModOCLOptConfigured(grid, dataset, configuration));
#else
        throw base::factory_exception(
            "Error creating function: the library wasn't compiled with OpenCL support");
#endif
      }
    }
  } else if (grid.getType() == base::GridType::Bspline) {
    if (configuration.getType() == datadriven::OperationMultipleEvalType::STREAMING) {
      if (configuration.getSubType() == sgpp::datadriven::OperationMultipleEvalSubType::OCL) {
#ifdef USE_OCL
        return std::unique_ptr<base::OperationMultipleEval>(
            datadriven::createStreamingBSplineOCLConfigured(grid, dataset, configuration));
#else
        throw base::factory_exception(
            "Error creating function: the library wasn't compiled with OpenCL support");
#endif
      }
    }
  }

  throw base::factory_exception("OperationMultiEval is not implemented for this grid type.");
}
}  // namespace op_factory
}  // namespace sgpp
