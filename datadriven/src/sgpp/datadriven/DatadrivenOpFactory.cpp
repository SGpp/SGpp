// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <sgpp/datadriven/DatadrivenOpFactory.hpp>

#include <cstring>

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

#ifdef __AVX__
//this use require either SSE3 or AVX
#include <sgpp/datadriven/operation/hash/OperationMultiEvalStreaming/OperationMultiEvalStreaming.hpp>
#include <sgpp/datadriven/operation/hash/OperationMultiEvalModMaskStreaming/OperationMultiEvalModMaskStreaming.hpp>
#include <sgpp/datadriven/operation/hash/OperationMultipleEvalSubspace/combined/OperationMultipleEvalSubspaceCombined.hpp>
#include <sgpp/datadriven/operation/hash/OperationMultipleEvalSubspace/simple/OperationMultipleEvalSubspaceSimple.hpp>
#endif

#ifdef USE_OCL
#include <sgpp/datadriven/operation/hash/OperationMultipleEvalStreamingOCL/StreamingOCLOperatorFactory.hpp>
#include <sgpp/datadriven/operation/hash/OperationMultipleEvalStreamingOCLMultiPlatform/StreamingOCLMultiPlatformOperatorFactory.hpp>
#include "operation/hash/OperationMultipleEvalStreamingModOCL/StreamingModOCLOperatorFactory.hpp"
#include "operation/hash/OperationMultipleEvalStreamingModOCLMask/StreamingModOCLMaskOperatorFactory.hpp"
#include "operation/hash/OperationMultipleEvalStreamingModOCLFast/StreamingModOCLFastOperatorFactory.hpp"
#include "operation/hash/OperationMultipleEvalStreamingModOCLFastMultiPlattform/StreamingModOCLFastMultiPlatformOperatorFactory.hpp"
#include "operation/hash/OperationMultipleEvalAdaptiveOCL/AdaptiveOCLOperatorFactory.hpp"
#include "operation/hash/OperationMultipleEvalStreamingBSplineOCL/StreamingBSplineOCLOperatorFactory.hpp"
#endif

#include <sgpp/base/operation/BaseOpFactory.hpp>

#include <sgpp/globaldef.hpp>

namespace SGPP {

namespace op_factory {

datadriven::OperationTest* createOperationTest(base::Grid& grid) {
    if (grid.getType() == base::GridType::Linear) {
        return new datadriven::OperationTestLinear(grid.getStorage());
    } else if (grid.getType() == base::GridType::LinearL0Boundary || grid.getType() == base::GridType::LinearBoundary) {
        return new datadriven::OperationTestLinearBoundary(grid.getStorage());
    } else if (grid.getType() == base::GridType::ModBspline) {
        return new datadriven::OperationTestModBspline(grid.getStorage(), ((base::ModBsplineGrid*) &grid)->getDegree());
    } else if (grid.getType() == base::GridType::ModLinear) {
        return new datadriven::OperationTestModLinear(grid.getStorage());
    } else if (grid.getType() == base::GridType::Poly) {
        return new datadriven::OperationTestPoly(grid.getStorage(), ((base::PolyGrid*) &grid)->getDegree());
    } else if (grid.getType() == base::GridType::ModPoly) {
        return new datadriven::OperationTestModPoly(grid.getStorage(), ((base::ModPolyGrid*) &grid)->getDegree());
    } else if (grid.getType() == base::GridType::ModWavelet) {
        return new datadriven::OperationTestModWavelet(grid.getStorage());
    } else if (grid.getType() == base::GridType::Prewavelet) {
        return new datadriven::OperationTestPrewavelet(grid.getStorage());
    } else if (grid.getType() == base::GridType::LinearStretched) {
        return new datadriven::OperationTestLinearStretched(grid.getStorage());
    } else if (grid.getType() == base::GridType::LinearStretchedBoundary) {
        return new datadriven::OperationTestLinearStretchedBoundary(grid.getStorage());
    }

    else
        throw base::factory_exception("OperationTest is not implemented for this grid type.");
}

base::OperationMatrix* createOperationRegularizationDiagonal(base::Grid& grid, int mode, float_t k) {
    if (grid.getType() == base::GridType::Linear || grid.getType() == base::GridType::LinearL0Boundary
            || grid.getType() == base::GridType::LinearBoundary || grid.getType() == base::GridType::ModLinear) {
        return new datadriven::OperationRegularizationDiagonalLinearBoundary(grid.getStorage(), mode, k);
    } else
        throw base::factory_exception("OperationRegularizationDiagonal is not implemented for this grid type.");
}

datadriven::OperationDensityMarginalize* createOperationDensityMarginalize(base::Grid& grid) {
    if (grid.getType() == base::GridType::Linear)
        return new datadriven::OperationDensityMarginalizeLinear(&grid);
    else
        throw base::factory_exception("OperationDensityMarginalize is not implemented for this grid type.");
}

datadriven::OperationDensityMargTo1D* createOperationDensityMargTo1D(base::Grid& grid) {
    if (grid.getType() == base::GridType::Linear)
        return new datadriven::OperationDensityMargTo1DLinear(&grid);
    else
        throw base::factory_exception("OperationDensityMargTo1D is not implemented for this grid type.");
}

datadriven::OperationDensitySampling1D* createOperationDensitySampling1D(base::Grid& grid) {
    if (grid.getType() == base::GridType::Linear)
        return new datadriven::OperationDensitySampling1DLinear(&grid);
    else
        throw base::factory_exception("OperationDensitySampling1D is not implemented for this grid type.");
}

datadriven::OperationDensitySampling* createOperationDensitySampling(base::Grid& grid) {
    if (grid.getType() == base::GridType::Linear)
        return new datadriven::OperationDensitySamplingLinear(&grid);
    else
        throw base::factory_exception("OperationDensitySampling is not implemented for this grid type.");
}

datadriven::OperationDensityRejectionSampling* createOperationDensityRejectionSampling(base::Grid& grid) {
    if (grid.getType() == base::GridType::Linear)
        return new datadriven::OperationDensityRejectionSamplingLinear(&grid);
    else
        throw base::factory_exception("OperationDensityRejectionSampling is not implemented for this grid type.");
}

datadriven::OperationDensityConditional* createOperationDensityConditional(base::Grid& grid) {
    if (grid.getType() == base::GridType::Linear)
        return new datadriven::OperationDensityConditionalLinear(&grid);
    else
        throw base::factory_exception("OperationDensityConditional is not implemented for this grid type.");
}

datadriven::OperationRosenblattTransformation* createOperationRosenblattTransformation(base::Grid& grid) {
    if (grid.getType() == base::GridType::Linear)
        return new datadriven::OperationRosenblattTransformationLinear(&grid);
    else
        throw base::factory_exception("OperationRosenblattTransformation is not implemented for this grid type.");
}

datadriven::OperationTransformation1D* createOperationRosenblattTransformation1D(base::Grid& grid) {
    if (grid.getType() == base::GridType::Linear)
        return new datadriven::OperationRosenblattTransformation1DLinear(&grid);
    else
        throw base::factory_exception("OperationRosenblattTransformation1D is not implemented for this grid type.");
}

datadriven::OperationInverseRosenblattTransformation* createOperationInverseRosenblattTransformation(base::Grid& grid) {
    if (grid.getType() == base::GridType::Linear)
        return new datadriven::OperationInverseRosenblattTransformationLinear(&grid);
    else
        throw base::factory_exception(
                "OperationInverseRosenblattTransformation is not implemented for this grid type.");
}

datadriven::OperationTransformation1D* createOperationInverseRosenblattTransformation1D(base::Grid& grid) {
    if (grid.getType() == base::GridType::Linear)
        return new datadriven::OperationInverseRosenblattTransformation1DLinear(&grid);
    else
        throw base::factory_exception(
                "OperationInverseRosenblattTransformation1D is not implemented for this grid type.");
}

datadriven::OperationInverseRosenblattTransformationKDE* createOperationInverseRosenblattTransformationKDE(
        datadriven::GaussianKDE& kde) {
    return new datadriven::OperationInverseRosenblattTransformationKDE(kde);
}

datadriven::OperationRosenblattTransformationKDE* createOperationRosenblattTransformationKDE(
        datadriven::GaussianKDE& kde) {
    return new datadriven::OperationRosenblattTransformationKDE(kde);
}

datadriven::OperationDensityMarginalizeKDE* createOperationDensityMarginalizeKDE(datadriven::GaussianKDE& kde) {
    return new datadriven::OperationDensityMarginalizeKDE(kde);
}

datadriven::OperationDensityConditionalKDE* createOperationDensityConditionalKDE(datadriven::GaussianKDE& kde) {
    return new datadriven::OperationDensityConditionalKDE(kde);
}

base::OperationMultipleEval* createOperationMultipleEval(base::Grid& grid, base::DataMatrix& dataset,
SGPP::datadriven::OperationMultipleEvalConfiguration &configuration) {
    if (configuration.getType() == SGPP::datadriven::OperationMultipleEvalType::DEFAULT) {
        return createOperationMultipleEval(grid, dataset);
    }

    if (grid.getType() == base::GridType::Linear) {
        if (configuration.getType() == datadriven::OperationMultipleEvalType::DEFAULT
                || configuration.getType() == datadriven::OperationMultipleEvalType::STREAMING) {
            if (configuration.getSubType() == SGPP::datadriven::OperationMultipleEvalSubType::DEFAULT) {
#ifdef __AVX__
                return new datadriven::OperationMultiEvalStreaming(grid, dataset);
#else
                throw base::factory_exception("Error creating function: the library wasn't compiled with AVX1/2/512");
#endif
            } else if (configuration.getSubType() == SGPP::datadriven::OperationMultipleEvalSubType::OCL) {
#ifdef USE_OCL
                return datadriven::createStreamingOCLConfigured(grid, dataset, configuration);
#else
                throw base::factory_exception(
                        "Error creating function: the library wasn't compiled with OpenCL support");
#endif
            } else if (configuration.getSubType() == SGPP::datadriven::OperationMultipleEvalSubType::OCLMP) {
#ifdef USE_OCL
                return datadriven::createStreamingOCLMultiPlatformConfigured(grid, dataset, configuration);
#else
                throw base::factory_exception(
                        "Error creating function: the library wasn't compiled with OpenCL support");
#endif
            }
        } else if (configuration.getType() == datadriven::OperationMultipleEvalType::SUBSPACELINEAR) {
            if (configuration.getSubType() == SGPP::datadriven::OperationMultipleEvalSubType::DEFAULT
                    || configuration.getSubType() == SGPP::datadriven::OperationMultipleEvalSubType::COMBINED) {
#ifdef __AVX__
                return new datadriven::OperationMultipleEvalSubspaceCombined(grid, dataset);
#else
                throw base::factory_exception("Error creating function: the library wasn't compiled with AVX");
#endif
            } else if (configuration.getSubType() == SGPP::datadriven::OperationMultipleEvalSubType::SIMPLE) {
#ifdef __AVX__
                return new datadriven::OperationMultipleEvalSubspaceSimple(grid, dataset);
#else
                throw base::factory_exception("Error creating function: the library wasn't compiled with AVX");
#endif
            }
        } else if (configuration.getType() == datadriven::OperationMultipleEvalType::ADAPTIVE) {
            if (configuration.getSubType() == SGPP::datadriven::OperationMultipleEvalSubType::DEFAULT) {
#ifdef USE_OCL
                return datadriven::createAdaptiveOCLConfigured(grid, dataset, configuration);
#else
                throw base::factory_exception(
                        "Error creating function: the library wasn't compiled with OpenCL support");
#endif
            }
        }
    } else if (grid.getType() == base::GridType::ModLinear) {
        if (configuration.getType() == datadriven::OperationMultipleEvalType::STREAMING) {
            if (configuration.getSubType() == SGPP::datadriven::OperationMultipleEvalSubType::DEFAULT) {
#ifdef __AVX__
                return new datadriven::OperationMultiEvalModMaskStreaming(grid, dataset);
#else
                throw base::factory_exception("Error creating function: the library wasn't compiled with AVX");
#endif
            } else if (configuration.getSubType() == SGPP::datadriven::OperationMultipleEvalSubType::OCL) {
#ifdef USE_OCL
                return datadriven::createStreamingModOCLConfigured(grid, dataset, configuration);
                throw base::factory_exception(
                        "Error creating function: the library wasn't compiled with OpenCL support");
#endif
            } else if (configuration.getSubType() == SGPP::datadriven::OperationMultipleEvalSubType::OCLFAST) {
#ifdef USE_OCL
                return datadriven::createStreamingModOCLFastConfigured(grid, dataset, configuration);
#else
                throw base::factory_exception(
                        "Error creating function: the library wasn't compiled with OpenCL support");
#endif
            } else if (configuration.getSubType()
                    == SGPP::datadriven::OperationMultipleEvalSubType::OCLFASTMULTIPLATFORM) {
#ifdef USE_OCL
                return datadriven::createStreamingModOCLFastMultiPlatformConfigured(grid, dataset, configuration);
#else
                throw base::factory_exception(
                        "Error creating function: the library wasn't compiled with OpenCL support");
#endif
            } else if (configuration.getSubType() == SGPP::datadriven::OperationMultipleEvalSubType::OCLMASK) {
#ifdef USE_OCL
                return datadriven::createStreamingModOCLMaskConfigured(grid, dataset, configuration);
#else
                throw base::factory_exception(
                        "Error creating function: the library wasn't compiled with OpenCL support");
#endif
            }
        }
    } else if (grid.getType() == base::GridType::Bspline) {
        if (configuration.getType() == datadriven::OperationMultipleEvalType::STREAMING) {
            if (configuration.getSubType() == SGPP::datadriven::OperationMultipleEvalSubType::OCL) {
#ifdef USE_OCL
                return datadriven::createStreamingBSplineOCLConfigured(grid, dataset, configuration);
#else
                throw base::factory_exception(
                        "Error creating function: the library wasn't compiled with OpenCL support");
#endif
            }
        }
    }
    throw base::factory_exception("OperationMultiEval is not implemented for this grid type.");
}

}
}
