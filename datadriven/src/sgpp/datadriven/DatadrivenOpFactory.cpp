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

#include <sgpp/datadriven/operation/hash/OperationTestLinear.hpp>
#include <sgpp/datadriven/operation/hash/OperationTestLinearBoundary.hpp>
#include <sgpp/datadriven/operation/hash/OperationTestModBspline.hpp>
#include <sgpp/datadriven/operation/hash/OperationTestModLinear.hpp>
#include <sgpp/datadriven/operation/hash/OperationTestPoly.hpp>
#include <sgpp/datadriven/operation/hash/OperationTestModPoly.hpp>
#include <sgpp/datadriven/operation/hash/OperationTestModWavelet.hpp>
#include <sgpp/datadriven/operation/hash/OperationTestPrewavelet.hpp>
#include <sgpp/datadriven/operation/hash/OperationTestLinearStretchedBoundary.hpp>
#include <sgpp/datadriven/operation/hash/OperationTestLinearStretched.hpp>
#include <sgpp/datadriven/operation/hash/OperationDensityMarginalizeLinear.hpp>
#include <sgpp/datadriven/operation/hash/OperationDensityMargTo1DLinear.hpp>
#include <sgpp/datadriven/operation/hash/OperationDensitySampling1DLinear.hpp>
#include <sgpp/datadriven/operation/hash/OperationDensitySamplingLinear.hpp>
#include <sgpp/datadriven/operation/hash/OperationDensityRejectionSamplingLinear.hpp>
#include <sgpp/datadriven/operation/hash/OperationDensityConditionalLinear.hpp>
#include <sgpp/datadriven/operation/hash/OperationRosenblattTransformationLinear.hpp>
#include <sgpp/datadriven/operation/hash/OperationRosenblattTransformation1DLinear.hpp>
#include <sgpp/datadriven/operation/hash/OperationInverseRosenblattTransformationLinear.hpp>
#include <sgpp/datadriven/operation/hash/OperationRegularizationDiagonalLinearBoundary.hpp>

#ifdef __AVX__
#include <sgpp/datadriven/operation/hash/OperationMultiEvalStreaming/OperationMultiEvalStreaming.hpp>
#include <sgpp/datadriven/operation/hash/OperationMultipleEvalSubspace/combined/OperationMultipleEvalSubspaceCombined.hpp>
#include <sgpp/datadriven/operation/hash/OperationMultipleEvalSubspace/simple/OperationMultipleEvalSubspaceSimple.hpp>
#endif

#ifdef USE_OCL
#include <sgpp/datadriven/operation/hash/OperationMultipleEvalStreamingOCL/StreamingOCLOperatorFactory.hpp>
#endif

#include <sgpp/base/operation/BaseOpFactory.hpp>

#include <sgpp/globaldef.hpp>

namespace SGPP {

  namespace op_factory {

    datadriven::OperationTest* createOperationTest(base::Grid& grid) {
      if (strcmp(grid.getType(), "linear") == 0) {
        return new datadriven::OperationTestLinear(grid.getStorage());
      } else if (strcmp(grid.getType(), "linearBoundary") == 0
                 || strcmp(grid.getType(), "linearTruncatedBoundary") == 0) {
        return new datadriven::OperationTestLinearBoundary(grid.getStorage());
      } else if (strcmp(grid.getType(), "modBspline") == 0) {
        return new datadriven::OperationTestModBspline(grid.getStorage(),
               ((base::ModBsplineGrid*) &grid)->getDegree());
      } else if (strcmp(grid.getType(), "modlinear") == 0) {
        return new datadriven::OperationTestModLinear(grid.getStorage());
      } else if (strcmp(grid.getType(), "poly") == 0) {
        return new datadriven::OperationTestPoly(grid.getStorage(),
               ((base::PolyGrid*) &grid)->getDegree());
      } else if (strcmp(grid.getType(), "modpoly") == 0) {
        return new datadriven::OperationTestModPoly(grid.getStorage(),
               ((base::ModPolyGrid*) &grid)->getDegree());
      } else if (strcmp(grid.getType(), "modWavelet") == 0) {
        return new datadriven::OperationTestModWavelet(grid.getStorage());
      } else if (strcmp(grid.getType(), "prewavelet") == 0) {
        return new datadriven::OperationTestPrewavelet(grid.getStorage());
      } else if (strcmp(grid.getType(), "linearStretched") == 0) {
        return new datadriven::OperationTestLinearStretched(grid.getStorage());
      } else if (strcmp(grid.getType(), "linearStretchedTruncatedBoundary")
                 == 0) {
        return new datadriven::OperationTestLinearStretchedBoundary(
                 grid.getStorage());
      }

      else
        throw base::factory_exception(
          "OperationTest is not implemented for this grid type.");
    }

    base::OperationMatrix* createOperationRegularizationDiagonal(base::Grid& grid,
        int mode, float_t k) {
      if (strcmp(grid.getType(), "linear") == 0
          || strcmp(grid.getType(), "linearBoundary") == 0
          || strcmp(grid.getType(), "linearTruncatedBoundary") == 0
          || strcmp(grid.getType(), "modlinear") == 0) {
        return new datadriven::OperationRegularizationDiagonalLinearBoundary(
                 grid.getStorage(), mode, k);
      } else
        throw base::factory_exception(
          "OperationRegularizationDiagonal is not implemented for this grid type.");
    }

    datadriven::OperationDensityMarginalize* createOperationDensityMarginalize(
      base::Grid& grid) {
      if (strcmp(grid.getType(), "linear") == 0)
        return new datadriven::OperationDensityMarginalizeLinear(&grid);
      else
        throw base::factory_exception(
          "OperationDensityMarginalize is not implemented for this grid type.");
    }

    datadriven::OperationDensityMargTo1D* createOperationDensityMargTo1D(
      base::Grid& grid) {
      if (strcmp(grid.getType(), "linear") == 0)
        return new datadriven::OperationDensityMargTo1DLinear(&grid);
      else
        throw base::factory_exception(
          "OperationDensityMargTo1D is not implemented for this grid type.");
    }

    datadriven::OperationDensitySampling1D* createOperationDensitySampling1D(
      base::Grid& grid) {
      if (strcmp(grid.getType(), "linear") == 0)
        return new datadriven::OperationDensitySampling1DLinear(&grid);
      else
        throw base::factory_exception(
          "OperationDensitySampling1D is not implemented for this grid type.");
    }

    datadriven::OperationDensitySampling* createOperationDensitySampling(
      base::Grid& grid) {
      if (strcmp(grid.getType(), "linear") == 0)
        return new datadriven::OperationDensitySamplingLinear(&grid);
      else
        throw base::factory_exception(
          "OperationDensitySampling is not implemented for this grid type.");
    }

    datadriven::OperationDensityRejectionSampling* createOperationDensityRejectionSampling(
      base::Grid& grid) {
      if (strcmp(grid.getType(), "linear") == 0)
        return new datadriven::OperationDensityRejectionSamplingLinear(&grid);
      else
        throw base::factory_exception(
          "OperationDensityRejectionSampling is not implemented for this grid type.");
    }

    datadriven::OperationDensityConditional* createOperationDensityConditional(
      base::Grid& grid) {
      if (strcmp(grid.getType(), "linear") == 0)
        return new datadriven::OperationDensityConditionalLinear(&grid);
      else
        throw base::factory_exception(
          "OperationDensityConditional is not implemented for this grid type.");
    }

    datadriven::OperationRosenblattTransformation* createOperationRosenblattTransformation(
      base::Grid& grid) {
      if (strcmp(grid.getType(), "linear") == 0)
        return new datadriven::OperationRosenblattTransformationLinear(&grid);
      else
        throw base::factory_exception(
          "OperationRosenblattTransformation is not implemented for this grid type.");
    }

    datadriven::OperationTransformation1D* createOperationRosenblattTransformation1D(
      base::Grid& grid) {
      if (strcmp(grid.getType(), "linear") == 0)
        return new datadriven::OperationRosenblattTransformation1DLinear(&grid);
      else
        throw base::factory_exception(
          "OperationRosenblattTransformation1D is not implemented for this grid type.");
    }

    datadriven::OperationInverseRosenblattTransformation* createOperationInverseRosenblattTransformation(
      base::Grid& grid) {
      if (strcmp(grid.getType(), "linear") == 0)
        return new datadriven::OperationInverseRosenblattTransformationLinear(
                 &grid);
      else
        throw base::factory_exception(
          "OperationInverseRosenblattTransformation is not implemented for this grid type.");
    }

    base::OperationMultipleEval* createOperationMultipleEval(base::Grid& grid,
        base::DataMatrix& dataset,
        SGPP::datadriven::OperationMultipleEvalConfiguration configuration) {
      if (configuration.type
          == SGPP::datadriven::OperationMultipleEvalType::DEFAULT) {
        return createOperationMultipleEval(grid, dataset);
      }

      if (strcmp(grid.getType(), "linear") == 0) {
        switch (configuration.type) {
          case datadriven::OperationMultipleEvalType::DEFAULT:

          //already handled
          case datadriven::OperationMultipleEvalType::STREAMING:
            if (configuration.subType
                == SGPP::datadriven::OperationMultipleEvalSubType::DEFAULT) {
#ifdef __AVX__
              return new datadriven::OperationMultiEvalStreaming(grid,
                     dataset);
#else
              throw base::factory_exception("Error creating function: the library wasn't compiled with AVX");
#endif
              break;
            } else if (configuration.subType
                       == SGPP::datadriven::OperationMultipleEvalSubType::OCL) {
#ifdef USE_OCL
              return datadriven::createStreamingOCLConfigured(grid, dataset);
              //        return new datadriven::OperationMultiEvalStreamingOCL<double>(grid,
              //            dataset);
#else
              throw base::factory_exception("Error creating function: the library wasn't compiled with OpenCL support");
#endif
            }

            break;

          case datadriven::OperationMultipleEvalType::SUBSPACELINEAR:
            switch (configuration.subType) {
              case SGPP::datadriven::OperationMultipleEvalSubType::DEFAULT:
              case SGPP::datadriven::OperationMultipleEvalSubType::COMBINED:
#ifdef __AVX__
                return new datadriven::OperationMultipleEvalSubspaceCombined(
                         grid, dataset);
#else
                throw base::factory_exception("Error creating function: the library wasn't compiled with AVX");
#endif
                break;

              case SGPP::datadriven::OperationMultipleEvalSubType::SIMPLE:
#ifdef __AVX__
                return new datadriven::OperationMultipleEvalSubspaceSimple(grid,
                       dataset);
#else
                throw base::factory_exception("Error creating function: the library wasn't compiled with AVX");
#endif
                break;

              case SGPP::datadriven::OperationMultipleEvalSubType::OCL:
              default:
                break;
            }

            break;

          default:
            break;
        }
      }

      throw base::factory_exception(
        "OperationMultiEval is not implemented for this grid type.");
    }

  }
}
