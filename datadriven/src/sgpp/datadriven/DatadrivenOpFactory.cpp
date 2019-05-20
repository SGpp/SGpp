// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <sgpp/datadriven/DatadrivenOpFactory.hpp>

#include <sgpp/base/exception/factory_exception.hpp>

#include <sgpp/base/grid/type/ModBsplineGrid.hpp>
#include <sgpp/base/grid/type/ModPolyGrid.hpp>
#include <sgpp/base/grid/type/PolyGrid.hpp>
#include <sgpp/base/grid/type/PrewaveletGrid.hpp>

#include <sgpp/datadriven/operation/hash/simple/OperationDensityConditional.hpp>
#include <sgpp/datadriven/operation/hash/simple/OperationDensityConditionalLinear.hpp>
#include <sgpp/datadriven/operation/hash/simple/OperationDensityMargTo1D.hpp>
#include <sgpp/datadriven/operation/hash/simple/OperationDensityMarginalize.hpp>
#include <sgpp/datadriven/operation/hash/simple/OperationDensityMarginalizeLinear.hpp>
#include <sgpp/datadriven/operation/hash/simple/OperationDensityRejectionSamplingLinear.hpp>
#include <sgpp/datadriven/operation/hash/simple/OperationDensitySampling1DLinear.hpp>
#include <sgpp/datadriven/operation/hash/simple/OperationDensitySamplingLinear.hpp>

#include <sgpp/datadriven/operation/hash/simple/OperationRegularizationDiagonalLinearBoundary.hpp>

#include <sgpp/datadriven/operation/hash/simple/OperationInverseRosenblattTransformation1DBspline.hpp>
#include <sgpp/datadriven/operation/hash/simple/OperationInverseRosenblattTransformation1DBsplineBoundary.hpp>
#include <sgpp/datadriven/operation/hash/simple/OperationInverseRosenblattTransformation1DBsplineClenshawCurtis.hpp>
#include <sgpp/datadriven/operation/hash/simple/OperationInverseRosenblattTransformation1DLinear.hpp>
#include <sgpp/datadriven/operation/hash/simple/OperationInverseRosenblattTransformation1DModBspline.hpp>
#include <sgpp/datadriven/operation/hash/simple/OperationInverseRosenblattTransformation1DModBsplineClenshawCurtis.hpp>
#include <sgpp/datadriven/operation/hash/simple/OperationInverseRosenblattTransformation1DModPoly.hpp>
#include <sgpp/datadriven/operation/hash/simple/OperationInverseRosenblattTransformation1DModPolyClenshawCurtis.hpp>
#include <sgpp/datadriven/operation/hash/simple/OperationInverseRosenblattTransformation1DPoly.hpp>
#include <sgpp/datadriven/operation/hash/simple/OperationInverseRosenblattTransformation1DPolyBoundary.hpp>
#include <sgpp/datadriven/operation/hash/simple/OperationInverseRosenblattTransformation1DPolyClenshawCurtis.hpp>
#include <sgpp/datadriven/operation/hash/simple/OperationInverseRosenblattTransformation1DPolyClenshawCurtisBoundary.hpp>

#include <sgpp/datadriven/operation/hash/simple/OperationInverseRosenblattTransformationBspline.hpp>
#include <sgpp/datadriven/operation/hash/simple/OperationInverseRosenblattTransformationBsplineBoundary.hpp>
#include <sgpp/datadriven/operation/hash/simple/OperationInverseRosenblattTransformationBsplineClenshawCurtis.hpp>
#include <sgpp/datadriven/operation/hash/simple/OperationInverseRosenblattTransformationLinear.hpp>
#include <sgpp/datadriven/operation/hash/simple/OperationInverseRosenblattTransformationModBspline.hpp>
#include <sgpp/datadriven/operation/hash/simple/OperationInverseRosenblattTransformationModBsplineClenshawCurtis.hpp>
#include <sgpp/datadriven/operation/hash/simple/OperationInverseRosenblattTransformationModPoly.hpp>
#include <sgpp/datadriven/operation/hash/simple/OperationInverseRosenblattTransformationModPolyClenshawCurtis.hpp>
#include <sgpp/datadriven/operation/hash/simple/OperationInverseRosenblattTransformationPoly.hpp>
#include <sgpp/datadriven/operation/hash/simple/OperationInverseRosenblattTransformationPolyBoundary.hpp>
#include <sgpp/datadriven/operation/hash/simple/OperationInverseRosenblattTransformationPolyClenshawCurtis.hpp>
#include <sgpp/datadriven/operation/hash/simple/OperationInverseRosenblattTransformationPolyClenshawCurtisBoundary.hpp>

#include <sgpp/datadriven/operation/hash/simple/OperationRosenblattTransformation1DBspline.hpp>
#include <sgpp/datadriven/operation/hash/simple/OperationRosenblattTransformation1DBsplineBoundary.hpp>
#include <sgpp/datadriven/operation/hash/simple/OperationRosenblattTransformation1DBsplineClenshawCurtis.hpp>
#include <sgpp/datadriven/operation/hash/simple/OperationRosenblattTransformation1DLinear.hpp>
#include <sgpp/datadriven/operation/hash/simple/OperationRosenblattTransformation1DModBspline.hpp>
#include <sgpp/datadriven/operation/hash/simple/OperationRosenblattTransformation1DModBsplineClenshawCurtis.hpp>
#include <sgpp/datadriven/operation/hash/simple/OperationRosenblattTransformation1DModPoly.hpp>
#include <sgpp/datadriven/operation/hash/simple/OperationRosenblattTransformation1DModPolyClenshawCurtis.hpp>
#include <sgpp/datadriven/operation/hash/simple/OperationRosenblattTransformation1DPoly.hpp>
#include <sgpp/datadriven/operation/hash/simple/OperationRosenblattTransformation1DPolyBoundary.hpp>
#include <sgpp/datadriven/operation/hash/simple/OperationRosenblattTransformation1DPolyClenshawCurtis.hpp>
#include <sgpp/datadriven/operation/hash/simple/OperationRosenblattTransformation1DPolyClenshawCurtisBoundary.hpp>

#include <sgpp/datadriven/operation/hash/simple/OperationRosenblattTransformationBspline.hpp>
#include <sgpp/datadriven/operation/hash/simple/OperationRosenblattTransformationBsplineBoundary.hpp>
#include <sgpp/datadriven/operation/hash/simple/OperationRosenblattTransformationBsplineClenshawCurtis.hpp>
#include <sgpp/datadriven/operation/hash/simple/OperationRosenblattTransformationLinear.hpp>
#include <sgpp/datadriven/operation/hash/simple/OperationRosenblattTransformationModBspline.hpp>
#include <sgpp/datadriven/operation/hash/simple/OperationRosenblattTransformationModBsplineClenshawCurtis.hpp>
#include <sgpp/datadriven/operation/hash/simple/OperationRosenblattTransformationModPoly.hpp>
#include <sgpp/datadriven/operation/hash/simple/OperationRosenblattTransformationModPolyClenshawCurtis.hpp>
#include <sgpp/datadriven/operation/hash/simple/OperationRosenblattTransformationPoly.hpp>
#include <sgpp/datadriven/operation/hash/simple/OperationRosenblattTransformationPolyBoundary.hpp>
#include <sgpp/datadriven/operation/hash/simple/OperationRosenblattTransformationPolyClenshawCurtis.hpp>
#include <sgpp/datadriven/operation/hash/simple/OperationRosenblattTransformationPolyClenshawCurtisBoundary.hpp>

#include <sgpp/datadriven/operation/hash/simple/OperationTestLinear.hpp>
#include <sgpp/datadriven/operation/hash/simple/OperationTestLinearBoundary.hpp>
#include <sgpp/datadriven/operation/hash/simple/OperationTestLinearStretched.hpp>
#include <sgpp/datadriven/operation/hash/simple/OperationTestLinearStretchedBoundary.hpp>
#include <sgpp/datadriven/operation/hash/simple/OperationTestModBspline.hpp>
#include <sgpp/datadriven/operation/hash/simple/OperationTestModLinear.hpp>
#include <sgpp/datadriven/operation/hash/simple/OperationTestModPoly.hpp>
#include <sgpp/datadriven/operation/hash/simple/OperationTestModWavelet.hpp>
#include <sgpp/datadriven/operation/hash/simple/OperationTestPoly.hpp>
#include <sgpp/datadriven/operation/hash/simple/OperationTestPrewavelet.hpp>

#include <sgpp/datadriven/operation/hash/OperationMultiEvalModMaskStreaming/OperationMultiEvalModMaskStreaming.hpp>
#include <sgpp/datadriven/operation/hash/OperationMultiEvalStreaming/OperationMultiEvalStreaming.hpp>

#ifdef __AVX__
#include <sgpp/datadriven/operation/hash/OperationMultipleEvalSubspace/combined/OperationMultipleEvalSubspaceCombined.hpp>
#include <sgpp/datadriven/operation/hash/OperationMultipleEvalSubspace/simple/OperationMultipleEvalSubspaceSimple.hpp>
#endif

#ifdef USE_OCL
#include "operation/hash/OperationMultipleEvalStreamingBSplineOCL/StreamingBSplineOCLOperatorFactory.hpp"
#include "operation/hash/OperationMultipleEvalStreamingModOCLFastMultiPlattform/OperatorFactory.hpp"
#include "operation/hash/OperationMultipleEvalStreamingModOCLMaskMultiPlatform/OperatorFactory.hpp"
#include "operation/hash/OperationMultipleEvalStreamingModOCLOpt/OperatorFactory.hpp"
#include "operation/hash/OperationMultipleEvalStreamingOCLMultiPlatform/OperatorFactory.hpp"

#include "operation/hash/OperationCreateGraphOCL/OpFactory.hpp"
#include "operation/hash/OperationDensityOCLMultiPlatform/OpFactory.hpp"
#include "operation/hash/OperationPruneGraphOCL/OpFactory.hpp"
#endif

#ifdef USE_MPI
#include "operation/hash/OperationMultiEvalMPI/OperationMultiEvalMPI.hpp"
#endif

#ifdef USE_HPX
#include "operation/hash/OperationMultiEvalHPX/OperationMultiEvalHPX.hpp"
#endif

#include <sgpp/optimization/function/scalar/ScalarFunction.hpp>

#ifdef USE_CUDA
#include "operation/hash/OperationMultiEvalCuda/OperationMultiEvalCuda.hpp"
#endif

#ifdef USE_SCALAPACK
#include <sgpp/datadriven/operation/hash/OperationMultipleEvalScalapack/OperationMultipleEvalLinearDistributed.hpp>
#include <sgpp/datadriven/operation/hash/OperationMultipleEvalScalapack/OperationMultipleEvalModLinearDistributed.hpp>
#endif

#include <sgpp/base/operation/BaseOpFactory.hpp>
#include <sgpp/globaldef.hpp>

#include <cstring>

namespace sgpp {
namespace op_factory {

datadriven::OperationTest* createOperationTest(base::Grid& grid) {
  if (grid.getType() == base::GridType::Linear) {
    return new datadriven::OperationTestLinear(&grid.getStorage());
  } else if (grid.getType() == base::GridType::LinearL0Boundary ||
             grid.getType() == base::GridType::LinearBoundary) {
    return new datadriven::OperationTestLinearBoundary(&grid.getStorage());
  } else if (grid.getType() == base::GridType::ModBspline) {
    return new datadriven::OperationTestModBspline(&grid.getStorage(),
                                                   ((base::ModBsplineGrid*)&grid)->getDegree());
  } else if (grid.getType() == base::GridType::ModLinear) {
    return new datadriven::OperationTestModLinear(&grid.getStorage());
  } else if (grid.getType() == base::GridType::Poly) {
    return new datadriven::OperationTestPoly(&grid.getStorage(),
                                             ((base::PolyGrid*)&grid)->getDegree());
  } else if (grid.getType() == base::GridType::ModPoly) {
    return new datadriven::OperationTestModPoly(&grid.getStorage(),
                                                ((base::ModPolyGrid*)&grid)->getDegree());
  } else if (grid.getType() == base::GridType::ModWavelet) {
    return new datadriven::OperationTestModWavelet(&grid.getStorage());
  } else if (grid.getType() == base::GridType::Prewavelet) {
    return new datadriven::OperationTestPrewavelet(&grid.getStorage());
  } else if (grid.getType() == base::GridType::LinearStretched) {
    return new datadriven::OperationTestLinearStretched(&grid.getStorage());
  } else if (grid.getType() == base::GridType::LinearStretchedBoundary) {
    return new datadriven::OperationTestLinearStretchedBoundary(&grid.getStorage());
  } else {
    throw base::factory_exception("OperationTest is not implemented for this grid type.");
  }
}

base::OperationMatrix* createOperationRegularizationDiagonal(base::Grid& grid, int mode, double k) {
  if (grid.getType() == base::GridType::Linear ||
      grid.getType() == base::GridType::LinearL0Boundary ||
      grid.getType() == base::GridType::LinearBoundary ||
      grid.getType() == base::GridType::ModLinear) {
    return new datadriven::OperationRegularizationDiagonalLinearBoundary(&grid.getStorage(), mode,
                                                                         k);
  } else {
    throw base::factory_exception(
        "OperationRegularizationDiagonal is not implemented for this grid type.");
  }
}

datadriven::OperationDensityMarginalize* createOperationDensityMarginalize(base::Grid& grid) {
  if (grid.getType() == base::GridType::Linear)
    return new datadriven::OperationDensityMarginalizeLinear(&grid);
  else if (grid.getType() == base::GridType::LinearBoundary ||
           grid.getType() == base::GridType::ModLinear || grid.getType() == base::GridType::Poly ||
           grid.getType() == base::GridType::ModPoly ||
           grid.getType() == base::GridType::PolyBoundary ||
           grid.getType() == base::GridType::PolyClenshawCurtis ||
           grid.getType() == base::GridType::PolyClenshawCurtisBoundary ||
           grid.getType() == base::GridType::ModPolyClenshawCurtis ||
           grid.getType() == base::GridType::Bspline ||
           grid.getType() == base::GridType::ModBspline ||
           grid.getType() == base::GridType::BsplineBoundary ||
           grid.getType() == base::GridType::BsplineClenshawCurtis ||
           grid.getType() == base::GridType::ModBsplineClenshawCurtis)
    return new datadriven::OperationDensityMarginalize(&grid);
  else
    throw base::factory_exception(
        "OperationDensityMarginalize is not implemented for this grid type.");
}

datadriven::OperationDensityMargTo1D* createOperationDensityMargTo1D(base::Grid& grid) {
  if (grid.getType() == base::GridType::Linear ||
      grid.getType() == base::GridType::LinearBoundary ||
      grid.getType() == base::GridType::ModLinear || grid.getType() == base::GridType::Poly ||
      grid.getType() == base::GridType::ModPoly || grid.getType() == base::GridType::PolyBoundary ||
      grid.getType() == base::GridType::PolyClenshawCurtis ||
      grid.getType() == base::GridType::PolyClenshawCurtisBoundary ||
      grid.getType() == base::GridType::ModPolyClenshawCurtis ||
      grid.getType() == base::GridType::Bspline || grid.getType() == base::GridType::ModBspline ||
      grid.getType() == base::GridType::BsplineBoundary ||
      grid.getType() == base::GridType::BsplineClenshawCurtis ||
      grid.getType() == base::GridType::ModBsplineClenshawCurtis)
    return new datadriven::OperationDensityMargTo1D(&grid);
  else
    throw base::factory_exception(
        "OperationDensityMargTo1D is not implemented for this grid type.");
}

datadriven::OperationDensitySampling1D* createOperationDensitySampling1D(base::Grid& grid) {
  if (grid.getType() == base::GridType::Linear)
    return new datadriven::OperationDensitySampling1DLinear(&grid);
  else
    throw base::factory_exception(
        "OperationDensitySampling1D is not implemented for this grid type.");
}

datadriven::OperationDensitySampling* createOperationDensitySampling(base::Grid& grid) {
  if (grid.getType() == base::GridType::Linear)
    return new datadriven::OperationDensitySamplingLinear(&grid);
  else
    throw base::factory_exception(
        "OperationDensitySampling is not implemented for this grid type.");
}

datadriven::OperationDensityRejectionSampling* createOperationDensityRejectionSampling(
    base::Grid& grid) {
  if (grid.getType() == base::GridType::Linear)
    return new datadriven::OperationDensityRejectionSamplingLinear(&grid);
  else
    throw base::factory_exception(
        "OperationDensityRejectionSampling is not implemented for this grid type.");
}

datadriven::OperationDensityConditional* createOperationDensityConditional(base::Grid& grid) {
  if (grid.getType() == base::GridType::Linear)
    return new datadriven::OperationDensityConditionalLinear(&grid);
  else if (grid.getType() == base::GridType::LinearBoundary ||
           grid.getType() == base::GridType::Poly || grid.getType() == base::GridType::ModPoly ||
           grid.getType() == base::GridType::PolyBoundary ||
           grid.getType() == base::GridType::PolyClenshawCurtis ||
           grid.getType() == base::GridType::PolyClenshawCurtisBoundary ||
           grid.getType() == base::GridType::ModPolyClenshawCurtis ||
           grid.getType() == base::GridType::Bspline ||
           grid.getType() == base::GridType::ModBspline ||
           grid.getType() == base::GridType::BsplineBoundary ||
           grid.getType() == base::GridType::BsplineClenshawCurtis ||
           grid.getType() == base::GridType::ModBsplineClenshawCurtis)
    return new datadriven::OperationDensityConditional(&grid);
  else
    throw base::factory_exception(
        "OperationDensityConditional is not implemented for this grid type.");
}

datadriven::OperationRosenblattTransformation* createOperationRosenblattTransformation(
    base::Grid& grid) {
  if (grid.getType() == base::GridType::Linear)
    return new datadriven::OperationRosenblattTransformationLinear(&grid);
  else if (grid.getType() == base::GridType::Poly)
    return new datadriven::OperationRosenblattTransformationPoly(&grid);
  else if (grid.getType() == base::GridType::ModPoly)
    return new datadriven::OperationRosenblattTransformationModPoly(&grid);
  else if (grid.getType() == base::GridType::PolyBoundary)
    return new datadriven::OperationRosenblattTransformationPolyBoundary(&grid);
  else if (grid.getType() == base::GridType::PolyClenshawCurtis)
    return new datadriven::OperationRosenblattTransformationPolyClenshawCurtis(&grid);
  else if (grid.getType() == base::GridType::ModPolyClenshawCurtis)
    return new datadriven::OperationRosenblattTransformationModPolyClenshawCurtis(&grid);
  else if (grid.getType() == base::GridType::PolyClenshawCurtisBoundary)
    return new datadriven::OperationRosenblattTransformationPolyClenshawCurtisBoundary(&grid);
  else if (grid.getType() == base::GridType::Bspline)
    return new datadriven::OperationRosenblattTransformationBspline(&grid);
  else if (grid.getType() == base::GridType::ModBspline)
    return new datadriven::OperationRosenblattTransformationModBspline(&grid);
  else if (grid.getType() == base::GridType::BsplineBoundary)
    return new datadriven::OperationRosenblattTransformationBsplineBoundary(&grid);
  else if (grid.getType() == base::GridType::BsplineClenshawCurtis)
    return new datadriven::OperationRosenblattTransformationBsplineClenshawCurtis(&grid);
  else if (grid.getType() == base::GridType::ModBsplineClenshawCurtis)
    return new datadriven::OperationRosenblattTransformationModBsplineClenshawCurtis(&grid);
  else
    throw base::factory_exception(
        "OperationRosenblattTransformation is not implemented for this grid type.");
}

datadriven::OperationTransformation1D* createOperationRosenblattTransformation1D(base::Grid& grid) {
  if (grid.getType() == base::GridType::Linear)
    return new datadriven::OperationRosenblattTransformation1DLinear(&grid);
  else if (grid.getType() == base::GridType::Poly)
    return new datadriven::OperationRosenblattTransformation1DPoly(&grid);
  else if (grid.getType() == base::GridType::ModPoly)
    return new datadriven::OperationRosenblattTransformation1DModPoly(&grid);
  else if (grid.getType() == base::GridType::PolyBoundary)
    return new datadriven::OperationRosenblattTransformation1DPolyBoundary(&grid);
  else if (grid.getType() == base::GridType::PolyClenshawCurtis)
    return new datadriven::OperationRosenblattTransformation1DPolyClenshawCurtis(&grid);
  else if (grid.getType() == base::GridType::ModPolyClenshawCurtis)
    return new datadriven::OperationRosenblattTransformation1DModPolyClenshawCurtis(&grid);
  else if (grid.getType() == base::GridType::PolyClenshawCurtisBoundary)
    return new datadriven::OperationRosenblattTransformation1DPolyClenshawCurtisBoundary(&grid);
  else if (grid.getType() == base::GridType::Bspline)
    return new datadriven::OperationRosenblattTransformation1DBspline(&grid);
  else if (grid.getType() == base::GridType::ModBspline)
    return new datadriven::OperationRosenblattTransformation1DModBspline(&grid);
  else if (grid.getType() == base::GridType::BsplineBoundary)
    return new datadriven::OperationRosenblattTransformation1DBsplineBoundary(&grid);
  else if (grid.getType() == base::GridType::BsplineClenshawCurtis)
    return new datadriven::OperationRosenblattTransformation1DBsplineClenshawCurtis(&grid);
  else if (grid.getType() == base::GridType::ModBsplineClenshawCurtis)
    return new datadriven::OperationRosenblattTransformation1DModBsplineClenshawCurtis(&grid);
  else
    throw base::factory_exception(
        "OperationRosenblattTransformation1D is not implemented for this grid type.");
}

datadriven::OperationInverseRosenblattTransformation*
createOperationInverseRosenblattTransformation(base::Grid& grid) {
  if (grid.getType() == base::GridType::Linear)
    return new datadriven::OperationInverseRosenblattTransformationLinear(&grid);
  else if (grid.getType() == base::GridType::Poly)
    return new datadriven::OperationInverseRosenblattTransformationPoly(&grid);
  else if (grid.getType() == base::GridType::ModPoly)
    return new datadriven::OperationInverseRosenblattTransformationModPoly(&grid);
  else if (grid.getType() == base::GridType::PolyBoundary)
    return new datadriven::OperationInverseRosenblattTransformationPolyBoundary(&grid);
  else if (grid.getType() == base::GridType::PolyClenshawCurtis)
    return new datadriven::OperationInverseRosenblattTransformationPolyClenshawCurtis(&grid);
  else if (grid.getType() == base::GridType::ModPolyClenshawCurtis)
    return new datadriven::OperationInverseRosenblattTransformationModPolyClenshawCurtis(&grid);
  else if (grid.getType() == base::GridType::PolyClenshawCurtisBoundary)
    return new datadriven::OperationInverseRosenblattTransformationPolyClenshawCurtisBoundary(
        &grid);
  else if (grid.getType() == base::GridType::Bspline)
    return new datadriven::OperationInverseRosenblattTransformationBspline(&grid);
  else if (grid.getType() == base::GridType::ModBspline)
    return new datadriven::OperationInverseRosenblattTransformationModBspline(&grid);
  else if (grid.getType() == base::GridType::BsplineBoundary)
    return new datadriven::OperationInverseRosenblattTransformationBsplineBoundary(&grid);
  else if (grid.getType() == base::GridType::BsplineClenshawCurtis)
    return new datadriven::OperationInverseRosenblattTransformationBsplineClenshawCurtis(&grid);
  else if (grid.getType() == base::GridType::ModBsplineClenshawCurtis)
    return new datadriven::OperationInverseRosenblattTransformationModBsplineClenshawCurtis(&grid);
  else
    throw base::factory_exception(
        "OperationInverseRosenblattTransformation is not implemented for this grid type.");
}

datadriven::OperationTransformation1D* createOperationInverseRosenblattTransformation1D(
    base::Grid& grid) {
  if (grid.getType() == base::GridType::Linear)
    return new datadriven::OperationInverseRosenblattTransformation1DLinear(&grid);
  else if (grid.getType() == base::GridType::Poly)
    return new datadriven::OperationInverseRosenblattTransformation1DPoly(&grid);
  else if (grid.getType() == base::GridType::ModPoly)
    return new datadriven::OperationInverseRosenblattTransformation1DModPoly(&grid);
  else if (grid.getType() == base::GridType::PolyBoundary)
    return new datadriven::OperationInverseRosenblattTransformation1DPolyBoundary(&grid);
  else if (grid.getType() == base::GridType::PolyClenshawCurtis)
    return new datadriven::OperationInverseRosenblattTransformation1DPolyClenshawCurtis(&grid);
  else if (grid.getType() == base::GridType::ModPolyClenshawCurtis)
    return new datadriven::OperationInverseRosenblattTransformation1DModPolyClenshawCurtis(&grid);
  else if (grid.getType() == base::GridType::PolyClenshawCurtisBoundary)
    return new datadriven::OperationInverseRosenblattTransformation1DPolyClenshawCurtisBoundary(
        &grid);
  else if (grid.getType() == base::GridType::Bspline)
    return new datadriven::OperationInverseRosenblattTransformation1DBspline(&grid);
  else if (grid.getType() == base::GridType::ModBspline)
    return new datadriven::OperationInverseRosenblattTransformation1DModBspline(&grid);
  else if (grid.getType() == base::GridType::BsplineBoundary)
    return new datadriven::OperationInverseRosenblattTransformation1DBsplineBoundary(&grid);
  else if (grid.getType() == base::GridType::BsplineClenshawCurtis)
    return new datadriven::OperationInverseRosenblattTransformation1DBsplineClenshawCurtis(&grid);
  else if (grid.getType() == base::GridType::ModBsplineClenshawCurtis)
    return new datadriven::OperationInverseRosenblattTransformation1DModBsplineClenshawCurtis(
        &grid);
  else
    throw base::factory_exception(
        "OperationInverseRosenblattTransformation1D is not implemented for this grid type.");
}

datadriven::OperationInverseRosenblattTransformationKDE*
createOperationInverseRosenblattTransformationKDE(datadriven::KernelDensityEstimator& kde) {
  return new datadriven::OperationInverseRosenblattTransformationKDE(kde);
}

datadriven::OperationRosenblattTransformationKDE* createOperationRosenblattTransformationKDE(
    datadriven::KernelDensityEstimator& kde) {
  return new datadriven::OperationRosenblattTransformationKDE(kde);
}

datadriven::OperationDensityMarginalizeKDE* createOperationDensityMarginalizeKDE(
    datadriven::KernelDensityEstimator& kde) {
  return new datadriven::OperationDensityMarginalizeKDE(kde);
}

datadriven::OperationDensityConditionalKDE* createOperationDensityConditionalKDE(
    datadriven::KernelDensityEstimator& kde) {
  return new datadriven::OperationDensityConditionalKDE(kde);
}

base::OperationMultipleEval* createOperationMultipleEval(
    base::Grid& grid, base::DataMatrix& dataset,
    sgpp::datadriven::OperationMultipleEvalConfiguration& configuration) {
  if (configuration.getMPIType() == sgpp::datadriven::OperationMultipleEvalMPIType::MASTERSLAVE) {
#ifdef USE_MPI
    if (grid.getType() == base::GridType::Linear) {
      return new datadriven::OperationMultiEvalMPI(
          grid, dataset, sgpp::datadriven::OperationMultipleEvalType::STREAMING,
          sgpp::datadriven::OperationMultipleEvalSubType::DEFAULT);
    }
#else
    throw base::factory_exception(
        "Error creating function: the library wasn't compiled with MPI support");
#endif
  } else if (configuration.getMPIType() == sgpp::datadriven::OperationMultipleEvalMPIType::HPX) {
#ifdef USE_HPX
    if (grid.getType() == base::GridType::Linear) {
      return new datadriven::OperationMultiEvalHPX(grid, dataset, configuration);
    }
#else
    throw base::factory_exception(
        "Error creating function: the library wasn't compiled with HPX support");
#endif
  }

  // can now assume that MPI type is NONE
  if (configuration.getType() == sgpp::datadriven::OperationMultipleEvalType::DEFAULT) {
    return createOperationMultipleEval(grid, dataset);
  }

  if (grid.getType() == base::GridType::Linear) {
    if (configuration.getType() == datadriven::OperationMultipleEvalType::DEFAULT ||
        configuration.getType() == datadriven::OperationMultipleEvalType::STREAMING) {
      if (configuration.getSubType() == sgpp::datadriven::OperationMultipleEvalSubType::DEFAULT) {
        return new datadriven::OperationMultiEvalStreaming(grid, dataset);
      }
      if (configuration.getSubType() == sgpp::datadriven::OperationMultipleEvalSubType::OCLMP) {
#ifdef USE_OCL
        return datadriven::createStreamingOCLMultiPlatformConfigured(grid, dataset, configuration);
#else
        throw base::factory_exception(
            "Error creating function: the library wasn't compiled with OpenCL support");
#endif
      }
    } else if (configuration.getType() == datadriven::OperationMultipleEvalType::SUBSPACELINEAR) {
      if (configuration.getSubType() == sgpp::datadriven::OperationMultipleEvalSubType::DEFAULT ||
          configuration.getSubType() == sgpp::datadriven::OperationMultipleEvalSubType::COMBINED) {
#ifdef __AVX__
        return new datadriven::OperationMultipleEvalSubspaceCombined(grid, dataset);
#else
        throw base::factory_exception(
            "Error creating function: the library wasn't compiled with AVX");
#endif
      } else if (configuration.getSubType() ==
                 sgpp::datadriven::OperationMultipleEvalSubType::SIMPLE) {
#ifdef __AVX__
        return new datadriven::OperationMultipleEvalSubspaceSimple(grid, dataset);
#else
        throw base::factory_exception(
            "Error creating function: the library wasn't compiled with AVX");
#endif
      }
    } else if (configuration.getType() == datadriven::OperationMultipleEvalType::SCALAPACK) {
#ifdef USE_SCALAPACK
      return new datadriven::OperationMultipleEvalLinearDistributed(grid, dataset);
#else
      throw base::factory_exception(
          "Error creating function: the library wasn't compiled with ScaLAPACK support");
#endif
    }
  } else if (grid.getType() == base::GridType::ModLinear) {
    if (configuration.getType() == datadriven::OperationMultipleEvalType::STREAMING) {
      if (configuration.getSubType() == sgpp::datadriven::OperationMultipleEvalSubType::DEFAULT) {
        return new datadriven::OperationMultiEvalModMaskStreaming(grid, dataset);
      }
      if (configuration.getSubType() == sgpp::datadriven::OperationMultipleEvalSubType::OCLFASTMP) {
#ifdef USE_OCL
        return datadriven::createStreamingModOCLFastMultiPlatformConfigured(grid, dataset,
                                                                            configuration);
#else
        throw base::factory_exception(
            "Error creating function: the library wasn't compiled with OpenCL support");
#endif
      } else if (configuration.getSubType() ==
                 sgpp::datadriven::OperationMultipleEvalSubType::OCLMASKMP) {
#ifdef USE_OCL
        return datadriven::createStreamingModOCLMaskMultiPlatformConfigured(grid, dataset,
                                                                            configuration);
#else
        throw base::factory_exception(
            "Error creating function: the library wasn't compiled with OpenCL support");
#endif
      } else if (configuration.getSubType() ==
                 sgpp::datadriven::OperationMultipleEvalSubType::OCLOPT) {
#ifdef USE_OCL
        return datadriven::createStreamingModOCLOptConfigured(grid, dataset, configuration);
#else
        throw base::factory_exception(
            "Error creating function: the library wasn't compiled with OpenCL support");
#endif
      }
    } else if (configuration.getType() == datadriven::OperationMultipleEvalType::SCALAPACK) {
#ifdef USE_SCALAPACK
      return new datadriven::OperationMultipleEvalModLinearDistributed(grid, dataset);
#else
      throw base::factory_exception(
          "Error creating function: the library wasn't compiled with ScaLAPACK support");
#endif
    }
  } else if (grid.getType() == base::GridType::Bspline) {
    if (configuration.getType() == datadriven::OperationMultipleEvalType::STREAMING) {
      if (configuration.getSubType() == sgpp::datadriven::OperationMultipleEvalSubType::OCL) {
#ifdef USE_OCL
        return datadriven::createStreamingBSplineOCLConfigured(grid, dataset, configuration);
#else
        throw base::factory_exception(
            "Error creating function: the library wasn't compiled with OpenCL support");
#endif
      }
    }
  } else if (grid.getType() == base::GridType::Poly) {
    if (configuration.getType() == datadriven::OperationMultipleEvalType::DEFAULT) {
      if (configuration.getSubType() == sgpp::datadriven::OperationMultipleEvalSubType::CUDA) {
#ifdef USE_CUDA
        return new datadriven::OperationMultiEvalCuda(grid, dataset, grid.getDegree(), false);
#else
        throw base::factory_exception(
            "Error creating function: the library wasn't compiled with CUDA support");
#endif
      }
    } else if (configuration.getType() == datadriven::OperationMultipleEvalType::MORTONORDER) {
      if (configuration.getSubType() == sgpp::datadriven::OperationMultipleEvalSubType::CUDA) {
#ifdef USE_CUDA
        return new datadriven::OperationMultiEvalCuda(grid, dataset, grid.getDegree(), true);
#else
        throw base::factory_exception(
            "Error creating function: the library wasn't compiled with CUDA support");
#endif
      }
    }
  }

  throw base::factory_exception("OperationMultiEval is not implemented for this grid type.");
}

datadriven::OperationMakePositive* createOperationMakePositive(
    datadriven::MakePositiveCandidateSearchAlgorithm candidateSearchAlgorithm,
    datadriven::MakePositiveInterpolationAlgorithm interpolationAlgorithm,
    bool generateConsistentGrid, bool verbose, sgpp::optimization::ScalarFunction* f) {
  return new datadriven::OperationMakePositive(candidateSearchAlgorithm, interpolationAlgorithm,
                                               generateConsistentGrid, verbose, f);
}

datadriven::OperationLimitFunctionValueRange* createOperationLimitFunctionValueRange(
    datadriven::MakePositiveCandidateSearchAlgorithm candidateSearchAlgorithm,
    datadriven::MakePositiveInterpolationAlgorithm interpolationAlgorithm, bool verbose,
    sgpp::optimization::ScalarFunction* f) {
  return new datadriven::OperationLimitFunctionValueRange(candidateSearchAlgorithm,
                                                          interpolationAlgorithm, verbose, f);
}

datadriven::OperationCovariance* createOperationCovariance(base::Grid& grid) {
  return new datadriven::OperationCovariance(grid);
}

}  // namespace op_factory
}  // namespace sgpp
