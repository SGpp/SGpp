// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#pragma once

#include <sgpp/base/grid/Grid.hpp>

#include <sgpp/datadriven/operation/hash/simple/OperationTest.hpp>
#include <sgpp/base/operation/hash/OperationMatrix.hpp>
#include <sgpp/datadriven/operation/hash/simple/OperationDensityMarginalize.hpp>
#include <sgpp/datadriven/operation/hash/simple/OperationDensityMargTo1D.hpp>
#include <sgpp/datadriven/operation/hash/simple/OperationDensityConditional.hpp>
#include <sgpp/datadriven/operation/hash/simple/OperationDensitySampling1D.hpp>
#include <sgpp/datadriven/operation/hash/simple/OperationDensitySampling.hpp>
#include <sgpp/datadriven/operation/hash/simple/OperationDensityRejectionSampling.hpp>
#include <sgpp/datadriven/operation/hash/simple/OperationTransformation1D.hpp>
#include <sgpp/datadriven/operation/hash/simple/OperationRosenblattTransformation.hpp>
#include <sgpp/datadriven/operation/hash/simple/OperationInverseRosenblattTransformation.hpp>
#include <sgpp/datadriven/operation/hash/simple/OperationCovariance.hpp>
#include <sgpp/base/operation/hash/OperationMultipleEval.hpp>
#include <sgpp/datadriven/operation/hash/DatadrivenOperationCommon.hpp>

#include <sgpp/datadriven/operation/hash/simple/OperationLimitFunctionValueRange.hpp>
#include <sgpp/datadriven/operation/hash/simple/OperationMakePositive.hpp>

#include <sgpp/datadriven/operation/hash/simple/OperationRosenblattTransformationKDE.hpp>
#include <sgpp/datadriven/operation/hash/simple/OperationInverseRosenblattTransformationKDE.hpp>
#include <sgpp/datadriven/operation/hash/simple/OperationDensityMarginalizeKDE.hpp>
#include <sgpp/datadriven/operation/hash/simple/OperationDensityConditionalKDE.hpp>

#include <sgpp/globaldef.hpp>
#include <sgpp/base/function/scalar/ScalarFunction.hpp>

namespace sgpp {
namespace op_factory {
/**
 * Factory method, returning an OperationTest for the grid at hand.
 *
 * @param grid Grid which is to be used for OperationTest
 * @return Pointer to the new OperationTest object for the Grid grid
 */
datadriven::OperationTest* createOperationTest(base::Grid& grid);

/**
 * Factory method, returning an OperationRegularizationDiagonal for the grid at hand.
 *
 * @param grid Grid which is to be used for OperationRegularizationDiagonal
 * @param mode Mode, specifying which regularization to use. Example:
 * OperationRegularizationDiagonal::HKMIX.
 * @param k Parameter for @f$H^k@f$
 * @return Pointer to the new OperationRegularizationDiagonal object for the Grid grid
 */
base::OperationMatrix* createOperationRegularizationDiagonal(base::Grid& grid, int mode, double k);

/**
 * Factory method, returning an OperationDensityMarginalize for the grid.
 *
 * @param grid Grid which is to be used for the operation
 * @return Pointer to new OperationDensityMarginalize for the Grid grid
 */
datadriven::OperationDensityMarginalize* createOperationDensityMarginalize(base::Grid& grid);

/**
 * Factory method, returning an OperationDensityMargTo1D for the grid.
 *
 * @param grid Grid which is to be used for the operation
 * @return Pointer to new OperationDensityMarginalize for the Grid grid
 */
datadriven::OperationDensityMargTo1D* createOperationDensityMargTo1D(base::Grid& grid);

/**
 * Factory method, returning an OperationDensitySampling1D for the grid.
 *
 * @param grid Grid which is to be used for the operation
 * @return Pointer to new OperationDensitySampling1D for the Grid grid
 */
datadriven::OperationDensitySampling1D* createOperationDensitySampling1D(base::Grid& grid);

/**
 * Factory method, returning an OperationDensitySampling for the grid.
 *
 * @param grid Grid which is to be used for the operation
 * @return Pointer to new OperationDensitySampling1D for the Grid grid
 */
datadriven::OperationDensitySampling* createOperationDensitySampling(base::Grid& grid);

/**
 * Factory method, returning an OperationDensityRejectionSampling for the grid.
 *
 * @param grid Grid which is to be used for the operation
 * @return Pointer to new OperationDensitySampling1D for the Grid grid
 */
datadriven::OperationDensityRejectionSampling* createOperationDensityRejectionSampling(
    base::Grid& grid);

/**
 * Factory method, returning an OperationDensityConditional for the grid.
 *
 * @param grid Grid which is to be used for the operation
 * @return Pointer to new OperationDensityConditional for the Grid grid
 */
datadriven::OperationDensityConditional* createOperationDensityConditional(base::Grid& grid);

/**
 * Factory method, returning an OperationRosenblattTransformation for the grid.
 *
 * @param grid Grid which is to be used for the operation
 * @return Pointer to new OperationRosenblattTransformation for the Grid grid
 */
datadriven::OperationRosenblattTransformation* createOperationRosenblattTransformation(
    base::Grid& grid);

/**
 * Factory method, returning an OperationRosenblattTransformation for the grid.
 *
 * @param grid Grid which is to be used for the operation
 * @return Pointer to new OperationRosenblattTransformation1D for the Grid grid
 */
datadriven::OperationTransformation1D* createOperationRosenblattTransformation1D(base::Grid& grid);

/**
 * Factory method, returning an OperationInverseRosenblattTransformation for the grid.
 *
 * @param grid Grid which is to be used for the operation
 * @return Pointer to new OperationInverseRosenblattTransformation for the Grid grid
 */
datadriven::OperationInverseRosenblattTransformation*
createOperationInverseRosenblattTransformation(base::Grid& grid);

/**
 * Factory method, returning an OperationInverseRosenblattTransformation1D for the grid.
 * @param grid Grid which is to be used for the operation
 * @return Pointer to new OperationInverseRosenblattTransformation1D for the Grid grid
 */
datadriven::OperationTransformation1D* createOperationInverseRosenblattTransformation1D(
    base::Grid& grid);

/**
 * Factory method, returning an OperationRosenblattTransformationKDE for the kde.
 *
 * @param kde KernelDensityEstimator for which the Rosenblatt transformation should be computed
 * @return Pointer to new OperationRosenblattTransformationKDE for the kde
 */
datadriven::OperationRosenblattTransformationKDE* createOperationRosenblattTransformationKDE(
    datadriven::KernelDensityEstimator& kde);

/**
 * Factory method, returning an OperationInverseRosenblattTransformationKDE for the kde.
 *
 * @param kde KernelDensityEstimator for which the inverse Rosenblatt transformation should be
 * computed
 * @return Pointer to new OperationInverseRosenblattTransformationKDE for the kde
 */
datadriven::OperationInverseRosenblattTransformationKDE*
createOperationInverseRosenblattTransformationKDE(datadriven::KernelDensityEstimator& kde);

/**
 * Factory method, returning an OperationDensityMarginalizeKDE for the kernel density.
 *
 * @param kde kernel density which is to be used for the operation
 * @return Pointer to new OperationDensityMarginalizeKDE
 */
datadriven::OperationDensityMarginalizeKDE* createOperationDensityMarginalizeKDE(
    datadriven::KernelDensityEstimator& kde);

/**
 * Factory method, returning an OperationDensityConditionalKDE for the kernel density.
 *
 * @param kde kernel density which is to be used for the operation
 * @return Pointer to new OperationDensityConditionalKDE
 */
datadriven::OperationDensityConditionalKDE* createOperationDensityConditionalKDE(
    datadriven::KernelDensityEstimator& kde);

/**
 * Factory method, returning an OperationMultipleEval for the grid.
 *
 * @param grid Grid which is to be used for the operation
 * @param dataset dataset to be evaluated
 * @param configuration configuration to be used (evalType and evalSubType)
 * @return Pointer to new OperationMultipleEval for the Grid grid
 */
base::OperationMultipleEval* createOperationMultipleEval(
    base::Grid& grid, base::DataMatrix& dataset,
    sgpp::datadriven::OperationMultipleEvalConfiguration& configuration);

/**
 * Factory method, returning an OperationMakePositive for an arbitrary function f or some
 * sparse grid, which is yet to be defined.
 * Note: object has to be freed after use.
 *
 * @param candidateSearchAlgorithm defines algorithm for candidate set enumeration
 * @param interpolationAlgorithm defines algorithm for coefficient estimation of extension set
 * @param generateConsistentGrid if set to true, all hierarchical ancestors are available in the
 * resulting grid
 * @param verbose verbosity
 * @param f function to be approximated (as an alternative to a sparse grid function)
 *
 * @return Pointer to the new OperationMakePositive object for the Grid grid
 */
datadriven::OperationMakePositive* createOperationMakePositive(
    datadriven::MakePositiveCandidateSearchAlgorithm candidateSearchAlgorithm =
        datadriven::MakePositiveCandidateSearchAlgorithm::Intersections,
    datadriven::MakePositiveInterpolationAlgorithm interpolationAlgorithm =
        datadriven::MakePositiveInterpolationAlgorithm::SetToZero,
    bool generateConsistentGrid = true, sgpp::base::ScalarFunction* f = nullptr);

/**
 * Factory method, returning an OperationLimitFunctionValueRange for an arbitrary function f or some
 * sparse grid, which is yet to be defined.
 * Note: object has to be freed after use.
 *
 * @param candidateSearchAlgorithm defines algorithm for candidate set enumeration
 * @param interpolationAlgorithm defines algorithm for coefficient estimation of extension set
 * @param verbose verbosity
 * @param f function to be approximated (as an alternative to a sparse grid function)
 *
 * @return Pointer to the new OperationLimitFunctionValueRange object for the Grid grid
 */
datadriven::OperationLimitFunctionValueRange* createOperationLimitFunctionValueRange(
    datadriven::MakePositiveCandidateSearchAlgorithm candidateSearchAlgorithm =
        datadriven::MakePositiveCandidateSearchAlgorithm::Intersections,
    datadriven::MakePositiveInterpolationAlgorithm interpolationAlgorithm =
        datadriven::MakePositiveInterpolationAlgorithm::SetToZero,
    sgpp::base::ScalarFunction* f = nullptr);

/**
 * Factory method, returning an OperationCovariance for the grid at hand.
 * Note: object has to be freed after use.
 *
 * @param grid Grid which is to be used
 * @return Pointer to the new OperationCovariance object for the Grid grid
 */
datadriven::OperationCovariance* createOperationCovariance(base::Grid& grid);

}  // namespace op_factory
}  // namespace sgpp
