// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#ifndef BASE_OP_FACTORY_HPP
#define BASE_OP_FACTORY_HPP

#include <sgpp/base/grid/Grid.hpp>

#include <sgpp/base/operation/hash/OperationConvert.hpp>
#include <sgpp/base/operation/hash/OperationDiagonal.hpp>
#include <sgpp/base/operation/hash/OperationEval.hpp>
#include <sgpp/base/operation/hash/OperationEvalGradient.hpp>
#include <sgpp/base/operation/hash/OperationEvalHessian.hpp>
#include <sgpp/base/operation/hash/OperationEvalPartialDerivative.hpp>
#include <sgpp/base/operation/hash/OperationFirstMoment.hpp>
#include <sgpp/base/operation/hash/OperationHierarchisation.hpp>
#include <sgpp/base/operation/hash/OperationIdentity.hpp>
#include <sgpp/base/operation/hash/OperationMatrix.hpp>
#include <sgpp/base/operation/hash/OperationMultipleEval.hpp>
#include <sgpp/base/operation/hash/OperationQuadrature.hpp>
#include <sgpp/base/operation/hash/OperationSecondMoment.hpp>
#include <sgpp/base/operation/hash/OperationStencilHierarchisation.hpp>
#include <sgpp/base/operation/hash/OperationWeightedQuadrature.hpp>
#include <sgpp/base/operation/hash/OperationWeightedSecondMoment.hpp>

/*
 * This file contains factory methods for operations.
 */

#include <sgpp/globaldef.hpp>

#include <set>
#include <vector>

namespace sgpp {

namespace op_factory {
/**
 * Factory method, returning an OperationDiagonal (OperationMatrix) for the grid at hand.
 *
 * @param grid Grid which is to be used
 * @param multiplicationFactor MultiplicationFactor which is to be used
 * @return Pointer to the new OperationMatrix object for the Grid grid
 */
base::OperationMatrix* createOperationDiagonal(base::Grid& grid,
                                               double multiplicationFactor = 0.25);
/**
 * Factory method, returning an OperationHierarchisation for the grid at hand.
 * Note: object has to be freed after use.
 *
 * @param grid Grid which is to be used for hierarchisation
 * @return Pointer to the new OperationHierarchisation object for the Grid grid
 */
base::OperationHierarchisation* createOperationHierarchisation(base::Grid& grid);
/**
 * Factory method, returning an OperationArbitraryBoundaryHierarchisation for the grid at hand.
 * Note: object has to be freed after use. This operation should be used if the boundary level
 * of your grid is larger than 1.
 *
 * @param grid Grid which is to be used for hierarchisation
 * @return Pointer to the new OperationArbitraryBoundaryHierarchisation object for the Grid grid
 */
base::OperationHierarchisation* createOperationArbitraryBoundaryHierarchisation(base::Grid& grid);
/**
 * Factory method, returning an OperationQuadrature for the grid at hand.
 * Note: object has to be freed after use.
 *
 * @param grid Grid which is to be used for quadrature
 * @return Pointer to the new OperationQuadrature for the Grid grid
 */
base::OperationQuadrature* createOperationQuadrature(base::Grid& grid);
/**
 * Factory method, returning an OperationWeightedQuadrature for the grid at hand.
 * Note: object has to be freed after use.
 *
 * @param grid Grid which is to be used for quadrature
 * @param quadOrder quadrature order
 * @return Pointer to the new OperationWeightedQuadrature for the Grid grid
 */
base::OperationWeightedQuadrature* createOperationWeightedQuadrature(base::Grid& grid,
                                                                     size_t quadOrder);
/**
 * Factory method, returning an OperationFirstMoment for the grid at hand.
 * Note: object has to be freed after use.
 *
 * @param grid Grid which is to be used for quadrature
 * @return Pointer to the new OperationFirstMoment for the Grid grid
 */
base::OperationFirstMoment* createOperationFirstMoment(base::Grid& grid);
/**
 * Factory method, returning an OperationSecondMoment for the grid at hand.
 * Note: object has to be freed after use.
 *
 * @param grid Grid which is to be used for quadrature
 * @return Pointer to the new OperationSecondMoment for the Grid grid
 */
base::OperationSecondMoment* createOperationSecondMoment(base::Grid& grid);
/**
 * Factory method, returning an OperationWeightedSecondMoment for the grid at hand.
 * Note: object has to be freed after use.
 *
 * @param grid Grid which is to be used for quadrature
 * @param quadOrder quadrature order
 * @return Pointer to the new OperationSecondMoment for the Grid grid
 */
base::OperationWeightedSecondMoment* createOperationWeightedSecondMoment(base::Grid& grid,
                                                                         size_t quadOrder);
/**
 * Factory method, returning an OperationConvert for the grid at hand.
 * Note: object has to be freed after use.
 *
 * @param grid Grid which is to be used for conversion
 * @return Pointer to the new OperationConvert object for the Grid grid
 */
base::OperationConvert* createOperationConvert(base::Grid& grid);
/**
 * Factory method, returning an OperationIdentity for the grid at hand.
 * Note: object has to be freed after use.
 * Just calls OperationIdentity() independent of grid; factory method
 * provided for uniform use.
 *
 * @return Pointer to the new OperationIdentity object
 */
base::OperationMatrix* createOperationIdentity(base::Grid& grid);
/**
 * Factory method, returning an OperationEval for the grid at hand.
 * Note: object has to be freed after use.
 *
 * @param grid Grid which is to be used
 * @return Pointer to the new OperationEval object for the Grid grid
 */
base::OperationEval* createOperationEval(base::Grid& grid);
/**
 * Factory method, returning an OperationMultipleEval for the grid at hand.
 * Note: object has to be freed after use.
 *
 * @param grid Grid which is to be used
 * @param dataset The dataset (DataMatrix, one datapoint per row) that is to be evaluated for
 * the sparse grid function
 * @return Pointer to the new OperationMultipleEval object for the Grid grid
 */

base::OperationMultipleEval* createOperationMultipleEval(base::Grid& grid,
                                                         base::DataMatrix& dataset);
/**
 * Similar to createOperationMultipleEval, but makes use of interaction terms during evaluation
 *
 * @param grid Grid which is to be used
 * @param dataset The dataset (DataMatrix, one datapoint per row) that is to be evaluated for
 * the sparse grid function
 * @param interactions A list of Interaction the SG is reduced to
 * @return Pointer to the new OperationMultipleEval object for the Grid grid
 */
base::OperationMultipleEval* createOperationMultipleEvalInter(
    base::Grid& grid, base::DataMatrix& dataset, std::set<std::set<size_t>> interactions);

/**
 * Factory method, returning an OperationMultipleEvalNaive for the grid at hand.
 * Note: object has to be freed after use.
 *
 * @param grid Grid which is to be used
 * @param dataset The dataset (DataMatrix, one datapoint per row) that is to be evaluated for
 * the sparse grid function
 * @return Pointer to the new OperationMultipleEval object for the Grid grid
 */
base::OperationMultipleEval* createOperationMultipleEvalNaive(base::Grid& grid,
                                                              base::DataMatrix& dataset);

/**
 * Factory method, returning an OperationEval for the grid at hand.
 * In contrast to OperationEval, implementations of OperationEval
 * returned by this function should
 * use a "naive" method for evaluating sparse grid functions, e.g. evaluate
 * all basis functions by brute force.
 * Note: object has to be freed after use.
 *
 * @param grid Grid which is to be used
 * @return Pointer to the new OperationEval object for the Grid grid
 */
base::OperationEval* createOperationEvalNaive(base::Grid& grid);
/**
 * Factory method, returning an OperationEvalGradient for the grid at hand.
 * Implementations of OperationEvalGradientNaive returned by this function should
 * use a "naive" method for evaluating sparse grid function gradients, e.g. evaluate
 * all basis functions by brute force.
 * Note: object has to be freed after use.
 *
 * @param grid Grid which is to be used
 * @return Pointer to the new OperationEvalGradient object for the Grid grid
 */

base::OperationEvalGradient* createOperationEvalGradientNaive(base::Grid& grid);
/**
 * Factory method, returning an OperationEvalHessian for the grid at hand.
 * Implementations of OperationEvalHessianNaive returned by this function should
 * use a "naive" method for evaluating sparse grid function Hessians, e.g. evaluate
 * all basis functions by brute force.
 * Note: object has to be freed after use.
 *
 * @param grid Grid which is to be used
 * @return Pointer to the new OperationEvalHessian object for the Grid grid
 */

base::OperationEvalHessian* createOperationEvalHessianNaive(base::Grid& grid);
/**
 * Factory method, returning an OperationEvalPartialDerivative for the grid at hand.
 * Implementations of OperationEvalPartialDerivativeNaive returned by this function should
 * use a "naive" method for evaluating sparse grid function partial derivatives, e.g. evaluate
 * all basis functions by brute force.
 * Note: object has to be freed after use.
 *
 * @param grid Grid which is to be used
 * @return Pointer to the new OperationEvalPartialDerivative object for the Grid grid
 */
base::OperationEvalPartialDerivative* createOperationEvalPartialDerivativeNaive(base::Grid& grid);

}  // namespace op_factory
}  // namespace sgpp

#endif /*BASE_OP_FACTORY_HPP*/
