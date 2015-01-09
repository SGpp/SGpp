/* ****************************************************************************
 * Copyright (C) 2009 Technische Universitaet Muenchen                         *
 * This file is part of the SG++ project. For conditions of distribution and   *
 * use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
 **************************************************************************** */
// @author Valeriy Khakhutskyy (khakhutv@in.tum.de), Dirk Pflueger (pflueged@in.tum.de)
#ifndef DATADRIVEN_OP_FACTORY_HPP
#define DATADRIVEN_OP_FACTORY_HPP

#include "base/grid/Grid.hpp"

#include "datadriven/operation/OperationTest.hpp"
#include "base/operation/OperationMatrix.hpp"
#include "datadriven/operation/OperationDensityMarginalize.hpp"
#include "datadriven/operation/OperationDensityMargTo1D.hpp"
#include "datadriven/operation/OperationDensityConditional.hpp"
#include "datadriven/operation/OperationDensitySampling1D.hpp"
#include "datadriven/operation/OperationDensitySampling.hpp"
#include "datadriven/operation/OperationDensityRejectionSampling.hpp"
#include "datadriven/operation/OperationRosenblattTransformation.hpp"
#include "datadriven/operation/OperationTransformation1D.hpp"
#include "datadriven/operation/OperationInverseRosenblattTransformation.hpp"
#include "base/operation/OperationMultipleEval.hpp"
#include "datadriven/operation/DatadrivenOperationCommon.hpp"

/*
 * This file contains factory methods for operations.
 */

namespace sg {
namespace op_factory {
/**
 * Factory method, returning an OperationTest for the grid at hand.
 * Note: object has to be freed after use.
 *
 * @param grid Grid which is to be used for OperationTest
 * @return Pointer to the new OperationTest object for the Grid grid
 */
datadriven::OperationTest* createOperationTest(base::Grid& grid);

/**
 * Factory method, returning an OperationRegularizationDiagonal for the grid at hand.
 * Note: object has to be freed after use.
 *
 * @param grid Grid which is to be used for OperationRegularizationDiagonal
 * @param mode Mode, specifying which regularization to use. Example: OperationRegularizationDiagonal::HKMIX.
 * @param k Parameter for @f$H^k@f$
 * @return Pointer to the new OperationRegularizationDiagonal object for the Grid grid
 */
base::OperationMatrix* createOperationRegularizationDiagonal(base::Grid& grid, int mode, double k);

/**
 * Factory method, returning an OperationDensityMarginalize for the grid.
 * Note: object has to be freed after use.
 *
 * @param grid Grid which is to be used for the operation
 * @return Pointer to new OperationDensityMarginalize for the Grid grid
 */
datadriven::OperationDensityMarginalize* createOperationDensityMarginalize(base::Grid& grid);

/**
 * Factory method, returning an OperationDensityMargTo1D for the grid.
 * Note: object has to be freed after use.
 *
 * @param grid Grid which is to be used for the operation
 * @return Pointer to new OperationDensityMarginalize for the Grid grid
 */
datadriven::OperationDensityMargTo1D* createOperationDensityMargTo1D(base::Grid& grid);

/**
 * Factory method, returning an OperationDensitySampling1D for the grid.
 * Note: object has to be freed after use.
 *
 * @param grid Grid which is to be used for the operation
 * @return Pointer to new OperationDensitySampling1D for the Grid grid
 */
datadriven::OperationDensitySampling1D* createOperationDensitySampling1D(base::Grid& grid);

/**
 * Factory method, returning an OperationDensitySampling for the grid.
 * Note: object has to be freed after use.
 *
 * @param grid Grid which is to be used for the operation
 * @return Pointer to new OperationDensitySampling1D for the Grid grid
 */
datadriven::OperationDensitySampling* createOperationDensitySampling(base::Grid& grid);

/**
 * Factory method, returning an OperationDensityRejectionSampling for the grid.
 * Note: object has to be freed after use.
 *
 * @param grid Grid which is to be used for the operation
 * @return Pointer to new OperationDensitySampling1D for the Grid grid
 */
datadriven::OperationDensityRejectionSampling* createOperationDensityRejectionSampling(base::Grid& grid);

/**
 * Factory method, returning an OperationDensityConditional for the grid.
 * Note: object has to be freed after use.
 *
 * @param grid Grid which is to be used for the operation
 * @return Pointer to new OperationDensityConditional for the Grid grid
 */
datadriven::OperationDensityConditional* createOperationDensityConditional(base::Grid& grid);

/**
 * Factory method, returning an OperationRosenblattTransformation for the grid.
 * Note: object has to be freed after use.
 *
 * @param grid Grid which is to be used for the operation
 * @return Pointer to new OperationRosenblattTransformation for the Grid grid
 */
datadriven::OperationRosenblattTransformation* createOperationRosenblattTransformation(base::Grid& grid);

/**
 * Factory method, returning an OperationInverseRosenblattTransformation for the grid.
 * Note: object has to be freed after use.
 *
 * @param grid Grid which is to be used for the operation
 * @return Pointer to new OperationRosenblattTransformation for the Grid grid
 */
datadriven::OperationInverseRosenblattTransformation* createOperationInverseRosenblattTransformation(base::Grid &grid);

/**
 * Factory method, returning an OperationRosenblattTransformation1D for the grid.
 * Note: object has to be freed after use.
 *
 * @param grid Grid which is to be used for the operation
 * @return Pointer to new OperationRosenblattTransformation1D for the Grid grid
 */
datadriven::OperationTransformation1D* createOperationRosenblattTransformation1D(base::Grid &grid);

base::OperationMultipleEval* createOperationMultipleEval(base::Grid &grid, base::DataMatrix &dataset,
		sg::datadriven::OperationMultipleEvalConfiguration configuration);

}
}

#endif /*DATADRIVEN_OP_FACTORY_HPP*/
