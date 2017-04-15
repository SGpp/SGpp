/* Copyright (C) 2008-today The SG++ project
 * This file is part of the SG++ project. For conditions of distribution and
 * use, please see the copyright notice provided with SG++ or at
 * sgpp.sparsegrids.org
 *
 * IChol.hpp
 *
 *  Created on: Nov 26, 2016
 *      Author: Michael Lettrich
 */

#pragma once

#include <sgpp/base/datatypes/DataVector.hpp>
#include <sgpp/datadriven/algorithm/SparseDataMatrix.hpp>

namespace sgpp {
namespace datadriven {

using sgpp::base::DataVector;

namespace IChol {

/**
 * decompose in place
 */
void decompose(const DataMatrix& matrix, SparseDataMatrix& result, size_t sweeps);

void decompose(const DataMatrix& matrix, DataMatrix& result, size_t sweeps);

/**
 * do a cholesky update on the last n rows
 */
void updateLastNRows(SparseDataMatrix& matrix, size_t numRows, size_t sweeps);

void normToUnitDiagonal(SparseDataMatrix& matrix, DataVector& norms);

void reaplyDiagonal(SparseDataMatrix& matrix, const DataVector& norms);

void reaplyDiagonalLowerTriangular(SparseDataMatrix& matrix, const DataVector& norms);

void normToUnitDiagonal(DataMatrix& matrix, DataVector& norms);

void reaplyDiagonal(DataMatrix& matrix, const DataVector& norms);

void reaplyDiagonalLowerTriangular(DataMatrix& matrix, const DataVector& norms);

} /* namespace IChol */
} /* namespace datadriven */
} /* namespace sgpp */
