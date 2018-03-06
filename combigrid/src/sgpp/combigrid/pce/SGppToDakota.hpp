// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#pragma once

#include <sgpp/base/datatypes/DataVector.hpp>
#include <sgpp/base/datatypes/DataMatrix.hpp>
#include <sgpp/combigrid/algebraic/FloatTensorVector.hpp>

#ifdef USE_DAKOTA
#include <pecos_data_types.hpp>
#endif

namespace sgpp {
namespace combigrid {

#ifdef USE_DAKOTA
namespace sgppToDakota {

void convertDataVectorToRealVector(sgpp::base::DataVector& data_vec, Pecos::RealVector& real_vec);
void convertDataMatrixToRealMatrix(sgpp::base::DataMatrix& data_mat, Pecos::RealMatrix& real_mat);
void convertMultiIndexToUShortArray(MultiIndex& multiIndex, Pecos::UShortArray& ushort_vec);
void convertTensorToIndexArrayAndCoefficients(sgpp::combigrid::FloatTensorVector& tensorResult,
                                              Pecos::UShort2DArray& multiIndices,
                                              Pecos::RealVector& expansionCoefficients);
} /* namespace sgppToDakota */

#endif
} /* namespace combigrid */
} /* namespace sgpp */
