// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <sgpp/combigrid/pce/SGppToDakota.hpp>

#ifdef USE_DAKOTA
namespace sgpp {
namespace combigrid {
namespace sgppToDakota {

void convertDataVectorToRealVector(sgpp::base::DataVector& data_vec, Pecos::RealVector& real_vec) {
  for (size_t i = 0; i < data_vec.getSize(); i++) {
    real_vec[static_cast<int>(i)] = data_vec[i];
  }
}

void convertDataMatrixToRealMatrix(sgpp::base::DataMatrix& data_mat, Pecos::RealMatrix& real_mat) {
  for (size_t i = 0; i < data_mat.getNcols(); i++) {
    for (size_t j = 0; j < data_mat.getNrows(); j++) {
      real_mat(static_cast<int>(i), static_cast<int>(j)) = data_mat.get(i, j);
    }
  }
}

void convertTensorToIndexArrayAndCoefficients(sgpp::combigrid::FloatTensorVector& tensorResult,
                                              Pecos::UShort2DArray& multiIndices,
                                              Pecos::RealVector& expansionCoefficients) {
  auto it = tensorResult.getValues()->getStoredDataIterator();
  if (it->isValid()) {
    size_t numDims = tensorResult.getValues()->getNumDimensions();
    for (int n = 0; it->isValid(); it->moveToNext(), n++) {
      auto multiIndex = it->getMultiIndex();
      Pecos::UShortArray multiIndex_ushort(numDims);
      convertMultiIndexToUShortArray(multiIndex, multiIndex_ushort);
      multiIndices.push_back(multiIndex_ushort);
      expansionCoefficients.resize(n + 1);
      expansionCoefficients[n] = it->value().value();
    }
  }
}

void convertMultiIndexToUShortArray(MultiIndex& multiIndex, Pecos::UShortArray& ushort_vec) {
  for (size_t idim = 0; idim < multiIndex.size(); idim++) {
    ushort_vec[idim] = multiIndex[idim];
  }
}

} /* namespace sgppToDakota */
} /* namespace combigrid */
} /* namespace sgpp */

#endif
