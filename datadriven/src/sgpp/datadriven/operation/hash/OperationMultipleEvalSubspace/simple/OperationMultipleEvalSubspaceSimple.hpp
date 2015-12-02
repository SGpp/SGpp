// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#pragma once

#include <iostream>
#include <vector>
#include <map>

//////////////////////////////////////////////////////////////////////
// Caution: Subspace-skipping is disabled by default for this kernel
//////////////////////////////////////////////////////////////////////

#include <omp.h>

#include <sgpp/datadriven/operation/hash/OperationMultipleEvalSubspace/AbstractOperationMultipleEvalSubspace.hpp>
#include <sgpp/datadriven/operation/hash/OperationMultipleEvalSubspace/simple/SubspaceNodeSimple.hpp>
#include <sgpp/datadriven/operation/hash/OperationMultipleEvalSubspace/simple/OperationMultipleEvalSubspaceSimpleParameters.hpp>

#include <sgpp/globaldef.hpp>


namespace SGPP {
  namespace datadriven {

    class OperationMultipleEvalSubspaceSimple: public AbstractOperationMultipleEvalSubspace {
      private:

        size_t dim = -1;
        size_t maxLevel = 0;

        size_t* allSubspaces = NULL;
        size_t subspaceCount = -1;
        size_t subspaceSize = -1;

        size_t totalGridPoints = 0;
        float_t* allSurplusses = nullptr;
        std::map<uint32_t, uint32_t> allSurplussesIndexMap;

        void prepareSubspaceIterator();

        void createFlatStorage();

        void setSurplus(std::vector<size_t>&  level, std::vector<size_t>& maxIndices, std::vector<size_t>& index, float_t value);

        void getSurplus(std::vector<size_t>& level, std::vector<size_t>& maxIndices, std::vector<size_t>& index, float_t& value, bool& isVirtual);

        void setCoefficients(base::DataVector& surplusVector);

        void unflatten(base::DataVector& result);

        size_t flattenIndex(size_t dim, std::vector<size_t>& maxIndices, std::vector<size_t>& index);

        size_t flattenIndex(size_t* intermediates, size_t dim, size_t* maxIndicesPtr, size_t* indexPtr,
                            size_t toRecalc);

        size_t flattenLevel(size_t dim, size_t maxLevel, std::vector<size_t>& level);

        static inline size_t calculateIndexComponent(size_t dim, float_t unadjusted) {
          //implies flooring
          size_t rounded = static_cast<size_t>(unadjusted);

          size_t mask = 0x1;
          size_t sign = mask ^ (mask & rounded);

          size_t componentIndex = rounded + sign;
          return componentIndex;
        }

      public:

        OperationMultipleEvalSubspaceSimple(base::Grid& grid, base::DataMatrix& dataset);

        ~OperationMultipleEvalSubspaceSimple();

        void prepare() override;

        void multTransposeImpl(SGPP::base::DataVector& alpha, SGPP::base::DataVector& result, const size_t start_index_data,
                               const size_t end_index_data) override;

        void multImpl(SGPP::base::DataVector& source, SGPP::base::DataVector& result, const size_t start_index_data,
                      const size_t end_index_data) override;

        size_t getAlignment() override;

        std::string getImplementationName() override;
    };
  }
}
