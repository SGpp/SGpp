// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#pragma once

#include "sgpp/base/grid/GridStorage.hpp"
#include "sgpp/base/grid/storage/hashmap/HashGridIndex.hpp"
#include "sgpp/base/grid/storage/hashmap/HashGridIterator.hpp"
#include "sgpp/base/grid/storage/hashmap/SerializationVersion.hpp"

#include "sgpp/base/grid/common/BoundingBox.hpp"
#include "sgpp/base/grid/common/Stretching.hpp"

#include "sgpp/base/datatypes/DataMatrix.hpp"
#include "sgpp/base/datatypes/DataMatrixSP.hpp"

#include "sgpp/parallel/tools/TypesParallel.hpp"
#include "sgpp/parallel/datadriven/tools/DMVectorizationPaddingAssistant.hpp"

#include <memory>
#include <vector>
#include <string>
#include <sstream>
#include <exception>
#include <list>
#include <typeinfo>
#include <stdint.h>


namespace SGPP{
namespace parallel{
/**
 * Generic hash table based index storage converter.
 * used to convert HashGridStorage to AOS (array of structures type)
 */
class HashGridStorageConverter{
public:
	/**
	 * Converts this storage from AOS (array of structures) to SOA (structure of array)
	 * with modification to speed up iterative function evaluation. The Level
	 * array won't contain the levels, it contains the level to the power of two
	*
	* The generated arrays are made in format optimized for minimizing page faults
	 *
	 * @param storage GridStorage, which should be converted
	 * @param level DataMatrix to store the grid's level to the power of two
	 * @param index DataMatrix to store the grid's indices
	 * @param vectorizationType Vectorization type
	 * @param blocking_length parameter for an additional blocking length to avoid TLB misses
	 */
	static void getLevelIndexArraysForEvalTLBOptimized(SGPP::base::GridStorage* storage,
			SGPP::base::DataMatrix& level, SGPP::base::DataMatrix& index,
			sg::parallel::VectorizationType vectorizationType, size_t blocking_length) {
	  typename SGPP::base::HashGridStorage::index_type::level_type curLevel;
	  typename SGPP::base::HashGridStorage::index_type::level_type curIndex;

	  //pad datasets
	  sg::parallel::DMVectorizationPaddingAssistant::padDataset(level, vectorizationType);
	  sg::parallel::DMVectorizationPaddingAssistant::padDataset(index, vectorizationType);

	  level.setAll(0.0);
	  index.setAll(0.0);

	  //transpose
	  level.transpose();
	  index.transpose();

	  //make optimized for reducing page faults

	  double* level_ptr = level.getPointer();
	  double* index_ptr = index.getPointer();

	  for (size_t i = 0; i < storage->getSize(); i += blocking_length) {
		for (size_t current_dim = 0; current_dim < storage->getDimension(); current_dim++) {
		  for (size_t t = i; t < i + blocking_length; ++t) {
			if (t < storage->getSize()) {
			  (*storage)[t]->get(current_dim, curLevel, curIndex);
			  *level_ptr = static_cast<double>(1 << curLevel);
			  *index_ptr = static_cast<double>(curIndex);
			}

			++level_ptr;
			++index_ptr;
		  }
		}
	  }
	}

    /**
     * Converts the storage from AOS (array of structures) to SOA (structure of array)
     * with modification to speed up iterative Laplace Calculations: the level
     * won't contain the levels, it contains 2 to the neagative power of the level.
     * Additional blocking for better TLB usage is provided.
     *
     * @param storage GridStorage, which should be converted
     * @param level DataMatrix to store the grid's modified level
     * @param vectorizationType Vectorization type
     * @param blocking_length parameter for an additional blocking length to avoid TLB misses
     */
	static void getLevelForIntegralTLBOptimized(SGPP::base::GridStorage* storage,
    		SGPP::base::DataMatrix& level,
    		sg::parallel::VectorizationType vectorizationType,
			size_t blocking_length) {
      typename SGPP::base::HashGridStorage::index_type::level_type curLevel;
      typename SGPP::base::HashGridStorage::index_type::level_type curIndex;

      //pad datasets
      sg::parallel::DMVectorizationPaddingAssistant::padDataset(level, vectorizationType);

      level.setAll(0.0);

      //transpose
      level.transpose();

      //make optimized for reducing page faults

      double* level_ptr = level.getPointer();

      for (size_t i = 0; i < storage->getSize(); i += blocking_length) {
        for (size_t current_dim = 0; current_dim < storage->getDimension(); current_dim++) {
          for (size_t t = i; t < i + blocking_length; ++t) {
            if (t < storage->getSize()) {
              (*storage)[t]->get(current_dim, curLevel, curIndex);
              *level_ptr = pow(2.0, static_cast<int>(-curLevel));
            }

            ++level_ptr;
          }
        }
      }
    }
};

}  // namespace parallel
}  // namespace SGPP
