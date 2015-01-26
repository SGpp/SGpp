/* ****************************************************************************
* Copyright (C) 2010 Technische Universitaet Muenchen                         *
* This file is part of the SG++ project. For conditions of distribution and   *
* use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
**************************************************************************** */
// @author Richard Röttger

#include <sgpp/base/grid/generation/PrewaveletGridGenerator.hpp>
#include <sgpp/base/grid/GridStorage.hpp>

#include <sgpp/base/exception/generation_exception.hpp>

#include <sgpp/base/grid/generation/hashmap/HashCoarsening.hpp>
#include <sgpp/base/grid/generation/hashmap/HashRefinement.hpp>
#include <sgpp/base/grid/generation/hashmap/HashGenerator.hpp>

#include <iostream>
//
#include <sgpp/globaldef.hpp>


namespace SGPP {
  namespace base {


    PrewaveletGridGenerator::PrewaveletGridGenerator(GridStorage* storage,
        GridStorage* shadowstorage) :
      storage(storage), shadowstorage(shadowstorage) {
    }

    PrewaveletGridGenerator::~PrewaveletGridGenerator() {
    }

    void PrewaveletGridGenerator::regular(int level) {
      if (level < 0) {
        throw generation_exception("Grid level value is negative");
      }

      HashGenerator gen;
      gen.regular(this->storage, level);
    }

    void PrewaveletGridGenerator::cliques(int level, size_t clique_size) {
	  if (level < 0) {
		throw generation_exception("Grid level value is negative");
	  }

	  HashGenerator gen;
	  gen.cliques(this->storage, static_cast<HashGenerator::level_t>(level), clique_size);
	}

    void PrewaveletGridGenerator::full(int level) {
      if (level < 0) {
        throw generation_exception("Grid level value is negative");
      }

      HashGenerator gen;
      gen.full(this->storage, level);
    }

    /**
     * Refines the grid and updates the shadow storage.
     */
    void PrewaveletGridGenerator::refine(RefinementFunctor* func) {
      HashRefinement refine;
      size_t start = this->storage->size();
      refine.free_refine(this->storage, func);
      size_t end = this->storage->size();
      //All added gridpoint are between [start,end[

      //Check if a gridpoint within the shadow storage is now part of the actual grid!
      for (size_t i = start; i < end; i++) {
        if (shadowstorage->find(storage->get(i)) != shadowstorage->end()) {
          consolidateShadow();
          break;
        }
      }

      //Now add all missing neigbours to the shadowStorage
      for (size_t i = start; i < end; i++) {
        GridStorage::index_pointer index = this->storage->get(i);

        level_t sum = 0;

        for (size_t d = 0; d < storage->dim(); ++d) {
          index_t current_index;
          level_t current_level;
          index->get(d, current_level, current_index);
          sum += current_level;
        }

        GridStorage::grid_iterator iter(storage);
        GridStorage::grid_iterator shadowIter(shadowstorage);
        addNeighbours(*this->storage->get(i), 0, sum, iter, shadowIter);
      }

    }

    size_t PrewaveletGridGenerator::getNumberOfRefinablePoints() {
      HashRefinement refine;
      return refine.getNumberOfRefinablePoints(this->storage);
    }

    void PrewaveletGridGenerator::insertParents(GridStorage::grid_iterator& iter,
        GridStorage::grid_iterator& shadowIter) {
      // Call parents in every dimension
      for (size_t d = 0; d < storage->dim(); d++) {
        index_t current_index;
        level_t current_level;
        iter.get(d, current_level, current_index);

        if (current_level == 1)
          continue;

        iter.up(d);
        shadowIter.up(d);
        this->storage->get(iter.seq())->setLeaf(false);

        //Ok, point is neither in storage, nor in shadowstorage ...
        if (storage->end(iter.seq()) && shadowstorage->end(shadowIter.seq())) {
          GridStorage::index_pointer new_index = new GridStorage::index_type(
            storage->dim());

          for (size_t dim = 0; dim < storage->dim(); ++dim) {
            index_t target_index;
            level_t target_level;

            iter.get(dim, target_level, target_index);

            new_index->set(dim, target_level, target_index);

          }

          shadowstorage->insert(*new_index);
          delete new_index;
          insertParents(iter, shadowIter);
        }

        iter.set(d, current_level, current_index);
        shadowIter.set(d, current_level, current_index);
      }
    }

    void PrewaveletGridGenerator::addNeighbours(index_type& index,
        size_t current_dim, level_t target_level,
        GridStorage::grid_iterator& iter,
        GridStorage::grid_iterator& shadowIter) {

      level_t sum = 0;

      for (size_t d = 0; d < storage->dim(); ++d) {
        index_t current_index;
        level_t current_level;
        iter.get(d, current_level, current_index);
        sum += current_level;
      }

      if (sum == target_level) {
        GridStorage::index_pointer new_index = new GridStorage::index_type(
          storage->dim());

        if (storage->end(iter.seq()) && shadowstorage->end(shadowIter.seq())) {
          //Ok, point is neither in storage, nor in shadowstorage ...
          //check if the border of index and iter touching
          for (size_t d = 0; d < storage->dim(); ++d) {
            index_t target_index;
            level_t target_level;

            index_t current_index;
            level_t current_level;

            iter.get(d, current_level, current_index);
            index.get(d, target_level, target_index);

            new_index->set(d, current_level, current_index);

            // The index cast to int is required to allow a negative index
            int target_left = static_cast<int>((1.0 / (1 << target_level))
                                               * static_cast<double> (static_cast<int>(target_index) - 3));
            int target_right = static_cast<int>((1.0 / (1 << target_level))
                                                * static_cast<double> (static_cast<int>(target_index) + 3));
            int current_left = static_cast<int>((1.0 / (1 << current_index))
                                                * static_cast<double> (static_cast<int>(current_level) + 3));
            int current_right = static_cast<int>((1.0 / (1 << current_index))
                                                 * static_cast<double> (static_cast<int>(current_level) + 3));

            if (!(current_right > target_left || current_left
                  < target_right)) {
              return;
              delete new_index;
            }

          }

          //Yepp, the supports touching each other! Add point to shadow!
          shadowstorage->insert(*new_index);
          delete new_index;
          //Call for parents
          insertParents(iter, shadowIter);
        }

        return;
      } else if (sum > target_level) {
        return;
      }

      for (size_t d = current_dim; d < storage->dim(); d++) {
        index_t save_index;
        level_t save_level;
        iter.get(d, save_level, save_index); //Save current index

        iter.left_child(d);
        shadowIter.left_child(d);
        addNeighbours(index, d, target_level, iter, shadowIter);
        iter.set(d, save_level, save_index); //reset index
        shadowIter.set(d, save_level, save_index); //reset index

        iter.right_child(d);
        shadowIter.right_child(d);
        addNeighbours(index, d, target_level, iter, shadowIter);
        iter.set(d, save_level, save_index); //reset index
        shadowIter.set(d, save_level, save_index);
      }

    }

    /**
     * If during the refinement one or more points of the shadow register are added
     * to the actual grid then we have to remove these points from the shadow storage.
     */
    void PrewaveletGridGenerator::consolidateShadow() {

      GridStorage* temp = new GridStorage(storage->dim());

      for (size_t i = 0; i < shadowstorage->size(); i++) {
        temp->insert(*shadowstorage->get(i));
      }

      shadowstorage->emptyStorage();

      for (size_t i = 0; i < temp->size(); i++) {
        if (storage->find(temp->get(i)) == storage->end()) {
          shadowstorage->insert(*temp->get(i));
        }
      }

      delete temp;
    }

    void PrewaveletGridGenerator::coarsen(CoarseningFunctor* func,
                                          DataVector* alpha) {
      HashCoarsening coarsen;
      coarsen.free_coarsen(this->storage, func, alpha);
    }

    void PrewaveletGridGenerator::coarsenNFirstOnly(CoarseningFunctor* func,
        DataVector* alpha, size_t numFirstOnly) {
      HashCoarsening coarsen;
      coarsen.free_coarsen_NFirstOnly(this->storage, func, alpha, numFirstOnly);
    }

    size_t PrewaveletGridGenerator::getNumberOfRemovablePoints() {
      HashCoarsening coarsen;
      return coarsen.getNumberOfRemovablePoints(this->storage);
    }

    void PrewaveletGridGenerator::refineMaxLevel(RefinementFunctor* func,
        int maxLevel) {
      if (maxLevel < 0) {
        throw generation_exception("Grid level value is negative");
      }

      throw generation_exception("PrewaveletGridGenerator::refineMaxLevel is not implemented");
    }

    size_t PrewaveletGridGenerator::getNumberOfRefinablePointsToMaxLevel(
      int maxLevel) {
      if (maxLevel < 0) {
        throw generation_exception("Grid level value is negative");
      }

      throw generation_exception("PrewaveletGridGenerator::getNumberOfRefinablePointsToMaxLevel is not implemented");
      return 0;
    }
  }
}
