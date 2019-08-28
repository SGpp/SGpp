// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#ifndef MULTIPLECLASSPOINT_HPP
#define MULTIPLECLASSPOINT_HPP

#include <sgpp/base/grid/Grid.hpp>

#include <vector>
#include <tuple>
#include <array>
#include <string>

namespace sgpp {
namespace base {
/**
 * Multiple Class Point provides a structure to save additional data
 * for each GridPoint in a Sparse Grid.
 *
 * The Multiple Class Point is used by the MultipleClassRefinementFunctor
 * to provide the needed data.
 */
class MultipleClassPoint{
 public:
        /**
        * Constructor.
        *
        * @param gp GridPoint the MultipleClassPoint is associated with
        * @param grids Vector of grids
        * @param alphas Vector of surpluses related to the grids
        */
        MultipleClassPoint(base::HashGridPoint& gp, std::vector<base::Grid*> grids,
                    std::vector<base::DataVector*> alphas);
        virtual ~MultipleClassPoint() {}

        /**
         * Get the index of the class with the highest denstity at this point
         * @return the index of the doinating class
         */
        size_t getDominateClass() const;
        /**
         * Get the density of the class at this point
         * @param classId the index of the class
         * @return the density of the class with given index
         */
        double getDensity(size_t classId) const;
        /**
         * Adds the information for a found neighbor to the point
         * @param neighbor the index of the neighbor
         * @param dim the dimensions of the neighbor
         * @param isLeft the direction of the neighbor
         */
        void addNeighbor(size_t neighbor, size_t dim, bool isLeft);
        /**
         * Gets the list of all Neighbors
         * @return The list of all neigbors
         */
        std::vector<std::tuple<size_t, size_t, bool>> getNeighbors() const;
        /**
         * Adds the information for a found boundary to the point
         * @param dim the dimensions of the boundary
         * @param level the level of the point in given dimension
         * @param isLeft the direction of the boundary
         */
        void addBorder(size_t dim, size_t level, bool isLeft);
        /**
         * Gets the list of all Boundaries
         * @return The list of all boundaries
         */
        std::vector<std::tuple<size_t, size_t, bool>> getBorders() const;
        /**
         * Scores the boundaries dependent on the level.
         * @return A score based on the distance to the boundaries of the domain
         */
        double getBorderScore() const;
        /**
         * Sorts the classes by the densities.
         */
        void resortClasses();
        /**
         * Checks if a point exits in the sparse grid of the given class
         * @param classId the index of the class to check for
         * @return true, if point exists in grid with index classId
         */
        bool isClassSet(size_t classId) const;
        /**
         * Gets all classes with a density in the range of
         * percent to the dominating class
         * @param percent the range in which to return the classes
         * @return list of calees with their densities
         */
        std::vector<std::tuple<double, size_t, bool>> getTopClasses(double percent) const;

 private:
    // tuple: density, classId, points exits in class
    std::vector<std::tuple<double, size_t, bool>> classById;
    std::vector<std::tuple<double, size_t, bool>> classByDensity;
    // sequence number and dimension of neighbors with a change in dominate classes, isLeft
    std::vector<std::tuple<size_t, size_t, bool>> neighbors;
    // dimension, level of point, direction isLeft
    std::vector<std::tuple<size_t, size_t, bool>> borders;

    void insertDensitySorted(std::tuple<double, size_t, bool>* density);
};


} /* namespace base */
} /* namespace sgpp */

#endif /* MULTIPLECLASSPOINT_HPP */
