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
class MultipleClassPoint{
 public:
        explicit MultipleClassPoint(size_t classes);
        MultipleClassPoint(base::HashGridPoint& gp, std::vector<base::Grid*> grids,
                    std::vector<base::DataVector*> alphas);
        virtual ~MultipleClassPoint() {}

        size_t getDominateClass() const;
        double getDensity(size_t classId) const;
        void updateClass(size_t classId, double newDen, bool hasPoint);

        void addNeighbor(size_t neighbor, size_t dim, bool isLeft);
        std::vector<std::tuple<size_t, size_t, bool>> getNeighbors() const;

        void addBorder(size_t dim, size_t level, bool isLeft);
        double getBorderScore() const;
        std::vector<std::tuple<size_t, size_t, bool>> getBorders() const;

        void resortClasses();
        bool isClassSet(size_t classId) const;
        std::vector<std::tuple<double, size_t, bool>> getTopClasses(double percent) const;

        std::string toString();

 private:
    int classes;
    // tuple: density, classId, points exits in class
    std::vector<std::tuple<double, size_t, bool> *> classById;
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
