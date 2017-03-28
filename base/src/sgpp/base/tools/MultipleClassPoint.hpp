// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#ifndef MULTIPLECLASSPOINT_HPP
#define MULTIPLECLASSPOINT_HPP

#include <vector>
#include <tuple>
#include <array>
#include <sgpp/base/grid/Grid.hpp>

namespace sgpp {
namespace base {
class MultipleClassPoint{
    public:
        explicit MultipleClassPoint(int classes);
        MultipleClassPoint(base::HashGridPoint& gp, std::vector<base::Grid*> grids,
                    std::vector<base::DataVector*> alphas);
        virtual ~MultipleClassPoint() {};

        int getDominateClass() const;
        double getDensity(int classId) const;
        void updateClass(int classId, double newDen, bool hasPoint);

        void addNeighbor(int neighbor, size_t dim, bool isLeft);
        std::vector<std::tuple<int, size_t, bool>> getNeighbors();

        void resortClasses();
        bool isClassSet(int classId);
        std::vector<std::tuple<double, int, bool>> getTopClasses(double percent);

        std::string toString();

    private:
        int classes;
        // tuple: density, classId, points exits in class
        std::vector<std::tuple<double, int, bool> *> classById;
        std::vector<std::tuple<double, int, bool>> classByDensity;
        // sequence number and dimension of neighbors with a change in dominate classes
        std::vector<std::tuple<int, size_t, bool>> neighbors;

        void insertDensitySorted(std::tuple<double, int, bool>* density);
};


} /* namespace base */
} /* namespace sgpp */

#endif /* MULTIPLECLASSPOINT_HPP */
