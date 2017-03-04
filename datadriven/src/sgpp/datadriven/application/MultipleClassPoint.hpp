// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#ifndef MULIPLECLASSHASHGRIDPOINT_HPP
#define MULIPLECLASSHASHGRIDPOINT_HPP

#include <vector>
#include <unordered_set>
#include <tuple>
#include <array>
#include <sgpp/datadriven/tools/SimpleMultiClassGenerator.hpp>

namespace sgpp {
namespace datadriven {
class MultipleClassPoint{
public:
    MultipleClassPoint(int classes);
    MultipleClassPoint(base::HashGridPoint& gp, std::vector<base::Grid*> grids,
            std::vector<base::DataVector*> alphas);
    virtual ~MultipleClassPoint();

    int getDominateClass() const;
    void updateClass(int classId, double newDen, bool hasPoint);
    double getDensity(int classId) const;

    void addNeighbor(int neighbor, int dim);
    std::vector<std::tuple<int, int>> getNeighbors();

    void resortClasses();

    std::vector<std::tuple<double, int, bool>> getTopClasses(double percent);

    std::string toString();

private:
    int classes;
    // tuple: density, classId, points exits in class
    std::vector<std::tuple<double, int, bool> *> classById;
    std::vector<std::tuple<double, int, bool>> classByDensity;
    // sequence number and dimension of neighbors with a change in dominate classes
    std::vector<std::tuple<int, int>> neighbors;

    void insertDensitySorted(std::tuple<double, int, bool>* density);
};


} /* namespace datadriven */
} /* namespace sgpp */

#endif /* MULIPLECLASSHASHGRIDPOINT_HPP */
