// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#ifndef MULIPLECLASSHASHGRIDPOINT_HPP
#define MULIPLECLASSHASHGRIDPOINT_HPP

#include <vector>
#include <set>
#include <tuple>
#include <array>
#include <sgpp/datadriven/tools/SimpleMultiClassGenerator.hpp>

namespace sgpp {
namespace datadriven {
class MultipleClassPoint{
public:
    MultipleClassPoint(int classes);
    MultipleClassPoint(int classes, int seq, sgpp::datadriven::SimpleMultiClassGenerator gen);
    virtual ~MultipleClassPoint();

    int getDominateClass();
    void updateClass(int classId, double newDen);
    double getDensity(int classId);

    void addNeighbor(int neighbor);
    std::vector<int> getNeighbors();

    std::vector<std::tuple<double, int, bool>> getTopClasses(double percent);

    std::string toString();

private:
    int classes;
    // tuple: density, classId, points exits in class
    std::vector<std::tuple<double, int, bool> *> classById;
    std::vector<std::tuple<double, int, bool>> classByDensity;
    std::vector<int> neighbors;

    void insertDensitySorted(std::tuple<double, int, bool>* density);
};


} /* namespace datadriven */
} /* namespace sgpp */

#endif /* MULIPLECLASSHASHGRIDPOINT_HPP */
