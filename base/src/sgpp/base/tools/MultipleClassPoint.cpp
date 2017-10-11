// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <sgpp/base/tools/MultipleClassPoint.hpp>

#include <sgpp/base/operation/hash/OperationEval.hpp>
#include <sgpp/base/operation/BaseOpFactory.hpp>
#include <vector>
#include <cmath>
#include <tuple>
#include <string>
#include <algorithm>

namespace sgpp {
namespace base {
MultipleClassPoint::MultipleClassPoint(base::HashGridPoint& gp,
                                    std::vector<base::Grid*> grids,
                                    std::vector<base::DataVector*> alphas) {
    // Init all classes
    for (size_t t =  0 ; t < grids.size() ; t++) {
        base::DataVector coords(grids.at(t)->getDimension());
        std::unique_ptr<base::OperationEval>
              opEval(op_factory::createOperationEval(*grids.at(t)));
        gp.getStandardCoordinates(coords);
        double eval = opEval->eval(*alphas.at(t), coords);
        std::tuple<double, size_t, bool> c1 { eval, t,
                    grids.at(t)->getStorage().isContaining(gp) };
        insertDensitySorted(&c1);
        classById.push_back(c1);
    }
}

size_t MultipleClassPoint::getDominateClass() const {
    return std::get<1>(classByDensity.at(0));
}

double MultipleClassPoint::getDensity(size_t classId) const {
    return std::get<0>(classById.at(classId));
}

void MultipleClassPoint::addNeighbor(size_t neighbor, size_t dim, bool isLeft) {
    std::tuple<size_t, size_t, bool> neigh { neighbor, dim, isLeft };
    if ( std::find(neighbors.begin(), neighbors.end(), neigh) == neighbors.end() ) {
        neighbors.push_back(neigh);
    }
}

std::vector<std::tuple<size_t, size_t, bool>> MultipleClassPoint::getNeighbors() const {
    return neighbors;
}

void MultipleClassPoint::addBorder(size_t dim, size_t level, bool isLeft) {
    std::tuple<size_t, size_t, bool> b { dim, level, isLeft };
    if ( std::find(borders.begin(), borders.end(), b) == borders.end() ) {
        borders.push_back(b);
    }
}

std::vector<std::tuple<size_t, size_t, bool>> MultipleClassPoint::getBorders() const {
    return borders;
}

double MultipleClassPoint::getBorderScore() const {
    double result = 0.0;
    for ( size_t i = 0 ; i < borders.size() ; i++ ) {
        result += 1 / pow(2.0, static_cast<double>(std::get<1>(borders.at(i))));
    }
    return result;
}

void MultipleClassPoint::resortClasses() {
    // Resort classByDensity vector
    struct ClassCompare {
        bool operator() (const std::tuple<double, size_t, bool> & t1,
                         const std::tuple<double, size_t, bool> & t2) {
            return std::get<0>(t1) > std::get<0>(t2);
        }
    };
    std::sort(classByDensity.begin(), classByDensity.end(), ClassCompare());
}

bool MultipleClassPoint::isClassSet(size_t classId) const {
    return std::get<2>(classById.at(classId));
}

std::vector<std::tuple<double, size_t, bool>> MultipleClassPoint::getTopClasses(
                double percent) const {
    std::vector<std::tuple<double, size_t, bool>> result;
    double minDenNeeded = (1.0 - percent) * std::get<0>(classByDensity.at(0));
    for ( size_t i = 0 ; i < classByDensity.size() &&
                std::get<0>(classByDensity.at(i)) > minDenNeeded ; i++ ) {
        result.push_back(classByDensity.at(i));
    }
    return result;
}

void MultipleClassPoint::insertDensitySorted(std::tuple<double, size_t, bool>* ins) {
    struct ClassCompare {
        bool operator() (const std::tuple<double, size_t, bool> & t1,
                         const std::tuple<double, size_t, bool> & t2) {
            return std::get<0>(t1) > std::get<0>(t2);
        }
    };
    std::vector<std::tuple<double, size_t, bool>>::iterator iter = std::lower_bound(
            classByDensity.begin(), classByDensity.end(), *ins, ClassCompare());
    classByDensity.insert(iter, *ins);
}

} /* namespace base */
} /* namespace sgpp */
