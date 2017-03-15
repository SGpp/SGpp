// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <sgpp/base/operation/hash/OperationEval.hpp>
#include <sgpp/base/operation/BaseOpFactory.hpp>
#include "MultipleClassPoint.hpp"

#include <vector>
#include <unordered_set>
#include <cmath>
#include <sstream>
#include <iomanip>
#include <tuple>
#include <string>
#include <algorithm>

namespace sgpp {
namespace datadriven {
MultipleClassPoint::MultipleClassPoint(int classes) {
    // init all classes as 0.0
    for (int i = 0 ; i < classes ; i++) {
        std::tuple<double, int, bool>* c1 =
                new std::tuple<double, int, bool> { 0.0 , i , false };
        insertDensitySorted(c1);
        classById.push_back(c1);
    }
}

MultipleClassPoint::MultipleClassPoint(base::HashGridPoint& gp,
                                    std::vector<base::Grid*> grids,
                                    std::vector<base::DataVector*> alphas) {
    // init all classes
    for (size_t t =  0 ; t < grids.size() ; t++) {
        base::DataVector coords(grids.at(t)->getDimension());
        std::unique_ptr<base::OperationEval>
              opEval(op_factory::createOperationEval(*grids.at(t)));
        gp.getStandardCoordinates(coords);
        double eval = opEval->eval(*alphas.at(t), coords);
        std::tuple<double, int, bool>* c1 =
               new std::tuple<double, int, bool> { eval, t,
                    grids.at(t)->getStorage().isContaining(gp) };
        insertDensitySorted(c1);
        classById.push_back(c1);
    }
}
/*
MultipleClassPoint::~MultipleClassPoint() {
    // TODO(degelkn): Auto-generated destructor stub
}
*/
int MultipleClassPoint::getDominateClass() const {
    return std::get<1>(classByDensity.at(0));
}

double MultipleClassPoint::getDensity(int classId) const {
    return std::get<0>(*classById.at(classId));
}

void MultipleClassPoint::updateClass(int classId, double newDen, bool hasPoint) {
    std::tuple<double, int, bool>* c1 =
           new std::tuple<double, int, bool> { newDen , classId , hasPoint };
    struct ClassCompare {
     std::tuple<double, int, bool> oldClass;
     public: ClassCompare(std::tuple<double, int, bool> i):oldClass(i) {}
     bool operator()(const std::tuple<double, int, bool>& t1) {
       return std::get<1>(t1) == std::get<1>(oldClass);
     }
    };
    std::replace_if(classByDensity.begin(),
            classByDensity.end(), ClassCompare(*c1), *c1);
    classById.at(classId) = c1;
}

void MultipleClassPoint::addNeighbor(int neighbor, size_t dim, bool isLeft) {
    std::tuple<int, size_t, bool> neigh = { neighbor, dim, isLeft };
    if ( std::find(neighbors.begin(), neighbors.end(), neigh) == neighbors.end() ) {
        neighbors.push_back(neigh);
    }
}

std::vector<std::tuple<int, size_t, bool>> MultipleClassPoint::getNeighbors() {
    return neighbors;
}

void MultipleClassPoint::resortClasses() {
    // resort vector
    struct ClassCompare {
        bool operator() (const std::tuple<double, int, bool> & t1,
                         const std::tuple<double, int, bool> & t2) {
            return std::get<0>(t1) > std::get<0>(t2);
        }
    };
    std::sort(classByDensity.begin(), classByDensity.end(), ClassCompare());
}

bool MultipleClassPoint::isClassSet(int classId) {
    return std::get<2>(*(classById.at(classId)));
}

std::vector<std::tuple<double, int, bool>> MultipleClassPoint::getTopClasses(double percent) {
    std::vector<std::tuple<double, int, bool>> result;
    double minDenNeeded = (1.0 - percent) * std::get<0>(classByDensity.at(0));
    for ( unsigned int i = 0 ; i < classByDensity.size() &&
                std::get<0>(classByDensity.at(i)) > minDenNeeded ; i++ ) {
        result.push_back(classByDensity.at(i));
    }
    return result;
}

std::string MultipleClassPoint::toString() {
    std::string s = "-> ";
    for ( unsigned int i = 0 ; i < neighbors.size() ; i++ ) {
        s += std::to_string(std::get<0>(neighbors.at(i)))+ ", ";
    }
    s += "\n";
    for ( unsigned int i = 0 ; i < classById.size() ; i++ ) {
        std::tuple<double, int, bool> tmp = *(classById.at(i));
        s += " - (" + std::to_string(std::get<0>(tmp)) + ",";
        s += std::to_string(std::get<1>(tmp)) + ",";
        s += std::to_string(std::get<2>(tmp)) + ")";
    }
    s += "\n";
    for ( unsigned int i = 0 ; i < classByDensity.size() ; i++ ) {
        s += " - (" + std::to_string(std::get<0>(classByDensity.at(i))) + ",";
        s += std::to_string(std::get<1>(classByDensity.at(i))) + ")";
    }
    return s;
}

void MultipleClassPoint::insertDensitySorted(std::tuple<double, int, bool>* ins) {
    struct ClassCompare {
        bool operator() (const std::tuple<double, int, bool> & t1,
                         const std::tuple<double, int, bool> & t2) {
            return std::get<0>(t1) > std::get<0>(t2);
        }
    };
    std::vector<std::tuple<double, int, bool>>::iterator iter = std::lower_bound(
            classByDensity.begin(), classByDensity.end(), *ins, ClassCompare());
    classByDensity.insert(iter, *ins);
}

} /* namespace datadriven */
} /* namespace sgpp */
