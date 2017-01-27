// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include "MultipleClassPoint.hpp"

#include <vector>
#include <cmath>
#include <string>
#include <sstream>
#include <iomanip>
#include <tuple>
#include <algorithm>

namespace sgpp {
namespace datadriven {
MultipleClassPoint::MultipleClassPoint(int classes) {
	// init all classes as 0.0
    for (int i = 0 ; i < classes ; i++) {
    	std::tuple<double, int, bool>* c1 = new std::tuple<double, int, bool> { 0.0 , i , false };
	    insertDensitySorted(c1);
	    classById.push_back(c1);
	}
}

MultipleClassPoint::MultipleClassPoint(int classes, int seq, 
                                       sgpp::datadriven::SimpleMultiClassGenerator gen) {
	// init all classes
    for (int i = 0 ; i < classes ; i++) {
    	std::tuple<double, int, bool>* c1 =
    	       new std::tuple<double, int, bool> { gen.getEval(i, seq) , i , true };
	    insertDensitySorted(c1);
	    classById.push_back(c1);
	}

}

MultipleClassPoint::~MultipleClassPoint() {
    // TODO Auto-generated destructor stub
}

int MultipleClassPoint::getDominateClass() {
    return std::get<1>(classByDensity.at(0));
}

double MultipleClassPoint::getDensity(int classId) {
    return std::get<0>(*classById.at(classId));
}

void MultipleClassPoint::updateClass(int classId, double newDen) {
    // TODO
}

void MultipleClassPoint::addNeighbor(int neighbor) {
    neighbors.push_back(neighbor);
}

std::vector<int> MultipleClassPoint::getNeighbors() {
    return neighbors;
}

std::string MultipleClassPoint::toString() {
    std::string s = "";
    for (unsigned int i = 0 ; i < classById.size() ; i++ ) {
        std::tuple<double, int, bool> tmp = *(classById.at(i));
        s += " - (" + std::to_string(std::get<0>(tmp)) + ",";
        s += std::to_string(std::get<1>(tmp)) + ")";
    }
    s += "\n";
    for (unsigned int i = 0 ; i < classByDensity.size() ; i++ ) {
        s += " - (" + std::to_string(std::get<0>(classByDensity.at(i))) + ",";
        s += std::to_string(std::get<1>(classByDensity.at(i))) + ")";
    }
    s += "\n";
    for (unsigned int i = 0 ; i < neighbors.size() ; i++ ) {
        size_t tmp = neighbors.at(i);
        s += " - " + std::to_string(tmp);
    }
    return s;
}

std::vector<std::tuple<double, int, bool>> MultipleClassPoint::getTopClasses(double percent) {
	std::vector<std::tuple<double, int, bool>> result;
	double minDenNeeded = (1.0 - percent) * std::get<0>(classByDensity.at(0));
	for (unsigned int i = 0 ; std::get<0>(classByDensity.at(i)) > minDenNeeded ; i++ ) {
        result.push_back(classByDensity.at(i));
    }
	return result;
	
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