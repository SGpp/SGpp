// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <sgpp/datadriven/tools/NearestNeighbors.hpp>
#include <iostream>
#include <cmath>

/**
 * \page example_nearestNeighborsTest_cpp Nearest Neighbors
 * This example calculates all feature-interactions that arise
 * from an image with 64 pixels, when one only considers pixels
 * whose \f$ L_2 \f$ distance is not larger than \f$ \sqrt{2} \f$.
 */

int main(void) {
    /**
     * First create the neighbors of all pixels.
    */
    const auto neigh = sgpp::datadriven::NearestNeighbors(8, 8);
    /**
     * Then create all arising interaction terms up to an order of 3.
    */
    const auto combs = neigh.getAllInteractions(3, std::sqrt(2));
    for (const auto& comb : combs) {
        for (const auto term : comb) {
            std::cout << term << ' ';
        }
        std::cout << std::endl;
    }
}
