// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <sgpp/datadriven/tools/NearestNeighbors.hpp>
#include <iostream>
#include <cmath>

/**
 * @brief main prints out all interaction terms for a 64-dimensional
 * image dataset that are in a \f$ 3 \times 3 \f$ neighborhood.
 */
int main(void) {
    const auto neigh = sgpp::datadriven::NearestNeighbors(8, 8);
    const auto combs = neigh.getAllInteractions(3, std::sqrt(2));
    for (const auto& comb : combs) {
        for (const auto term : comb) {
            std::cout << term << ' ';
        }
        std::cout << std::endl;
    }
}
