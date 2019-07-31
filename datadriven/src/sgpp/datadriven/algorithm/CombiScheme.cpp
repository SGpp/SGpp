/*
 * Copyright (C) 2008-today The SG++ project
 * This file is part of the SG++ project. For conditions of distribution and
 * use, please see the copyright notice provided with SG++ or at
 * sgpp.sparsegrids.org
 *
 * CombiScheme.cpp
 *
 * Created on: Jul 25, 2019
 *     Author: Kilian Röhner
 */

#include <sgpp/datadriven/algorithm/CombiScheme.hpp>

#include <algorithm>
#include <cstdlib>
#include <functional>
#include <map>
#include <set>
#include <numeric>
#include <unordered_set>
#include <utility>
#include <vector>

namespace sgpp {
namespace datadriven {

void CombiScheme::initialize(size_t dim, size_t level) {
    dimension = dim;
    this->level = level;
    init_index_set();
}

bool CombiScheme::isRefinable(std::vector<size_t> levelvec) {
    if (index_set.find(levelvec) != index_set.end()) {
        return index_set[levelvec];
    }
    return false;
}

bool CombiScheme::refineComponent(std::vector<size_t> levelvec) {
    if (!isRefinable(levelvec)) {
        return false;
    }
    index_set[levelvec] = false;
    for (size_t dim = 0; dim < dimension; dim++) {
        refine_scheme(dim, levelvec);
    }
    return true;
}

void CombiScheme::refine_scheme(size_t dim, std::vector<size_t> levelvec) {
    levelvec[dim]++;
    for (size_t d = 0; d < dimension; d++) {
        std::vector<size_t> levelvec_copy = levelvec;
        levelvec_copy[d]--;
        // check if parents are all there
        if (index_set.find(levelvec_copy) == index_set.end() && levelvec_copy[d] > 0) {
            return;
        }
    }
    index_set[levelvec] = true;
}

void CombiScheme::init_index_set() {
    for (const std::vector<size_t>& v : getGrids(dimension, level)) {
        index_set[v] = true;
    }
    for (size_t q = 1; q < std::min(dimension, level); q++) {
        for (const std::vector<size_t>& v : getGrids(dimension, level - q)) {
            index_set[v] = false;
        }
    }
}

std::unordered_set<std::vector<size_t>, VectorHash> CombiScheme::getGrids(
    size_t dim, size_t values) {
    std::unordered_set<std::vector<size_t>, VectorHash> grids;
    if (dim == 1) {
         grids.insert(std::vector<size_t>(1, values));
         return grids;
    }
    for (size_t index = 0; index < values; index++) {
        for (std::vector<size_t> v : getGrids(dim - 1, values - index)) {
            v.insert(v.begin(), index + 1);
            grids.insert(v);
        }
    }
    return grids;
}

std::vector<std::pair<std::vector<size_t>, int>> CombiScheme::getCombiScheme() {
    std::vector<std::pair<std::vector<size_t>, int>> grid_array;
    std::map<std::vector<size_t>, int> grid_dict;

    for (auto const& index : index_set) {
        std::vector<std::vector<int>> stencil_elements;
        for (size_t d = 0; d < dimension; d++) {
            for (std::vector<int>& stencil_element : stencil_elements) {
                stencil_element.emplace_back(0);
            }
            if (stencil_elements.empty()) {
                std::vector<int> first{0};
                stencil_elements.push_back(first);
            }
            if (index.first[d] > 1) {
                // hier verdoppeln, zweite hälfte letztes element -1 statt 0
                std::vector<std::vector<int>> new_stencil_elements;
                for (std::vector<int>& stencil_element : stencil_elements) {
                    std::vector<int> new_stencil_element = stencil_element;
                    new_stencil_element.back() = -1;
                    new_stencil_elements.push_back(new_stencil_element);
                }
                for (std::vector<int>& stencil_element : new_stencil_elements) {
                    stencil_elements.push_back(stencil_element);
                }
            }
        }
        for (std::vector<int>& stencil_element : stencil_elements) {
            int update_coefficient = std::accumulate(stencil_element.begin(),
                                                     stencil_element.end(), 0);
            update_coefficient = (abs(update_coefficient - 1) % 2) - (abs(update_coefficient) % 2);
            std::vector<size_t> new_stencil_element(dimension, 0);
            std::transform(stencil_element.begin(), stencil_element.end(),
                           index.first.begin(), new_stencil_element.begin(), std::plus<int>());
            grid_dict[new_stencil_element] += update_coefficient;
        }
    }

    for (auto it = grid_dict.cbegin(); it != grid_dict.cend() ; it++) {
        if (it->second != 0) {
            grid_array.emplace_back(it->first, it->second);
        }
    }
    return grid_array;
}

} /* namespace datadriven */
} /* namespace sgpp */
