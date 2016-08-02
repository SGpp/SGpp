/*
 * Stopwatch.cpp
 *
 *  Created on: 05.11.2015
 *      Author: david
 */

#include "../utils/Stopwatch.hpp"

#include <iostream>

namespace sgpp{
namespace combigrid {

Stopwatch::Stopwatch()
    : startTime(std::chrono::high_resolution_clock::now()) {
}

void Stopwatch::start() {
    startTime = std::chrono::high_resolution_clock::now();
}

double Stopwatch::elapsedSeconds() {
    std::chrono::duration<double> diff = std::chrono::high_resolution_clock::now() - startTime;
    return diff.count();
}

void Stopwatch::log() {
    std::cout << "Time: " << elapsedSeconds() << "s." << std::endl;
}

}
} /* namespace sgpp*/
