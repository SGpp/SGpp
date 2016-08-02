/*
 * Stopwatch.hpp
 *
 *  Created on: 05.11.2015
 *      Author: david
 */

#ifndef STOPWATCH_HPP_
#define STOPWATCH_HPP_

#include <chrono>
#include <sgpp/globaldef.hpp>

namespace sgpp{
namespace combigrid {

class Stopwatch {
    std::chrono::time_point<std::chrono::high_resolution_clock> startTime;

public:
    Stopwatch();

    void start();

    double elapsedSeconds();
    void log();
};

}
} /* namespace sgpp*/

#endif /* STOPWATCH_HPP_ */
