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

namespace SGPP {
namespace combigrid {

class Stopwatch {
    std::chrono::time_point<std::chrono::high_resolution_clock> startTime;

public:
    Stopwatch();

    void start();

    SGPP::float_t elapsedSeconds();
    void log();
};

}
} /* namespace SGPP */

#endif /* STOPWATCH_HPP_ */
