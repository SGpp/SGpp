/*
 * KahanAdder.hpp
 *
 *  Created on: 05.11.2015
 *      Author: david
 */

#ifndef KAHANADDER_HPP_
#define KAHANADDER_HPP_
#include <sgpp/globaldef.hpp>

namespace SGPP {
namespace combigrid {

class KahanAdder {
    SGPP::float_t sum = 0.0;
    SGPP::float_t c = 0.0;

public:
    void add(SGPP::float_t x) {
        // taken from Wikipedia, Kahan summation algorithm
        SGPP::float_t y = x - c;
        SGPP::float_t t = sum + y;
        c = (t - sum) - y;
        sum = t;
    }

    SGPP::float_t value() const {
        return sum;
    }
};

} /* namespace combigrid */
} /* namespace SGPP */

#endif /* KAHANADDER_HPP_ */
