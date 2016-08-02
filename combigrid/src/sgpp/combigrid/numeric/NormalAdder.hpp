/*
 * NormalAdder.hpp
 *
 *  Created on: 05.11.2015
 *      Author: david
 */

#ifndef NORMALADDER_HPP_
#define NORMALADDER_HPP_
#include <sgpp/globaldef.hpp>

namespace SGPP {
namespace combigrid {

class NormalAdder {
    SGPP::float_t sum = 0.0;
public:
    void add(SGPP::float_t x) {
        sum += x;
    }

    SGPP::float_t value() const {
        return sum;
    }
};

}
} /* namespace SGPP */

#endif /* NORMALADDER_HPP_ */
