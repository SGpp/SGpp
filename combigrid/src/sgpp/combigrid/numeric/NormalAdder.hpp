/*
 * NormalAdder.hpp
 *
 *  Created on: 05.11.2015
 *      Author: david
 */

#ifndef NORMALADDER_HPP_
#define NORMALADDER_HPP_
#include <sgpp/globaldef.hpp>

namespace sgpp{
namespace combigrid {

class NormalAdder {
    double sum = 0.0;
public:
    void add(double x) {
        sum += x;
    }

    double value() const {
        return sum;
    }
};

}
} /* namespace sgpp*/

#endif /* NORMALADDER_HPP_ */
