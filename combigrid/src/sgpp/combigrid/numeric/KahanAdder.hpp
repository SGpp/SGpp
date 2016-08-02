/*
 * KahanAdder.hpp
 *
 *  Created on: 05.11.2015
 *      Author: david
 */

#ifndef KAHANADDER_HPP_
#define KAHANADDER_HPP_
#include <sgpp/globaldef.hpp>

namespace sgpp{
namespace combigrid {

class KahanAdder {
    double sum = 0.0;
    double c = 0.0;

public:
    void add(double x) {
        // taken from Wikipedia, Kahan summation algorithm
        double y = x - c;
        double t = sum + y;
        c = (t - sum) - y;
        sum = t;
    }

    double value() const {
        return sum;
    }
};

} /* namespace combigrid */
} /* namespace sgpp*/

#endif /* KAHANADDER_HPP_ */
