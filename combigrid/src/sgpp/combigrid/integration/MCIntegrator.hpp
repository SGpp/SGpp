/*
 * MCIntegrator.hpp
 *
 *  Created on: 04.11.2015
 *      Author: david
 */

#ifndef MCINTEGRATOR_HPP_
#define MCINTEGRATOR_HPP_

#include <vector>
#include <functional>
#include <sgpp/globaldef.hpp>
#include <sgpp/base/datatypes/DataVector.hpp>

namespace sgpp{
namespace combigrid {

class MCIntegrator {
    std::function<base::DataVector(std::vector<base::DataVector> const &)> func;

public:
    MCIntegrator(std::function<double(base::DataVector const &)> func);
    MCIntegrator(std::function<base::DataVector(std::vector<base::DataVector> const &)> func);

    double integrate(std::vector<std::pair<double, double>> domain, size_t num_samples);

    double average(std::vector<std::pair<double, double>> domain, size_t num_samples);
};

}
} /* namespace sgpp*/

#endif /* MCINTEGRATOR_HPP_ */
