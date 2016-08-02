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

namespace SGPP {
namespace combigrid {

class MCIntegrator {
    std::function<base::DataVector(std::vector<base::DataVector> const &)> func;

public:
    MCIntegrator(std::function<SGPP::float_t(base::DataVector const &)> func);
    MCIntegrator(std::function<base::DataVector(std::vector<base::DataVector> const &)> func);

    SGPP::float_t integrate(std::vector<std::pair<SGPP::float_t, SGPP::float_t>> domain, size_t num_samples);

    SGPP::float_t average(std::vector<std::pair<SGPP::float_t, SGPP::float_t>> domain, size_t num_samples);
};

}
} /* namespace SGPP */

#endif /* MCINTEGRATOR_HPP_ */
