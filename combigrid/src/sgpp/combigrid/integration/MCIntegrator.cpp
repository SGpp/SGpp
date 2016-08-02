/*
 * MCIntegrator.cpp
 *
 *  Created on: 04.11.2015
 *      Author: david
 */

#include "MCIntegrator.hpp"

#include <chrono>
#include <random>

namespace sgpp{
namespace combigrid {

MCIntegrator::MCIntegrator(std::function<double(const base::DataVector&)> func) : func([=](const std::vector<base::DataVector> &vec) -> base::DataVector {
	base::DataVector result(vec.size());
	for(size_t i = 0; i < vec.size(); ++i) {
		result[i] = func(vec[i]);
	}
	return result;
}) {
}

MCIntegrator::MCIntegrator(std::function<base::DataVector(const std::vector<base::DataVector> &)> func) : func(func) {
}

double MCIntegrator::average(std::vector<std::pair<double, double> > domain, size_t num_samples) {
    //static std::default_random_engine generator(std::chrono::system_clock::now().time_since_epoch().count());
    static std::default_random_engine generator(15680329860); // TODO: deterministic

    double sum = 0.0;

    std::vector<base::DataVector> evaluationPoints;

    for(size_t i = 0; i < num_samples; ++i) {
        base::DataVector coordinates(domain.size());
        for(size_t dim = 0; dim < domain.size(); ++dim) {
            std::uniform_real_distribution<double> distribution(domain[dim].first, domain[dim].second);
            coordinates[dim] = distribution(generator);
        }

        evaluationPoints.push_back(coordinates);
    }

    base::DataVector results = func(evaluationPoints);

    for(size_t i = 0; i < results.getSize(); ++i) {
    	sum += results[i];
    }

    return sum / static_cast<double>(num_samples);
}

double MCIntegrator::integrate(std::vector<std::pair<double, double> > domain, size_t num_samples) {
    double volume = 1.0;

    for(auto bounds : domain) {
        volume *= (bounds.second - bounds.first);
    }

    return volume * average(domain, num_samples);
}

}
} /* namespace sgpp*/
