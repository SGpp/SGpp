#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE SGPPCombigridModule
#include <boost/test/unit_test.hpp>
#include <sgpp/globaldef.hpp>
#include <sgpp/combigrid/operation/CombigridOperation.hpp>
#include <sgpp/combigrid/operation/CombigridMultiOperation.hpp>
#include <sgpp/combigrid/storage/tree/CombigridTreeStorage.hpp>
#include <sgpp/combigrid/integration/MCIntegrator.hpp>
#include <sgpp/combigrid/utils/Stopwatch.hpp>
#include <sgpp/combigrid/utils/Utils.hpp>

#include <memory>
#include <cmath>
#include <iostream>

using namespace SGPP::combigrid;
using namespace SGPP::base;

SGPP::float_t testFunction(DataVector const &coordinates) {
	SGPP::float_t prod = 1.0;

	for (size_t i = 0; i < coordinates.getSize(); ++i) {
		SGPP::float_t x = coordinates[i];
		prod *= exp(-x * x / float_t((i + 1) * (i + 1)));
	}

	return prod;
}

SGPP::float_t testFunction2(DataVector const &coordinates) {
	SGPP::float_t prod = 1.0;

	for (size_t i = 0; i < coordinates.getSize(); ++i) {
		SGPP::float_t x = coordinates[i];
		prod *= cos(-x * x / float_t((i + 1) * (i + 1)));
	}

	return prod;
}

SGPP::float_t testFunction3(DataVector const &coordinates) {
	SGPP::float_t prod = 1.0;

	for (size_t i = 0; i < coordinates.getSize(); ++i) {
		SGPP::float_t x = coordinates[i];
		prod *= (x * x) / (1 + x * x);
	}

	return prod;
}

SGPP::float_t testFunction4(DataVector const &coordinates) {
	SGPP::float_t prod = 1.0;

	for (size_t i = 0; i < coordinates.getSize(); ++i) {
		SGPP::float_t x = coordinates[i];
		prod *= exp(-x * x);
	}

	return prod;
}

SGPP::float_t testFunction5(DataVector const &x) {
	std::vector<double> k;
	for (size_t i = 0; i < x.getSize(); ++i) {
		SGPP::float_t root = sqrt(float_t(i) + 6405.0); // 6400 ist Wurzel von 80, die naechste Quadratzahl ist 6561, also 161 mal Nachkommastellen ungleich 0
		SGPP::float_t decimals = root - floor(root);
		k.push_back(float_t(i % 2 == 0 ? 1 : -1) * 100 * decimals * tan(float_t(i)));
	}

	SGPP::float_t sum = 0.0;

	for (size_t i = 0; i < x.getSize(); ++i) {
		SGPP::float_t tmp = 1;
		for (size_t j = 0; j < k.size(); ++j) {
			tmp *= sin((-1) * k[j] * (1 - x[i]) * float_t(i)) - SGPP::combigrid::pow(cos(k[k.size() - 1 - i] * x[i] * k[j]), 3)
					+ sin((1 - x[i]) * (1 - x[i]) + x[x.getSize() - 1 - i]);
		}
		sum += tmp;
	}

	return sum;
}

SGPP::float_t testFunction6(DataVector const &x) {
	return 1;
}

SGPP::float_t testFunction7(DataVector const &x) {
	return x[0];
}

SGPP::float_t testFunctionAtan(DataVector const &x) {
	return atan(50 * (x[0] - .35)) + M_PI / 2 + 4 * pow(x[1], 2); // + exp(x[0] * x[1] - 1);
}

/*
void printCTResults(size_t d, size_t q) {
	const size_t samples = 10;
	auto ctInterpolator = CombigridOperation::createExpClenshawCurtisPolynomialInterpolation(d, testFunction);
	auto domain = std::vector<std::pair<SGPP::float_t, SGPP::float_t>>(d, std::pair<SGPP::float_t, SGPP::float_t>(0.0, 1.0));

	MCIntegrator integrator([&](DataVector const &x) -> SGPP::float_t {
		double diff = testFunction(x) - ctInterpolator->evaluate(q, x);
		return diff * diff;
	});

	std::cout << "d = " << d << ", q = " << q << ": " << std::sqrt(integrator.average(domain, samples)) << std::endl;
}
*/

void printDifferences(size_t d, std::shared_ptr<AbstractMultiStorage<FloatArrayVector>> storage) {
	auto it = storage->getStoredDataIterator();

	std::cout << "Differences: \n";
	while(it->isValid()) {
		std::cout << "Level (";
		for(size_t i = 0; i < d-1; ++i) {
			std::cout << it->indexAt(i);
			std::cout << ", ";
		}
		std::cout << it->indexAt(d-1) << "): ";

		std::cout << it->value().norm() << "\n";

		it->moveToNext();
	}
}

void printCTResults(size_t d, size_t q) {
	auto func = testFunctionAtan;
	const size_t samples = 100;
	auto ctInterpolator = CombigridMultiOperation::createLinearLejaPolynomialInterpolation(d, func);
	auto domain = std::vector<std::pair<SGPP::float_t, SGPP::float_t>>(d, std::pair<SGPP::float_t, SGPP::float_t>(0.0, 1.0));

	MCIntegrator integrator([&](std::vector<DataVector> const &params) -> DataVector {
		//auto result = ctInterpolator->evaluate(q, params);
		auto result = ctInterpolator->evaluateAdaptive(q*5, params);

		//printDifferences(d, ctInterpolator->getDifferences());

		for(size_t i = 0; i < params.size(); ++i) {
			SGPP::float_t diff = func(params[i]) - result[i];
			result[i] = diff * diff;
		}

		return result;
	});

	std::cout << "d = " << d << ", q = " << q << ": " << std::sqrt(integrator.average(domain, samples)) << std::endl;
}

BOOST_AUTO_TEST_CASE(testInterpolation) {
	for (size_t d = 2; d <= 2; ++d) {
		for (size_t w = 2; w <= 8; ++w) {
			Stopwatch stopwatch;
			printCTResults(d, w);
			stopwatch.log();
		}
	}

	auto quadrature = CombigridMultiOperation::createLinearLejaQuadrature(3, MultiFunction([](SGPP::base::DataVector const &x){return 1.0;}));
	std::vector<SGPP::base::DataVector> input(1, SGPP::base::DataVector(0));
	auto result = quadrature->evaluate(3, input);
	std::cout << "Quadrature result: " << result[0] << "\n";
}


