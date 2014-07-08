#include "PerformanceTest.hpp"
#include <iostream>

int main() {
	size_t dim = 3;
	int level = 3;
	int trainSize = 1000;
	int seed = 1208108;

	double maxIterations = 1000;

	double eps = 0.001;
	double lambda = 0.00001;
	double gamma = 0.1;

	sg::test::PerformanceTest p(true);


	p.runFriedman(dim, level, trainSize, seed, maxIterations, eps, lambda, gamma);

	std::cout << "Error vector:" << std::endl;
	std::cout << p.getError().toString() << std::endl;
	std::cout << std::endl;
}
