#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE SGppQuadratureModule
#include <boost/test/unit_test.hpp>

#include <sgpp_quadrature.hpp>
#include <sgpp_base.hpp>
#include <sgpp/globaldef.hpp>

using namespace SGPP::base;
using namespace SGPP::quadrature;

SGPP::float_t f(DataVector x) {
    SGPP::float_t res = 1.0f;
    for (size_t i = 0; i < x.getSize(); i++) {
        res *= 4 * (1 - x[i]) * x[i];
    }
    return res;
}

void testSampler(SampleGenerator& sampler, size_t dim, size_t numSamples,
        SGPP::float_t analyticResult, double tol) {
    DataVector sample = DataVector(dim);
    SGPP::float_t sum = 0.0f;

    for (size_t i = 0; i < numSamples; i++) {
        sampler.getSample(sample);
        sum += f(sample);
    }
    sum /= static_cast<float_t>(numSamples);
	BOOST_CHECK_CLOSE(sum, analyticResult, tol * 1e2);
}

BOOST_AUTO_TEST_CASE(testSamplers) {
    size_t dim = 2;
    size_t numSamples = 100000;
    int seed = 1234567;
    SGPP::float_t analyticResult = std::pow(2. / 3., dim);

    NaiveSampleGenerator pNSampler(dim, seed);
    HaltonSampleGenerator pHSampler(dim);
    LatinHypercubeSampleGenerator pLHSampler(dim, numSamples, seed);
    std::vector<size_t> blockSize(dim);
    for (size_t i = 0; i < dim; i++) {
        blockSize[i] = 10;
    }
    StratifiedSampleGenerator pSSampler(blockSize, seed);

//    std::cout << "Testing samplers... " << std::endl;
//    std::cout << "Naive... " << std::endl;
    testSampler(pNSampler, dim, numSamples, analyticResult, 1e-2);
//    std::cout << "Halton... " << std::endl;
    testSampler(pHSampler, dim, numSamples, analyticResult, 1e-3);
//    std::cout << "LatinHypercube... " << std::endl;
    testSampler(pLHSampler, dim, numSamples, analyticResult, 1e-3);
//    std::cout << "Stratified... " << std::endl;
    testSampler(pSSampler, dim, numSamples, analyticResult, 1e-3);
//    std::cout << "Finished testing quadrature sample generators." << std::endl;
}
