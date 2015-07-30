#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE SGppQuadratureModule
#include <boost/test/unit_test.hpp>

#include <sgpp_quadrature.hpp>
#include <sgpp_base.hpp>
#include <sgpp/globaldef.hpp>

using namespace SGPP::base;
using namespace SGPP::quadrature;

sg::float_t f(DataVector x) {
	sg::float_t res = 1;
	for (size_t i = 0; i < x.getSize(); i++) {
		res *= 4 * (1 - x[i]) * x[i];
	}
	return res;
}

void testSampler (SampleGenerator& sampler, size_t dim, size_t numSamples)
{
	DataVector sample = DataVector(dim);
	sg::float_t sum = 0.0f;
	
	for (size_t i = 0; i < numSamples; i++)
	{
		sampler.getSample(sample);
		sum += f(sample);
	}
	sum /= static_cast<sg::float_t>(numSamples);
	
	BOOST_CHECK_CLOSE(sum, 2. / 3., 1e-3);
}

BOOST_AUTO_TEST_CASE(testSamplers) {
	size_t dim = 2;
	size_t numSamples = 10000;
	NaiveSampleGenerator pNSampler(dim);
	HaltonSampleGenerator pHSampler(dim);
	LatinHypercubeSampleGenerator pLHSampler(dim, numSamples);
	long long int blockSize[2];
	for (size_t i = 0; i < dim; i++)
	{
		blockSize[i] = 10;
	}
	StratifiedSampleGenerator pSSampler(dim, blockSize);
	
	printf("Testing samplers... \n");
	printf("Naive... \n");
	testSampler(pNSampler, dim, numSamples);
	printf("Halton... \n");
	testSampler(pHSampler, dim, numSamples);
	printf("LatinHypercube... \n");
	testSampler(pLHSampler, dim, numSamples);
	printf("Stratified... \n");
	testSampler(pSSampler, dim, numSamples);
	printf("Finished testing quadrature sample generators. \n");
}