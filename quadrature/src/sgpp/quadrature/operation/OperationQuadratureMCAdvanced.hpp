// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at 
// sgpp.sparsegrids.org

#ifndef OPERATIONQUADRATUREMCADVANCED_HPP
#define OPERATIONQUADRATUREMCADVANCED_HPP

#include <sgpp/base/operation/OperationQuadrature.hpp>
#include <sgpp/base/grid/Grid.hpp>
#include <sgpp/globaldef.hpp>
#include <sgpp/quadrature/sample/SampleGenerator.hpp>


namespace SGPP {
namespace quadrature {

/**
 * Typedef for general functions that can be passed to integration methods. Requires three parameters. First, the dimensionality, then dim-many coordinates, and then further client data for the function at hand.
 */
typedef double (*FUNC)(int, double*, void*);

/**
 * Quadrature on any sparse grid (that has OperationMultipleEval implemented)
 * using various Monte Carlo Methods (Advanced).
 */

class OperationQuadratureMCAdvanced: public SGPP::base::OperationQuadrature {

public:

	/**
	 * @brief Constructor of OperationQuadratureMCAdvanced, specifying a grid
	 * object and the number of samples to use.
	 *
	 * @param grid Reference to the grid object
	 * @param numberOfSamples Number of Monte Carlo samples
	 */
	OperationQuadratureMCAdvanced(SGPP::base::Grid& grid, int numberOfSamples);

	/**
	 * @brief Constructor of OperationQuadratureMCAdvanced, specifying dimensions
	 * and the number of samples to use.
	 *
	 * @param dimensions dimensionality of this problem
	 * @param numberOfSamples Number of Monte Carlo samples
	 */
	OperationQuadratureMCAdvanced(size_t dimensions, int numberOfSamples);

	virtual ~OperationQuadratureMCAdvanced() {
	}

	/**
	 * @brief Quadrature using advanced MC in @f$\Omega=[0,1]^d@f$.
	 *
	 * @param alpha Coefficient vector for current grid
	 */
	virtual double doQuadrature(SGPP::base::DataVector& alpha);

	/**
	 * @brief Quadrature of an arbitrary function using
	 * advanced MC in @f$\Omega=[0,1]^d@f$.
	 *
	 * @param func The function to integrate
	 * @param clientdata Optional data to pass to FUNC
	 */
	double doQuadratureFunc(FUNC func, void* clientdata);

	/**
	 * @brief Quadrature of the @f$L^2@f$-norm of the error,
	 * @f$ ||f(x)-u(x)||_{L^2} @f$, between a given function and the
	 * current sparse grid function using
	 * advanced MC in @f$\Omega=[0,1]^d@f$.
	 *
	 * @param func The function @f$f(x)@f$
	 * @param clientdata Optional data to pass to FUNC
	 * @param alpha Coefficient vector for current grid
	 */
	double doQuadratureL2Error(FUNC func, void* clientdata,
			SGPP::base::DataVector& alpha);

	/**
	 * @brief Initialize SampleGenerator for NaiveMC
	 */
	void useNaiveMonteCarlo();

	/**
	 * Initialize SampleGenerator for QuasiMC
	 */
	void useQuasiMonteCarlo();

	/**
	 * Initialize SampleGenerator for QuasiMC_Scrambled
	 */
	void useQuasiMonteCarloScrambled();

	/**
	 * @brief Initialize SampleGenerator for StratifiedMC
	 * 
	 * @param n Array of dimension proerties
	 */
	void useStratifiedMonteCarlo(long long int* n);
	void useStratifiedMonteCarlo(long long int* n, int dim) {
		useStratifiedMonteCarlo(n);
	}
	;

	/**
	 * @brief Initialize SampleGenerator for LatinHypercubeMC
	 */
	void useLatinHypercubeMonteCarlo();

	/**
	 * @brief Initialize SampleGenerator for SSobol
	 *
	 * @param scrambling The Type of scrambling to use, according to SSOBOL
	 */
	void useSSobol(int scrambling = 1);

	/**
	 * @brief Method returns the total number of samples which can be generated
	 * according to the sample generator settings (dimensions and subdivision into strata)
	 *
	 * @return size_t dimension of samples
	 */
	size_t getDimensions();

protected:
	// Pointer to the grid object
	SGPP::base::Grid* grid;
	// Number of MC samples
	size_t numberOfSamples;
	// number of dimensions (same as in Grid, if given)
	size_t dimensions;

	//SampleGenerator Instance 
	SGPP::quadrature::SampleGenerator* myGenerator;

};

}
}

#endif /* OPERATIONQUADRATUREMCADVANCED_HPP */