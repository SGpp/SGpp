/*
 * GeneralQuadrature.hpp
 *
 *  Created on: Jan 21, 2015
 *      Author: petz
 */

#ifndef GENERALQUADRATURE_HPP_
#define GENERALQUADRATURE_HPP_

#include <sgpp/combigrid/quadratures/AbstractQuadrature.hpp>

// need to include these headers to reference the corresponding "calculateCoefficients" static methods!
#include <sgpp/combigrid/quadratures/TrapezoidalRule.hpp>
#include <sgpp/combigrid/quadratures/ClenshawCurtisQuadrature.hpp>
#include <sgpp/combigrid/quadratures/GaussPattersonQuadrature.hpp>
#include <sgpp/combigrid/quadratures/BasuQuadrature.hpp>



namespace combigrid {

typedef struct{
	Stretching key; //
	int value; //
}QRtouple;


template<typename _Tp>
class QuadratureRule: public AbstractQuadratureRule<_Tp> {


public:

	/**
	 * A general quadrature rule class, performing integration over a combigrid with arbitrary stretching along different directions.
	 * For example, integration for a 2-D grid with equidistant stretching along x-direction and say Chebyshev-stretched grid along y, will be
	 * calculated with the coefficients for the trapezoidal rule and the clenshaw-curtis rule respectively.
	 * Based on the underlying stretching type, the class AUTOMATICALLY chooses the most optimal set of coefficients according to the following
	 * mapping:
	 *
	 * 		abscissas		<--------->   coefficients
	 *
	 * 1) EQUIDISTANT 	 	<--------->  Trapezoidal rule
	 * 2) CHEBYSHEV      	<--------->  Clenshaw-Curtis rule
	 * 3) LEGENDRE 			<--------->  Gauss-Patterson rule
	 * 4) BASU 				<--------->  Basu quadrature rule
	 * 5) OTHERS			<--------->  Trapezoidal rule
	 *
	 * @param max_lvl set the maximum level to pre-comupte the coefficients for!
	 * @return an instance of the class QuadratureRule with pre-compute coefficients for trapezoidal, clenshaw-curtis, Gauss-Patterson,
	 * and Basu quadrature rules up to (including) level = max_lvl
	 *
	 */
	QuadratureRule(int max_lvl = 9);

	/***
	 * @param max_lvls an integer vector setting the maximum level to pre-compute the coefficients for each of the (currently supported)
	 * 4 basic quadrature rules ( keeping the ordering ) {TrapezoidalRule, Clenshaw-Curtis,Gauss-Patterson, Basu}.
	 * Observing the ordering in brackets this vector should be of size 4 an
	 *
	 * @return an instance of the class QuadratureRule with pre-compute coefficients for trapezoidal, clenshaw-curtis, Gauss-Patterson,
	 * and Basu quadrature rules up to (including) level = max_lvl
	 *
	 */
	QuadratureRule(std::vector<int> max_lvls);

	virtual ~QuadratureRule();

	_Tp integrate(CombiGrid<_Tp>* grids, _Tp (*f)(std::vector<double>) = NULL);

	static void calculateCoefficients(int in_level, _Tp** out_coefs, QRtouple rule);

private:

	/**
	 * Do a full grid trapezoidal rule over a single dim-dimensional FULL GRID.
	 *
	 */
	_Tp quadrature_full_grid(int dim, _Tp (*f)(std::vector<double>),
			FGridContainer<_Tp>* gridContainer, bool interpolate);

	std::vector<int> MAX_LEVELS;
	/* Container for the pre-computed coefs for levels 1,2... MAX_LEVELS. coefficients[d] is a pointer to a _Tp array
	 * of size 2^(d+1) + 1 containing the trapezoidal rule coefficients for a 1-D grid of size 2^(d+1) + 1*/
	std::vector<_Tp**> coefficients;

	/**
	 *
	 * Returns the quadrature rule index corresponding to the str Stretching.
	 *
	 */
	int stretchingToIdx(Stretching str,QRtouple* qrtuples, int size);
	Stretching idxToStretching(int idx,QRtouple* qrtuples,int size);

};

}

#endif /* GENERALQUADRATURE_HPP_ */
