/*
 * Optimization.h
 *
 *  Created on: Jan 25, 2016
 *      Author: liedtkjn
 */

#ifndef OPTIMIZATION_H_
#define OPTIMIZATION_H_

#include <vector>
#include <sgpp/globaldef.hpp>

namespace SGPP {
namespace combigrid {

namespace optimize {

	// TODO destructor

	class func_base{
	public:
		virtual float_t operator() (float_t) = 0;
	};

	class monicPoly : public func_base {
	public:
		std::vector<float_t> coeff;
		virtual float_t operator() (float_t x);
		// constructors:
		monicPoly(const size_t degree)
			: coeff(degree) {}
		monicPoly(const std::vector<float_t>& v)
			: coeff(v) {}
		monicPoly(const float_t* c, size_t degree)
			: coeff(std::vector<float_t>(c, c + degree)) {}
	};

	class Poly : public func_base {
	public:
		std::vector<float_t> coeff;    // a vector of size nterms i.e. 1+degree
		virtual float_t operator() (float_t x);
		// constructors:
		Poly(const size_t degree)
			: coeff(1 + degree) {}
		Poly(const std::vector<float_t>& v)
			: coeff(v) {}
		Poly(const float_t* c, size_t degree)
			: coeff(std::vector<float_t>(c, 1 + c + degree)) {}
	};

	float_t glomin(float_t a, float_t b, float_t c, float_t m, float_t e, float_t t,
		func_base& f, float_t &x);
	float_t local_min(float_t a, float_t b, float_t t, func_base& f,
		float_t &x);
	float_t local_min_rc(float_t &a, float_t &b, int &status, float_t value);
	float_t r8_abs(float_t x);
	float_t r8_epsilon();
	float_t r8_max(float_t x, float_t y);
	float_t r8_sign(float_t x);
	void timestamp();

	// === simple wrapper functions
	// === for convenience and/or compatibility
	float_t glomin(float_t a, float_t b, float_t c, float_t m, float_t e, float_t t,
		float_t f(float_t x), float_t &x);
	float_t local_min(float_t a, float_t b, float_t t, float_t f(float_t x),
		float_t &x);

}

}
} /* namespace SGPP */
#endif /* OPTIMIZATION_H_ */
