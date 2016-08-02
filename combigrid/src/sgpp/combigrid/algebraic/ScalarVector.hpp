/*
 * ScalarVector.hpp
 *
 *  Created on: 11.12.2015
 *      Author: david
 */

#ifndef COMBIGRID_SRC_SGPP_COMBIGRID_ALGEBRAIC_SCALARVECTOR_HPP_
#define COMBIGRID_SRC_SGPP_COMBIGRID_ALGEBRAIC_SCALARVECTOR_HPP_

#include <sgpp/globaldef.hpp>
#include <cmath>

namespace SGPP {
namespace combigrid {

template<typename Scalar> class ScalarVector {
	Scalar val;
public:
	ScalarVector()
		: val(Scalar()) {

	}

	ScalarVector(Scalar const &value)
		: val(value) {

	}

	ScalarVector(ScalarVector const &other)
		: val(other.val) {

	}

	ScalarVector<Scalar> &operator=(ScalarVector<Scalar> const &other) {
		val = other.val;
		return *this;
	}

	Scalar const &value() const {
		return val;
	}

	Scalar &value() {
		return val;
	}

	Scalar getValue() const {
		return val;
	}

	void add(ScalarVector<Scalar> const &other) {
		val += other.val;
	}

	void sub(ScalarVector<Scalar> const &other) {
		val -= other.val;
	}

	void componentwiseMult(ScalarVector<Scalar> const &other) {
		val *= other.val;
	}

	void scalarMult(Scalar const &factor) {
		val *= factor;
	}

	Scalar norm() const {
		return std::abs(val);
	}

	static ScalarVector<Scalar> zero() {
		return ScalarVector(Scalar(0));
	}

	static ScalarVector<Scalar> one() {
		return ScalarVector(Scalar(1));
	}
};

typedef ScalarVector<SGPP::float_t> FloatScalarVector;

} /* namespace combigrid */
} /* namespace SGPP */

#endif /* COMBIGRID_SRC_SGPP_COMBIGRID_ALGEBRAIC_SCALARVECTOR_HPP_ */
