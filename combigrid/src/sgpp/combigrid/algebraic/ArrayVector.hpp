/*
 * ArrayVector.hpp
 *
 *  Created on: 11.12.2015
 *      Author: david
 */

#ifndef COMBIGRID_SRC_SGPP_COMBIGRID_ALGEBRAIC_ARRAYVECTOR_HPP_
#define COMBIGRID_SRC_SGPP_COMBIGRID_ALGEBRAIC_ARRAYVECTOR_HPP_

#include <vector>
#include <algorithm>
#include <sgpp/globaldef.hpp>
#include "ScalarVector.hpp"
#include <cmath>

namespace sgpp{
namespace combigrid {

template<typename Scalar, typename V> class ArrayVector {
	std::vector<V> values;

	void ensureMinimumSize(size_t minSize) {
		if (values.size() == 0) {
			values.push_back(V::zero());
		}
		while (values.size() < minSize) {
			values.push_back(values.back());
		}
	}

public:
	ArrayVector(std::vector<V> const &values) :
			values(values) {

	}

	ArrayVector(V value) :
			values(1, value) {

	}

	ArrayVector() : values() {

	}

	ArrayVector(ArrayVector<Scalar, V> const &other) :
			values(other.values) {

	}

	ArrayVector<Scalar, V> &operator=(ArrayVector<Scalar, V> const &other) {
		values = other.values;
		return *this;
	}

	std::vector<V> const &getValues() const {
		return values;
	}

	V get(size_t i) const {
		return (*this)[i];
	}

	V &at(size_t i) {
		ensureMinimumSize(i + 1);
		return values[i];
	}

	V &operator[](size_t i) {
		ensureMinimumSize(i + 1);
		return values[i];
	}

	/**
	 * This operator is unsafe because it does no range checking
	 */
	V const &operator[](size_t i) const {
		return values[i];
	}

	size_t size() const {
		return values.size();
	}

	void add(ArrayVector<Scalar, V> const &other) {
		ensureMinimumSize(other.size());

		size_t otherSize = other.size();
		for (size_t i = 0; i < otherSize; ++i) {
			values[i].add(other.values[i]);
		}

		for (size_t i = otherSize; i < size(); ++i) {
			values[i].add(other.values.back());
		}
	}

	void sub(ArrayVector<Scalar, V> const &other) {
		ensureMinimumSize(other.size());

		size_t otherSize = other.size();
		for (size_t i = 0; i < otherSize; ++i) {
			values[i].sub(other.values[i]);
		}

		for (size_t i = otherSize; i < size(); ++i) {
			values[i].sub(other.values.back());
		}
	}

	void componentwiseMult(ArrayVector<Scalar, V> const &other) {
		ensureMinimumSize(other.size());

		size_t otherSize = other.size();
		for (size_t i = 0; i < otherSize; ++i) {
			values[i].componentwiseMult(other.values[i]);
		}

		for (size_t i = otherSize; i < size(); ++i) {
			values[i].componentwiseMult(other.values.back());
		}
	}

	void scalarMult(Scalar const &factor) {
		for (size_t i = 0; i < values.size(); ++i) {
			values[i].scalarMult(factor);
		}
	}

	Scalar norm() const {
		Scalar result = Scalar();

		for (auto &val : values) {
			Scalar n = val.norm();
			result += n * n;
		}

		return sqrt(result);
	}

	static ArrayVector<Scalar, V> zero() {
		return ArrayVector<Scalar, V>(V::zero());
	}

	static ArrayVector<Scalar, V> one() {
		return ArrayVector<Scalar, V>(V::one());
	}
};

typedef ArrayVector<double, FloatScalarVector> FloatArrayVector;

} /* namespace combigrid */
} /* namespace sgpp*/

#endif /* COMBIGRID_SRC_SGPP_COMBIGRID_ALGEBRAIC_ARRAYVECTOR_HPP_ */
