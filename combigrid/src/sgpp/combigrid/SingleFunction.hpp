/*
 * SingleFunction.hpp
 *
 *  Created on: 27.02.2016
 *      Author: david
 */

#ifndef COMBIGRID_SRC_SGPP_COMBIGRID_SINGLEFUNCTION_HPP_
#define COMBIGRID_SRC_SGPP_COMBIGRID_SINGLEFUNCTION_HPP_

#include <sgpp/globaldef.hpp>
#include <functional>

namespace SGPP {
namespace combigrid {

class SingleFunction {
public:
	typedef std::function<SGPP::float_t(SGPP::float_t)> function_type;

private:
	function_type func;

public:
	/**
	 * for function pointers
	 */
	SingleFunction(SGPP::float_t (*ptr)(SGPP::float_t));

	/**
	 * for lambdas or function objects
	 */
	template<typename T> explicit SingleFunction(T f) : func(f) {

	}

	float_t operator()(float_t param);
	float_t call(float_t param);
};

} /* namespace combigrid */
} /* namespace SGPP */

#endif /* COMBIGRID_SRC_SGPP_COMBIGRID_SINGLEFUNCTION_HPP_ */
