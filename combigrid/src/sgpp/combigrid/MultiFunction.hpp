/*
 * MultiFunction.hpp
 *
 *  Created on: 13.01.2016
 *      Author: david
 */

#ifndef COMBIGRID_SRC_SGPP_COMBIGRID_MULTIFUNCTION_HPP_
#define COMBIGRID_SRC_SGPP_COMBIGRID_MULTIFUNCTION_HPP_

#include <sgpp/globaldef.hpp>
#include <sgpp/base/datatypes/DataVector.hpp>

#include <functional>

namespace SGPP {
namespace combigrid {

class MultiFunction {
public:
	typedef std::function<SGPP::float_t(base::DataVector const &)> function_type;

private:
	function_type func;

public:
	/**
	 * for function pointers
	 */
	MultiFunction(SGPP::float_t (*ptr)(base::DataVector const &));

	/**
	 * for lambdas or function objects
	 */
	template<typename T> explicit MultiFunction(T f) : func(f) {

	}

	SGPP::float_t operator()(base::DataVector const &vec);
	SGPP::float_t call(base::DataVector const &vec);
};

} /* namespace combigrid */
} /* namespace SGPP */

#endif /* COMBIGRID_SRC_SGPP_COMBIGRID_MULTIFUNCTION_HPP_ */
