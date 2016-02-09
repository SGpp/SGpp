/*
 * MSE.cpp
 *
 *  Created on: Feb 8, 2016
 *      Author: perun
 */

#include "MSE.hpp"

#include <sgpp/globaldef.hpp>


namespace SGPP {
namespace datadriven {

MSE::MSE() {
	// TODO Auto-generated constructor stub

}

MSE::~MSE() {
	// TODO Auto-generated destructor stub
}

double MSE::operator()(DataVector& predictedValues, DataVector& trueValues){
	DataVector tmp(predictedValues);
	tmp.sub(trueValues);
	double error =  tmp.l2Norm();
	return error*error/static_cast<double>(tmp.getSize());
}

} /* namespace datadriven */
} /* namespace SGPP */
