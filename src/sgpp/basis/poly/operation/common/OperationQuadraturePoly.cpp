/******************************************************************************
* Copyright (C) 2009 Technische Universitaet Muenchen                         *
* This file is part of the SG++ project. For conditions of distribution and   *
* use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
******************************************************************************/
// @author Kilian Roehner (roehner@tum.de)

#include "sgpp.hpp"

#include "basis/poly/operation/common/OperationQuadraturePoly.hpp"

#include "basis/basis.hpp"
#include "data/DataVector.hpp"
#include "grid/GridStorage.hpp"
#include <sstream>
#include <cmath>

namespace sg
{
namespace base
{

double OperationQuadraturePoly::doQuadrature(DataVector& alpha)
{
	double res = 0;
	GridStorage::index_type index;
	GridStorage::grid_map_iterator end_iter = storage->end();
	
	for(GridStorage::grid_map_iterator iter = storage->begin(); iter != end_iter; iter++) {
		//    index = *(iter->first);
		//    std::cout << iter->second << " " << iter->first->getLevelSum() << " " << pow(2.0, -static_cast<double>(iter->first->getLevelSum())) << std::endl;
		if(base.getDegree() == 2) {
			res += pow(2.0, (2.0-log(3.0)/log(2.0))*static_cast<double>(iter->first->dim())-static_cast<double>(iter->first->getLevelSum()))*alpha.get(iter->second);
		}
		/*else if(base->degree == 3) {
			res += pow(2.0, -static_cast<double>(iter->first->getLevelSum()))*alpha.get(iter->second);
		}*/
		else {
			res += 0; // if this case occurs, something has gone very wrong...
		}
	}
	return res;
}

}
}
