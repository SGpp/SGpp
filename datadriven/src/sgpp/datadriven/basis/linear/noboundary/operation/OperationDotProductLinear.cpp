/******************************************************************************
* Copyright (C) 2009 Technische Universitaet Muenchen                         *
* This file is part of the SG++ project. For conditions of distribution and   *
* use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
******************************************************************************/
// @author Alexander Heinecke (Alexander.Heinecke@mytum.de), Dirk Pflueger (pflueged@in.tum.de)


#include <sgpp/base/basis/linear/noboundary/LinearBasis.hpp>

#include <sgpp/datadriven/basis/linear/noboundary/operation/OperationDotProductLinear.hpp>

using namespace sg::base;
namespace sg {
  namespace datadriven {

    double OperationDotProductLinear::eval(base::DataVector& x1, base::DataVector& x2) {
      base::LinearBasis<unsigned int, unsigned int> base;
      GridStorage::index_type::level_type work_level = 1;
      GridStorage::index_type::index_type work_index;
	  GridStorage::index_type::level_type temp;
      double result = 0;
      //GridStorage::grid_iterator working;
      //for (GridStorage::grid_iterator working = storage->begin(); working != storage->end(); working++){
      for (size_t i = 0; i < storage->size(); i++){
    	  GridStorage::index_type working = storage->get(i);
    	  double value1 = 1.0;
    	  double value2 = 1.0;
    	  for (size_t d = 0; d < storage->dim(); d++){

                    working.get(d, temp, work_index);

                    value1 *= base.eval(work_level, work_index,
                                                  x1[d]);
                    value2 *= base.eval(work_level, work_index,
                                                                      x2[d]);

    	  //}
    	  }
    	  result += value1*value2;
      }
      return result;
    }

  }
}

