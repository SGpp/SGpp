/******************************************************************************
* Copyright (C) 2009 Technische Universitaet Muenchen                         *
* This file is part of the SG++ project. For conditions of distribution and   *
* use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
******************************************************************************/
// @author Julian Valentin (julian.valentin@stud.mathematik.uni-stuttgart.de)

#include "base/basis/modbspline/operation/OperationEvalGradientModBspline.hpp"

#include "base/exception/operation_exception.hpp"
#include "base/algorithm/GetAffectedBasisFunctionsGradient.hpp"

#include <tuple>

namespace sg
{
namespace base
{

double OperationEvalGradientModBspline::evalGradient(
        DataVector& alpha, const std::vector<double> &point, DataVector &gradient)
{
    typedef std::vector<std::tuple<size_t, double, DataVector> > IndexValVector;
    
    IndexValVector vec;
    GetAffectedBasisFunctionsGradient<SModBsplineBase> ga(storage);
    
    ga(base, point, vec);
    
    double result = 0.0;
    
    gradient = DataVector(point.size());
    gradient.setAll(0.0);
    
    for (IndexValVector::iterator iter = vec.begin(); iter != vec.end(); iter++)
    {
        double coeff = alpha[std::get<0>(*iter)];
        DataVector &current_gradient = std::get<2>(*iter);
        
        current_gradient.mult(coeff);
        gradient.add(current_gradient);
        
        result += coeff * std::get<1>(*iter);
    }
    
    return result;
}

}
}
