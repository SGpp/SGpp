/******************************************************************************
* Copyright (C) 2009 Technische Universitaet Muenchen                         *
* This file is part of the SG++ project. For conditions of distribution and   *
* use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
******************************************************************************/
// @author Julian Valentin (julian.valentin@stud.mathematik.uni-stuttgart.de)

#include "base/basis/wavelet/boundary/operation/OperationEvalHessianWaveletBoundary.hpp"

#include "base/exception/operation_exception.hpp"
#include "base/algorithm/GetAffectedBasisFunctionsHessian.hpp"

#include <tuple>

namespace sg
{
namespace base
{

double OperationEvalHessianWaveletBoundary::evalHessian(
        DataVector &alpha, const std::vector<double> &point,
        DataVector &gradient, DataMatrix &hessian)
{
    typedef std::vector<std::tuple<size_t, double, DataVector, DataMatrix> > IndexValVector;
    
    IndexValVector vec;
    GetAffectedBasisFunctionsHessian<SWaveletBoundaryBase > ga(storage);
    
    ga(base, point, vec);
    
    double result = 0.0;
    
    gradient = DataVector(point.size());
    gradient.setAll(0.0);
    
    hessian = DataMatrix(point.size(), point.size());
    hessian.setAll(0.0);
    
    for (IndexValVector::iterator iter = vec.begin(); iter != vec.end(); iter++)
    {
        double coeff = alpha[std::get<0>(*iter)];
        DataVector &current_gradient = std::get<2>(*iter);
        DataMatrix &current_hessian = std::get<3>(*iter);
        
        current_gradient.mult(coeff);
        current_hessian.mult(coeff);
        
        gradient.add(current_gradient);
        hessian.add(current_hessian);
        
        result += coeff * std::get<1>(*iter);
    }
    
    return result;
}

}
}
