/******************************************************************************
* Copyright (C) 2009 Technische Universitaet Muenchen                         *
* This file is part of the SG++ project. For conditions of distribution and   *
* use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
******************************************************************************/
// @author Julian Valentin (julian.valentin@stud.mathematik.uni-stuttgart.de)

#include "base/basis/wavelet/noboundary/operation/OperationEvalWavelet.hpp"

#include "base/algorithm/GetAffectedBasisFunctions.hpp"

namespace sg
{
namespace base
{

double OperationEvalWavelet::eval(DataVector &alpha, std::vector<double> &point)
{
    typedef std::vector<std::pair<size_t, double> > IndexValVector;
    
    IndexValVector vec;
    GetAffectedBasisFunctions<SWaveletBase > ga(storage);
    
    ga(base, point, vec);
    
    double result = 0.0;
    
    for (IndexValVector::iterator iter = vec.begin(); iter != vec.end(); iter++)
    {
        result += iter->second * alpha[iter->first];
    }
    
    return result;
}

}
}
