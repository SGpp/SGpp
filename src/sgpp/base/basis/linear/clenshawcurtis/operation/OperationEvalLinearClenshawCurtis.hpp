/* ****************************************************************************
* Copyright (C) 2009 Technische Universitaet Muenchen                         *
* This file is part of the SG++ project. For conditions of distribution and   *
* use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
**************************************************************************** */
// @author Julian Valentin (julian.valentin@stud.mathematik.uni-stuttgart.de)

#ifndef OPERATIONEVALLINEARCLENSHAWCURTIS_HPP
#define OPERATIONEVALLINEARCLENSHAWCURTIS_HPP

#include "base/operation/OperationEval.hpp"
#include "base/grid/GridStorage.hpp"

namespace sg
{
namespace base
{

class OperationEvalLinearClenshawCurtis : public OperationEval
{
public:
    OperationEvalLinearClenshawCurtis(GridStorage *storage) : storage(storage) {}
    virtual ~OperationEvalLinearClenshawCurtis() {}
    
    virtual double eval(DataVector &alpha, std::vector<double> &point);
    
protected:
    GridStorage *storage;
};

}
}

#endif /* OPERATIONEVALLINEARCLENSHAWCURTIS_HPP */
