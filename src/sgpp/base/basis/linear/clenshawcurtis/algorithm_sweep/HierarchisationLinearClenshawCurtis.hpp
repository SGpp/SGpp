/* ****************************************************************************
* Copyright (C) 2009 Technische Universitaet Muenchen                         *
* This file is part of the SG++ project. For conditions of distribution and   *
* use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
**************************************************************************** */
// @author Julian Valentin (julian.valentin@stud.mathematik.uni-stuttgart.de)

#ifndef HIERARCHISATIONLINEARCLENSHAWCURTIS_HPP
#define HIERARCHISATIONLINEARCLENSHAWCURTIS_HPP

#include "base/grid/GridStorage.hpp"
#include "base/datatypes/DataVector.hpp"
#include "base/basis/linear/clenshawcurtis/LinearClenshawCurtisBasis.hpp"

namespace sg
{
namespace base
{

class HierarchisationLinearClenshawCurtis
{
protected:
    typedef GridStorage::grid_iterator grid_iterator;
    
public:
    HierarchisationLinearClenshawCurtis(GridStorage *storage);
    HierarchisationLinearClenshawCurtis(GridStorage *storage, const CosineTable *cosine_table);
    ~HierarchisationLinearClenshawCurtis();
    
    void operator()(DataVector &source, DataVector &result, grid_iterator &index, size_t dim);
    
protected:
    GridStorage* storage;
    SLinearClenshawCurtisBase linear_cc_base;
    
    void rec(DataVector &source, DataVector &result, grid_iterator &index,
             size_t dim, double fl, double fr);
};
    
}
}

#endif /* HIERARCHISATIONLINEARCLENSHAWCURTIS_HPP */
