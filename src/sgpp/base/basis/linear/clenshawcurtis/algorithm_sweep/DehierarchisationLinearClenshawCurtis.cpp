/******************************************************************************
* Copyright (C) 2009 Technische Universitaet Muenchen                         *
* This file is part of the SG++ project. For conditions of distribution and   *
* use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
******************************************************************************/
// @author Alexander Heinecke (Alexander.Heinecke@mytum.de)

#include "base/basis/linear/clenshawcurtis/algorithm_sweep/DehierarchisationLinearClenshawCurtis.hpp"

namespace sg
{
namespace base
{

DehierarchisationLinearClenshawCurtis::DehierarchisationLinearClenshawCurtis(
        GridStorage *storage, const CosineTable *cosine_table) :
    storage(storage),
    linear_cc_base(cosine_table)
{
}

DehierarchisationLinearClenshawCurtis::~DehierarchisationLinearClenshawCurtis()
{
}

void DehierarchisationLinearClenshawCurtis::operator()(DataVector &source, DataVector &result,
                                                       grid_iterator &index, size_t dim)
{
    double left_boundary;
    double right_boundary;
    size_t seq;
    
    index.left_levelzero(dim);
    seq = index.seq();
    left_boundary = source[seq];
    
    index.right_levelzero(dim);
    seq = index.seq();
    right_boundary = source[seq];
    
    if (!index.hint())
    {
        index.top(dim);
        
        if (!storage->end(index.seq()))
        {
            rec(source, result, index, dim, left_boundary, right_boundary);
        }
        
        index.left_levelzero(dim);
    }
}

void DehierarchisationLinearClenshawCurtis::rec(DataVector& source, DataVector& result,
        grid_iterator& index, size_t dim, double fl, double fr)
{
    size_t seq = index.seq();
    double fm = source[seq];
    
    grid_iterator::index_t i;
    grid_iterator::level_t l;
    index.get(dim, l, i);
    
    double h = 1.0 / static_cast<double>(1 << l);
    double x0 = linear_cc_base.clenshawCurtisPoint(h, i - 1);
    double x1 = linear_cc_base.clenshawCurtisPoint(h, i);
    double x2 = linear_cc_base.clenshawCurtisPoint(h, i + 1);
    double t = (x1 - x0) / (x2 - x0);
    
    fm += ((1-t) * fl + t * fr);
    result[seq] = fm;
    
    if (index.hint() == false)
    {
        index.left_child(dim);
        
        if (!storage->end(index.seq()))
        {
            rec(source, result, index, dim, fl, fm);
        }
        
        index.step_right(dim);
        
        if (!storage->end(index.seq()))
        {
            rec(source, result, index, dim, fm, fr);
        }
        
        index.up(dim);
    }
}

}
}
