/* ****************************************************************************
* Copyright (C) 2009 Technische Universitaet Muenchen                         *
* This file is part of the SG++ project. For conditions of distribution and   *
* use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
**************************************************************************** */
// @author Julian Valentin (julian.valentin@stud.mathematik.uni-stuttgart.de)

#ifndef LINEARCLENSHAWCURTISBASE_HPP
#define LINEARCLENSHAWCURTISBASE_HPP

#include "base/basis/basis.hpp"
#include "base/tools/CosineTable.hpp"

#include <cmath>

namespace sg
{
namespace base
{

template<class LT, class IT>
class LinearClenshawCurtisBasis : public Basis<LT, IT>
{
protected:
    const CosineTable *cosine_table;
    
public:
    LinearClenshawCurtisBasis() : LinearClenshawCurtisBasis(NULL) {}
    LinearClenshawCurtisBasis(const CosineTable *cosine_table) : cosine_table(cosine_table) {}
    
    inline double clenshawCurtisPoint(double h, IT index) const
    {
        // h = 1.0 / static_cast<double>(1 << level)
        if (cosine_table == NULL)
        {
            return (cos(M_PI * (1.0 - static_cast<double>(index) * h)) + 1.0) / 2.0;
        } else
        {
            return (cosine_table->lookUp(
                    M_PI * (1.0 - static_cast<double>(index) * h)) + 1.0) / 2.0;
        }
    }
    
    inline double eval(LT level, IT index, double x)
    {
        if (level == 0)
        {
            if (index == 0)
            {
                return 1.0 - x;
            } else
            {
                return x;
            }
        } else
        {
            double h = 1.0 / static_cast<double>(1 << level);
            double x0 = clenshawCurtisPoint(h, index - 1);
            double x2 = clenshawCurtisPoint(h, index + 1);
            
            if ((x <= x0) || (x >= x2))
            {
                return 0.0;
            }
            
            double x1 = clenshawCurtisPoint(h, index);
            
            if (x < x1)
            {
                return 1.0 - (x1 - x) / (x1 - x0);
            } else
            {
                return (x2 - x) / (x2 - x1);
            }
        }
    }
};

typedef LinearClenshawCurtisBasis<unsigned int, unsigned int> SLinearClenshawCurtisBase;

}
}

#endif /* LINEARCLENSHAWCURTISBASE_HPP */
