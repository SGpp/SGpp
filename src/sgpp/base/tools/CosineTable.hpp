/* ****************************************************************************
* Copyright (C) 2010 Technische Universitaet Muenchen                         *
* This file is part of the SG++ project. For conditions of distribution and   *
* use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
**************************************************************************** */
// @author Julian Valentin (julian.valentin@stud.mathematik.uni-stuttgart.de)

#ifndef COSINETABLE_HPP
#define COSINETABLE_HPP

#include <vector>

namespace sg
{
namespace base
{

class CosineTable
{
public:
    static const size_t DEFAULT_SIZE = 1 << 20;
    
    CosineTable() : CosineTable(DEFAULT_SIZE) {}
    CosineTable(size_t size) : table(std::vector<double>(size+1, 0.0))
    {
        h = M_PI_2 / static_cast<double>(size);
        hinv = 1.0 / h;
        
        double x = 0.0;
        
        for (size_t i = 0; i <= size; i++)
        {
            table[i] = cos(x);
            // rounding errors can be neglected
            x += h;
        }
    }
    
    inline double lookUp(double x) const
    {
        /*if (x < 1e-5)
        {
            return cos(x);
        }*/
        
        double factor = 1.0;
        
        if (x > M_PI_2)
        {
            x = M_PI - x;
            factor = -1.0;
        } else if (x == M_PI_2)
        {
            return 0.0;
        }
        
        double xl = x * hinv;
        size_t i = static_cast<size_t>(xl);
        double t = xl - static_cast<double>(i);
        return factor * ((1.0-t) * table[i] + t * table[i+1]);
        
        /*size_t index = static_cast<size_t>(x * hinv);
        return factor * table[index];*/
        
        /*if ((x < 1e-5) || (x > M_PI - 1e-5))
        {
            return cos(x);
        } else if (x <= M_PI_2)
        {
            size_t index = static_cast<size_t>(x * hinv + 0.5);
            return table[index];
        } else
        {
            size_t index = static_cast<size_t>((M_PI - x) * hinv + 0.5);
            return -table[index];
        }*/
    }
    
protected:
    std::vector<double> table;
    size_t size;
    double h, hinv;
};

}
}

#endif
