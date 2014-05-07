/**
 * Note by Julian Valentin:
 * 
 * This code is (C) 2007 by Numerical Recipes Software.
 * It is not allowed to be distributed publicly in source code form.
 * Furthermore, it's not allowed to make shared libraries (e.g. DLLs) containing this code
 * or to make commercial binaries containing this code.
 * See http://www.nr.com/licenses/.
 */

#ifndef SGPP_OPT_OPTIMIZER_NUMRECIPES_BRENT_HPP
#define SGPP_OPT_OPTIMIZER_NUMRECIPES_BRENT_HPP

#include "opt/optimizer/numrecipes/Bracketmethod.hpp"

#include <cmath>
#include <cstddef>
#include <limits>
#include <iostream>

namespace sg
{
namespace opt
{
namespace optimizer
{
namespace numrecipes
{

struct Brent : Bracketmethod
{
    double xmin, fmin;
    const double tol;
    
    Brent(const double toll = 3.0e-8) : tol(toll)
    {
    }
    
    template <class T>
    double minimize(T &func)
    {
        const size_t ITMAX = 100;
        const double CGOLD = 0.3819660;
        const double ZEPS = std::numeric_limits<double>::epsilon() * 1.0e-3;
        double a, b, d = 0.0, etemp, fu, fv, fw, fx;
        double p, q, r, tol1, tol2, u, v, w, x, xm;
        double e = 0.0;
        
        a = ((ax < cx) ? ax : cx);
        b = ((ax > cx) ? ax : cx);
        
        x = bx;
        w = x;
        v = x;
        
        fx = func(x);
        fw = fx;
        fv = fx;
        
        for (size_t iter = 0; iter < ITMAX; iter++)
        {
            xm = 0.5 * (a + b);
            tol1 = tol * std::abs(x) + ZEPS;
            tol2 = 2.0 * tol1;
            
            if (std::abs(x - xm) <= tol2 - 0.5 * (b - a))
            {
                xmin = x;
                fmin = fx;
                
                return xmin;
            }
            
            if (std::abs(e) > tol1)
            {
                r = (x - w) * (fx - fv);
                q = (x - v) * (fx - fw);
                p = (x - v) * q - (x - w) * r;
                q = 2.0 * (q - r);
                
                if (q > 0.0)
                {
                    p = -p;
                }
                
                q = std::abs(q);
                etemp = e;
                e = d;
                
                if ((std::abs(p) >= std::abs(0.5*q*etemp)) || (p <= q * (a-x)) || (p >= q * (b-x)))
                {
                    e = ((x >= xm) ? (a - x) : (b - x));
                    d = CGOLD * e;
                } else
                {
                    d = p / q;
                    u = x + d;
                    if ((u - a < tol2) || (b - u < tol2))
                    {
                        d = tol1 * std::copysign(1.0, xm - x);
                    }
                }
            } else
            {
                e = ((x >= xm) ? (a - x) : (b - x));
                d = CGOLD * e;
            }
            
            u = ((std::abs(d) >= tol1) ? (x + d) : (x + tol1 * std::copysign(1.0, d)));
            fu = func(u);
            
            if (fu <= fx)
            {
                if (u >= x)
                {
                    a = x;
                } else
                {
                    b = x;
                }
                
                shft3(v, w, x, u);
                shft3(fv, fw, fx, fu);
            } else
            {
                if (u < x)
                {
                    a = u;
                } else
                {
                    b = u;
                }
                
                if ((fu <= fw) || (w == x))
                {
                    v = w;
                    w = u;
                    fv = fw;
                    fw = fu;
                } else if ((fu <= fv) || (v == x) || (v == w))
                {
                    v = u;
                    fv = fu;
                }
            }
        }
        
        std::cerr << "Too many iterations in brent\n";
        
        return 0.0;
    }
};

}
}
}
}

#endif
