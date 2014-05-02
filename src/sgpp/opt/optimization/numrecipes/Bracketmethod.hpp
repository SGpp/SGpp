#ifndef SGPP_OPT_OPTIMIZATION_NUMRECIPES_BRACKETMETHOD_HPP
#define SGPP_OPT_OPTIMIZATION_NUMRECIPES_BRACKETMETHOD_HPP

#include <cmath>

namespace sg
{
namespace opt
{
namespace optimization
{
namespace numrecipes
{

struct Bracketmethod
{
    double ax, bx, cx, fa, fb, fc;
    
    template <class T>
    void bracket(const double a, const double b, T &func)
    {
        const double GOLD = 1.618034, GLIMIT = 100.0, TINY = 1.0e-20;
        
        ax = a;
        bx = b;
        
        double fu;
        
        fa = func(ax);
        fb = func(bx);
        
        if (fb > fa)
        {
            std::swap(ax, bx);
            std::swap(fb, fa);
        }
        
        cx = bx + GOLD * (bx - ax);
        fc = func(cx);
        
        while (fb > fc)
        {
            double r = (bx - ax) * (fb - fc);
            double q = (bx - cx) * (fb - fa);
            double u = bx - ((bx-cx) * q - (bx-ax) * r) /
                       (2.0 * std::max(std::abs(q-r), TINY) * std::copysign(1.0, q-r));
            double ulim = bx + GLIMIT * (cx - bx);
            
            if ((bx - u) * (u - cx) > 0.0)
            {
                fu = func(u);
                
                if (fu < fc)
                {
                    ax = bx;
                    bx = u;
                    fa = fb;
                    fb = fu;
                    return;
                } else if (fu > fb)
                {
                    cx = u;
                    fc = fu;
                    return;
                }
                
                u = cx + GOLD * (cx - bx);
                fu = func(u);
            } else if ((cx - u) * (u - ulim) > 0.0)
            {
                fu = func(u);
                
                if (fu < fc)
                {
                    shft3(bx, cx, u, u + GOLD * (u - cx));
                    shft3(fb, fc, fu, func(u));
                }
            } else if ((u - ulim) * (ulim - cx) >= 0.0)
            {
                u = ulim;
                fu = func(u);
            } else
            {
                u = cx + GOLD * (cx - bx);
                fu = func(u);
            }
            
            shft3(ax, bx, cx, u);
            shft3(fa, fb, fc, fu);
        }
    }
    
    inline void shft2(double &a, double &b, const double c)
    {
        a = b;
        b = c;
    }
    
    inline void shft3(double &a, double &b, double &c, const double d)
    {
        a = b;
        b = c;
        c = d;
    }
    
    inline void mov3(double &a, double &b, double &c,
                     const double d, const double e, const double f)
    {
        a = d;
        b = e;
        c = f;
    }
};

}
}
}
}

#endif
