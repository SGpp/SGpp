#ifndef SGPP_OPT_OPTIMIZATION_NUMRECIPES_FRPRMN_HPP
#define SGPP_OPT_OPTIMIZATION_NUMRECIPES_FRPRMN_HPP

#include "opt/optimization/numrecipes/Linemethod.hpp"
#include "opt/tools/Printer.hpp"

#include <vector>
#include <cstddef>
#include <cmath>

namespace sg
{
namespace opt
{
namespace optimization
{
namespace numrecipes
{

template <class T>
struct Frprmn : Linemethod<T>
{
    size_t iter;
    double fret;
    using Linemethod<T>::func;
    using Linemethod<T>::linmin;
    using Linemethod<T>::p;
    using Linemethod<T>::xi;
    const double ftol;
    const size_t itmax;
    
    Frprmn(T &funcd, const double ftoll = 3.0e-8, const size_t itmaxx = 200) :
            Linemethod<T>(funcd), ftol(ftoll), itmax(itmaxx)
    {
    }
    
    std::vector<double> minimize(const std::vector<double> &pp)
    {
        /*if (verbosity_level >= 1)
        {
            Output::printStatusBegin("Optimizing (FRPR method)...");
        }*/
        
        const double EPS = 1.0e-18;
        const double GTOL = 1.0e-8;
        
        double gg, dgg;
        size_t n = pp.size();
        
        p = pp;
        
        std::vector<double> g(n, 0.0), h(n, 0.0);
        
        xi.resize(n);
        
        double fp = func(p);
        
        func.df(p, xi);
        
        for (int j = 0; j < n; j++)
        {
            g[j] = -xi[j];
            xi[j] = g[j];
            h[j] = g[j];
        }
        
        size_t its;
        
        for (its = 0; its < itmax; its++)
        {
            if (its % 10 == 0)
            {
                std::stringstream msg;
                msg << its << " steps, f(x) = " << fp;
                tools::printer.printStatusUpdate(msg.str());
            }
            
            iter = its;
            fret = linmin();
            
            if (2.0 * std::abs(fret - fp) <= ftol * (std::abs(fret) + std::abs(fp) + EPS))
            {
                break;
            }
            
            fp = fret;
            func.df(p, xi);
            
            /*std::cout << "\nits: " << its << "\n";
            std::cout << "p: " << p << "\n";
            std::cout << "fp: " << fp << "\n";
            std::cout << "grad_fp: " << xi << "\n";*/
            
            double test = 0.0;
            double den = std::max(fp, 1.0);
            
            for (size_t j = 0; j < n; j++)
            {
                double temp = std::abs(xi[j]) * std::max(std::abs(p[j]), 1.0) / den;
                
                if (temp > test)
                {
                    test = temp;
                }
            }
            
            if (test < GTOL)
            {
                break;
            }
            
            dgg = 0.0;
            gg = 0.0;
            
            for (size_t j = 0; j < n; j++)
            {
                gg += g[j] * g[j];
                //dgg += xi[j] * xi[j];
                dgg += (xi[j] + g[j]) * xi[j];
            }
            
            if (gg == 0.0)
            {
                break;
            }
            
            double gam = dgg / gg;
            
            for (size_t j = 0; j < n; j++)
            {
                g[j] = -xi[j];
                h[j] = g[j] + gam * h[j];
                xi[j] = h[j];
            }
        }
        
        //std::cerr << "Too many iterations in frprmn\n";
        
        /*if (verbosity_level >= 1)
        {
            Output::printStatusEnd();
        }*/
        
        std::stringstream msg;
        msg << its << " steps, f(x) = " << fp;
        tools::printer.printStatusUpdate(msg.str());
        
        return p;
    }
};

}
}
}
}

#endif
