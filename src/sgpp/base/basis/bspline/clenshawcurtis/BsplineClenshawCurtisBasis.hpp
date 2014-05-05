/* ****************************************************************************
* Copyright (C) 2009 Technische Universitaet Muenchen                         *
* This file is part of the SG++ project. For conditions of distribution and   *
* use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
**************************************************************************** */
// @author Julian Valentin (julian.valentin@stud.mathematik.uni-stuttgart.de)

#ifndef BSPLINE_CLENSHAWCURTIS_BASE_HPP
#define BSPLINE_CLENSHAWCURTIS_BASE_HPP

#include "base/basis/basis.hpp"
#include "base/basis/bspline/noboundary/BsplineBasis.hpp"
#include "base/tools/CosineTable.hpp"

#include <cmath>

namespace sg
{
namespace base
{

template<class LT, class IT>
class BsplineClenshawCurtisBasis : public Basis<LT, IT>
{
protected:
    BsplineBasis<LT, IT> bspline_basis;
    const CosineTable *cosine_table;
    std::vector<double> knots;
    
public:
    BsplineClenshawCurtisBasis() : bspline_basis(BsplineBasis<LT, IT>()) {}
    BsplineClenshawCurtisBasis(size_t degree) : BsplineClenshawCurtisBasis(degree, NULL) {}
    BsplineClenshawCurtisBasis(size_t degree, const CosineTable *cosine_table) :
            bspline_basis(BsplineBasis<LT, IT>(degree)),
            cosine_table(cosine_table),
            knots(std::vector<double>(degree+2, 0.0)) {}
    
    inline double nonUniformBSpline(double x, size_t p, size_t k,
                                    const std::vector<double> &xi) const
    {
        if (p == 0)
        {
            return (((x >= xi[k]) && (x < xi[k+1])) ? 1.0 : 0.0);
        } else if ((x < xi[k]) || (x >= xi[k+p+1]))
        {
            return 0.0;
        } else
        {
            return (x - xi[k]) / (xi[k+p] - xi[k]) * nonUniformBSpline(x, p-1, k, xi)
                    + (1.0 - (x - xi[k+1]) / (xi[k+p+1] - xi[k+1]))
                      * nonUniformBSpline(x, p-1, k+1, xi);
        }
    }
    
    inline double nonUniformBSplineDx(double x, size_t p, size_t k,
                                      const std::vector<double> &xi) const
    {
        if (p == 0)
        {
            return 0.0;
        } else if ((x < xi[k]) || (x >= xi[k+p+1]))
        {
            return 0.0;
        } else
        {
            double p_double = static_cast<double>(p);
            return p_double / (xi[k+p] - xi[k]) * nonUniformBSpline(x, p-1, k, xi)
                    - p_double / (xi[k+p+1] - xi[k+1]) * nonUniformBSpline(x, p-1, k+1, xi);
        }
    }
    
    inline double nonUniformBSplineDxDx(double x, size_t p, size_t k,
                                        const std::vector<double> &xi) const
    {
        if (p <= 1)
        {
            return 0.0;
        } else if ((x < xi[k]) || (x >= xi[k+p+1]))
        {
            return 0.0;
        } else
        {
            double p_double = static_cast<double>(p);
            double alpha_k_p = p_double / (xi[k+p] - xi[k]);
            double alpha_kp1_p = p_double / (xi[k+p+1] - xi[k+1]);
            double alpha_k_pm1 = (p_double - 1.0) / (xi[k+p-1] - xi[k]);
            double alpha_kp1_pm1 = (p_double - 1.0) / (xi[k+p] - xi[k+1]);
            double alpha_kp2_pm1 = (p_double - 1.0) / (xi[k+p+1] - xi[k+2]);
            
            return alpha_k_p * alpha_k_pm1 * nonUniformBSpline(x, p-2, k, xi)
                    - (alpha_k_p + alpha_kp1_p) * alpha_kp1_pm1 *
                      nonUniformBSpline(x, p-2, k+1, xi)
                    - alpha_kp1_p * alpha_kp2_pm1 * nonUniformBSpline(x, p-2, k+2, xi);
        }
    }
    
    inline double clenshawCurtisPoint(double h, IT index) const
    {
        /*double arg = M_PI * (1.0 - static_cast<double>(index) /
                             static_cast<double>(1 << level));
        if ((arg < 0.01) && (arg > 0.0))
        {
            std::cout << "cos(" << arg << ")\n";
        }*/
        
        // h = 1.0 / static_cast<double>(1 << level)
        if (cosine_table == NULL)
        {
            return (cos(M_PI * (1.0 - static_cast<double>(index) * h)) + 1.0) / 2.0;
        } else
        {
            return (cosine_table->lookUp(
                    M_PI * (1.0 - static_cast<double>(index) * h)) + 1.0) / 2.0;
            /*size_t hash = 0;
            while ((index % 2 == 0) && (level > 0))
            {
                level--;
                index /= 2;
            }
            boost::hash_combine(hash, level);
            boost::hash_combine(hash, index);
            
            double result = (cosine_table->lookUp(M_PI * (1.0 - static_cast<double>(index) /
                                         static_cast<double>(1 << level))) + 1.0) / 2.0;
            auto a = hash_map.find(result);
            if (a == hash_map.end())
            {
                hash_map[result] = hash;
                hash_map_l[result] = level;
                hash_map_i[result] = index;
            } else
            {
                if (a->second != hash)
                {
                    std::cout << "BOOM! l = " << level << ", i = " << index << ", other_l = " << hash_map_l[result] << ", other_i = " << hash_map_i[result] << ", result = " << result << "\n";
                }
            }
            return result;*/
        }
    }
    
    inline void constructKnots(LT level, IT index)
    {
        IT n = 1 << level;
        double h = 1.0 / static_cast<double>(n);
        size_t p = bspline_basis.getDegree();
        knots[(p+1)/2] = clenshawCurtisPoint(h, index);
        
        if (index < (p+1)/2)
        {
            size_t a = (p+1)/2 - index;
            
            for (size_t i = a; i < (p+1)/2; i++)
            {
                knots[i] = clenshawCurtisPoint(h, static_cast<IT>(i - a));
            }
            
            double h = knots[a+1] - knots[a];
            
            // equivalent to "for (int i = a-1; i >= 0; i--)"
            for (size_t i = a; i-- > 0; )
            {
                knots[i] = knots[i+1] - h;
            }
        } else
        {
            for (size_t i = 0; i < (p+1)/2; i++)
            {
                knots[i] = clenshawCurtisPoint(h, static_cast<IT>(index - (p+1)/2 + i));
            }
        }
        
        if (index + (p+1)/2 > n)
        {
            size_t b = n + (p+1)/2 - index;
            
            for (size_t i = (p+1)/2 + 1; i <= b; i++)
            {
                knots[i] = clenshawCurtisPoint(h, static_cast<IT>(index - (p+1)/2 + i));
                //std::cout << "clenshawCurtisPoint(" << level << ", " << static_cast<IT>(index - (p+1)/2 + i) << ") = " << knots[i] << "\n";
            }
            
            double h = knots[b] - knots[b-1];
            
            for (size_t i = b+1; i < p+2; i++)
            {
                knots[i] = knots[i-1] + h;
            }
        } else
        {
            for (size_t i = (p+1)/2 + 1; i < p+2; i++)
            {
                knots[i] = clenshawCurtisPoint(h, static_cast<IT>(index - (p+1)/2 + i));
            }
        }
    }
    
    inline double eval(LT level, IT index, double x)
    {
        //std::cout << "eval: l = " << level << ", i = " << index << ", x = " << x << "\n";
        
        if (level == 0)
        {
            return bspline_basis.uniformBSpline(
                    x + static_cast<double>(bspline_basis.getDegree() + 1) / 2 - index,
                    bspline_basis.getDegree());
        } else
        {
            constructKnots(level, index);
            
            /*if (level == 1)
            {*/
                /*std::cout << "eval: l = " << level << ", i = " << index << ", x = " << x << "\n";
                std::cout << "knots = [";
                for (size_t i = 0; i < knots.size(); i++)
                {
                    if (i > 0) std::cout << ", ";
                    std::cout << knots[i];
                }
                std::cout << "]\n";*/
            //}
            
            return nonUniformBSpline(x, bspline_basis.getDegree(), 0, knots);
        }
    }
    
    inline double evalDx(LT level, IT index, double x)
    {
        if (level == 0)
        {
            return bspline_basis.uniformBSplineDx(
                    x + static_cast<double>(bspline_basis.getDegree() + 1) / 2 - index,
                    bspline_basis.getDegree());
        } else
        {
            constructKnots(level, index);
            return nonUniformBSplineDx(x, bspline_basis.getDegree(), 0, knots);
        }
    }
    
    inline double evalDxDx(LT level, IT index, double x)
    {
        if (level == 0)
        {
            return bspline_basis.uniformBSplineDxDx(
                    x + static_cast<double>(bspline_basis.getDegree() + 1) / 2 - index,
                    bspline_basis.getDegree());
        } else
        {
            constructKnots(level, index);
            return nonUniformBSplineDxDx(x, bspline_basis.getDegree(), 0, knots);
        }
    }
    
    inline size_t getDegree() const { return bspline_basis.degree; }
};

typedef BsplineClenshawCurtisBasis<unsigned int, unsigned int> SBsplineClenshawCurtisBase;

}
}

#endif /* BSPLINE_CLENSHAWCURTIS_BASE_HPP */
