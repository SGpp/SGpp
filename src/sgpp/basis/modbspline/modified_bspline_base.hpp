/******************************************************************************
* Copyright (C) 2009 Technische Universitaet Muenchen                         *
* This file is part of the SG++ project. For conditions of distribution and   *
* use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
******************************************************************************/
// @author Dirk Pflueger (pflueged@in.tum.de)

#ifndef MODIFIED_BSPLINE_BASE_HPP
#define MODIFIED_BSPLINE_BASE_HPP

#include <cmath>

namespace sg
{
namespace base
{

/**
 * Modified B-Spline basis functions.
 *
 */
template<class LT, class IT>
class modified_bspline_base
{
    protected:
        size_t degree;


    public:
        modified_bspline_base():degree(0){}
        modified_bspline_base(size_t degree) : degree(degree){
            if (degree < 1)
                this->degree=1;
            if(degree%2 == 0)
                this->degree = degree-1;

        }
        double BoundaryBSpline(double x, int p, int level){
            double y=0;
            for(int k=0; k <= (p+1)/2; k++)
                y=y+(k+1)*UniformBSpline(x*(1<<level)+(p-1)/2+k,p);
            return y;
        }

        double UniformBSpline(double x,int p)
        {
            if (p==0)
            { if ( (0 < x) && (x <=1 ) )
                return 1;
              else
                return 0;
            }
            return (x/p)*UniformBSpline(x,p-1)+ ((p+1-x)/p)*UniformBSpline(x-1,p-1);

        }

        double eval(LT level, IT index, double p)
        {   //std::cout << "Level " << level <<" Index "<<index<<" Point "<<p<<" BasisValue ";
            if (static_cast<int>(level)==1)
            {
                //std::cout<<1<<std::endl;
                return 1;
            }
            if (static_cast<int>(index)==static_cast<int>(1<<level)-1)
            {
                p=1-p;
                index=1;
            }
            if (index==1)
            {   //std::cout<<BoundaryBSpline(p,this->degree,level)<<std::endl;
                return BoundaryBSpline(p,this->degree,level);
            }
            //std::cout<<UniformBSpline(p*(1<<level)+(this->degree+1)/2-index,this->degree)<<std::endl;
            return UniformBSpline(p*(1<<level)+(this->degree+1)/2-index,this->degree);
        }

};

}
}

#endif /* MODIFIED_BSPLINE_BASE_HPP */
