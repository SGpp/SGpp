/* ****************************************************************************
* Copyright (C) 2009 Technische Universitaet Muenchen                         *
* This file is part of the SG++ project. For conditions of distribution and   *
* use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
**************************************************************************** */
// @author Dirk Pflueger (pflueged@in.tum.de), Julian Valentin (julian.valentin@stud.mathematik.uni-stuttgart.de)

#ifndef MODIFIED_BSPLINE_BASE_HPP
#define MODIFIED_BSPLINE_BASE_HPP

#include "base/basis/basis.hpp"
#include "base/basis/bspline/noboundary/BsplineBasis.hpp"

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
class ModifiedBsplineBasis : public Basis<LT, IT>
{
protected:
    BsplineBasis<LT, IT> bspline_basis;
    
public:
    ModifiedBsplineBasis() : bspline_basis(BsplineBasis<LT, IT>()) {}
    ModifiedBsplineBasis(size_t degree) : bspline_basis(BsplineBasis<LT, IT>(degree)) {}
    
    /*inline double modifiedBSpline(double x, size_t p, int level) const
    {
        double y = 0.0;
        double x2 = x * (1 << level) + static_cast<double>(p + 1) / 2.0 - 1.0;
        
        if (x2 > static_cast<double>(p) + 1.0)
        {
            return 0.0;
        }
        
        // the upper summation bound is defined to be ceil((p + 1) / 2.0),
        // which is the same as (p + 2) / 2 written in C
        for (size_t k = 0; k <= (p + 2) / 2; k++)
        {
            // x2 is chosen such that it holds: x2 = x * (1 << level) + (p+1)/2 + k - 1
            y += static_cast<double>(k + 1) * bspline_basis.uniformBSpline(x2, p);
            // the rounding errors induced by this method can be neglected
            x2 += 1.0;
        }
        
        return y;
    }
    
    inline double modifiedBSplineDx(double x, size_t p, int level) const
    {
        double y = 0.0;
        double x2 = x * (1 << level) + static_cast<double>(p + 1) / 2.0 - 1.0;
        
        if (x2 > static_cast<double>(p) + 1.0)
        {
            return 0.0;
        }
        
        for (size_t k = 0; k <= (p + 2) / 2; k++)
        {
            y += static_cast<double>(k + 1) * bspline_basis.uniformBSplineDx(x2, p);
            x2 += 1.0;
        }
        
        return y;
    }
    
    inline double modifiedBSplineDxDx(double x, size_t p, int level) const
    {
        double y = 0.0;
        double x2 = x * (1 << level) + static_cast<double>(p + 1) / 2.0 - 1.0;
        
        if (x2 > static_cast<double>(p) + 1.0)
        {
            return 0.0;
        }
        
        for (size_t k = 0; k <= (p + 2) / 2; k++)
        {
            y += static_cast<double>(k + 1) * bspline_basis.uniformBSplineDxDx(x2, p);
            x2 += 1.0;
        }
        
        return y;
    }*/
    
    inline double modifiedBSpline(double x, size_t p) const
    {
        switch (p)
        {
        case 1:
            if ((x < 0.0) || (x >= 2.0))
            {
                return 0.0;
            } else
            {
                return -x + 2.0;
            }
            break;
        case 2:
            if ((x < 0.0) || (x >= 5.0/2.0))
            {
                return 0.0;
            } else if (x < 3.0/2.0)
            {
                return -x + 2.0;
            } else
            {
                return 1.0/2.0*x*x - 5.0/2.0*x + 25.0/8.0;
            }
            break;
        case 3:
            if ((x < 0.0) || (x >= 3.0))
            {
                return 0.0;
            } else if (x < 1.0)
            {
                return -x + 2.0;
            } else if (x < 2.0)
            {
                return 1.0/6.0*x*x*x - 1.0/2.0*x*x - 1.0/2.0*x + 11.0/6.0;
            } else
            {
                return -1.0/6.0*x*x*x + 3.0/2.0*x*x - 9.0/2.0*x + 9.0/2.0;
            }
            break;
        case 4:
            if ((x < 0.0) || (x >= 7.0/2.0))
            {
                return 0.0;
            } else if (x < 1.0/2.0)
            {
                return -x + 2.0;
            } else if (x < 3.0/2.0)
            {
                double result = 1.0/24.0;
                result = -1.0/12.0 + result * x;
                result = 1.0/16.0 + result * x;
                result = -49.0/48.0 + result * x;
                result = 769.0/384.0 + result * x;
                return result;
            } else if (x < 5.0/2.0)
            {
                double result = -1.0/12.0;
                result = 2.0/3.0 + result * x;
                result = -13.0/8.0 + result * x;
                result = 2.0/3.0 + result * x;
                result = 263.0/192.0 + result * x;
                return result;
            } else
            {
                double result = 1.0/24.0;
                result = -7.0/12.0 + result * x;
                result = 49.0/16.0 + result * x;
                result = -343.0/48.0 + result * x;
                result = 2401.0/384.0 + result * x;
                return result;
            }
            break;
        case 5:
            if ((x < 0.0) || (x >= 4.0))
            {
                return 0.0;
            } else if (x < 1.0)
            {
                double result = 1.0/120.0;
                result *= x;
                result *= x;
                result *= x;
                result = -1.0 + result * x;
                result = 2.0 + result * x;
                return result;
            } else if (x < 2.0)
            {
                double result = -1.0/40.0;
                result = 1.0/6.0 + result * x;
                result = -1.0/3.0 + result * x;
                result = 1.0/3.0 + result * x;
                result = -7.0/6.0 + result * x;
                result = 61.0/30.0 + result * x;
                return result;
            } else if (x < 3.0)
            {
                double result = 1.0/40.0;
                result = -1.0/3.0 + result * x;
                result = 5.0/3.0 + result * x;
                result = -11.0/3.0 + result * x;
                result = 17.0/6.0 + result * x;
                result = 13.0/30.0 + result * x;
                return result;
            } else
            {
                double result = -1.0/120.0;
                result = 1.0/6.0 + result * x;
                result = -4.0/3.0 + result * x;
                result = 16.0/3.0 + result * x;
                result = -32.0/3.0 + result * x;
                result = 128.0/15.0 + result * x;
                return result;
            }
            break;
        case 6:
            if ((x < 0.0) || (x >= 9.0/2.0))
            {
                return 0.0;
            } else if (x < 1.0/2.0)
            {
                double result = 1.0/720.0;
                result = 1.0/240.0 + result * x;
                result = 1.0/192.0 + result * x;
                result = 1.0/288.0 + result * x;
                result = 1.0/768.0 + result * x;
                result = -3839.0/3840.0 + result * x;
                result = 92161.0/46080.0 + result * x;
                return result;
            } else if (x < 3.0/2.0)
            {
                double result = -1.0/180.0;
                result = 1.0/40.0 + result * x;
                result = -1.0/48.0 + result * x;
                result = 1.0/48.0 + result * x;
                result = -1.0/192.0 + result * x;
                result = -639.0/640.0 + result * x;
                result = 23039.0/11520.0 + result * x;
                return result;
            } else if (x < 5.0/2.0)
            {
                double result = 1.0/120.0;
                result = -1.0/10.0 + result * x;
                result = 43.0/96.0 + result * x;
                result = -11.0/12.0 + result * x;
                result = 403.0/384.0 + result * x;
                result = -261.0/160.0 + result * x;
                result = 49723.0/23040.0 + result * x;
                return result;
            } else if (x < 7.0/2.0)
            {
                double result = -1.0/180.0;
                result = 13.0/120.0 + result * x;
                result = -41.0/48.0 + result * x;
                result = 493.0/144.0 + result * x;
                result = -1361.0/192.0 + result * x;
                result = 12493.0/1920.0 + result * x;
                result = -14201.0/11520.0 + result * x;
                return result;
            } else
            {
                double result = 1.0/720.0;
                result = -3.0/80.0 + result * x;
                result = 27.0/64.0 + result * x;
                result = -81.0/32.0 + result * x;
                result = 2187.0/256.0 + result * x;
                result = -19683.0/1280.0 + result * x;
                result = 59049.0/5120.0 + result * x;
                return result;
            }
            break;
        case 7:
            if ((x < 0.0) || (x >= 5.0))
            {
                return 0.0;
            } else if (x < 1.0)
            {
                double result = -1.0/1008.0;
                result = 1.0/720.0 + result * x;
                result = 1.0/240.0 + result * x;
                result = 1.0/144.0 + result * x;
                result = 1.0/144.0 + result * x;
                result = 1.0/240.0 + result * x;
                result = -719.0/720.0 + result * x;
                result = 10081.0/5040.0 + result * x;
                return result;
            } else if (x < 2.0)
            {
                double result = 1.0/504.0;
                result = -7.0/360.0 + result * x;
                result = 1.0/15.0 + result * x;
                result = -7.0/72.0 + result * x;
                result = 1.0/9.0 + result * x;
                result = -7.0/120.0 + result * x;
                result = -44.0/45.0 + result * x;
                result = 719.0/360.0 + result * x;
                return result;
            } else if (x < 3.0)
            {
                double result = -1.0/504.0;
                result = 13.0/360.0 + result * x;
                result = -4.0/15.0 + result * x;
                result = 73.0/72.0 + result * x;
                result = -19.0/9.0 + result * x;
                result = 313.0/120.0 + result * x;
                result = -124.0/45.0 + result * x;
                result = 6313.0/2520.0 + result * x;
                return result;
            } else if (x < 4.0)
            {
                double result = 1.0/1008.0;
                result = -19.0/720.0 + result * x;
                result = 71.0/240.0 + result * x;
                result = -259.0/144.0 + result * x;
                result = 911.0/144.0 + result * x;
                result = -3019.0/240.0 + result * x;
                result = 8951.0/720.0 + result * x;
                result = -20179.0/5040.0 + result * x;
                return result;
            } else
            {
                double result = -1.0/5040.0;
                result = 1.0/144.0 + result * x;
                result = -5.0/48.0 + result * x;
                result = 125.0/144.0 + result * x;
                result = -625.0/144.0 + result * x;
                result = 625.0/48.0 + result * x;
                result = -3125.0/144.0 + result * x;
                result = 15625.0/1008.0 + result * x;
                return result;
            }
            break;
        case 8:
            if ((x < 0.0) || (x >= 11.0/2.0))
            {
                return 0.0;
            } else if (x < 1.0/2.0)
            {
                double result = -1.0/6720.0;
                result = -1.0/2520.0 + result * x;
                result = 1.0/2880.0 + result * x;
                result = 1.0/288.0 + result * x;
                result = 37.0/4608.0 + result * x;
                result = 59.0/5760.0 + result * x;
                result = 361.0/46080.0 + result * x;
                result = -32147.0/32256.0 + result * x;
                result = 10325197.0/5160960.0 + result * x;
                return result;
            } else if (x < 3.0/2.0)
            {
                double result = 1.0/2688.0;
                result = -5.0/2016.0 + result * x;
                result = 23.0/5760.0 + result * x;
                result = -1.0/5760.0 + result * x;
                result = 95.0/9216.0 + result * x;
                result = 43.0/4608.0 + result * x;
                result = 743.0/92160.0 + result * x;
                result = -642961.0/645120.0 + result * x;
                result = 4130083.0/2064384.0 + result * x;
                return result;
            } else if (x < 5.0/2.0)
            {
                double result = -1.0/2016.0;
                result = 1.0/126.0 + result * x;
                result = -73.0/1440.0 + result * x;
                result = 59.0/360.0 + result * x;
                result = -685.0/2304.0 + result * x;
                result = 109.0/288.0 + result * x;
                result = -6193.0/23040.0 + result * x;
                result = -35401.0/40320.0 + result * x;
                result = 1021039.0/516096.0 + result * x;
                return result;
            } else if (x < 7.0/2.0)
            {
                double result = 1.0/2688.0;
                result = -19.0/2016.0 + result * x;
                result = 583.0/5760.0 + result * x;
                result = -3431.0/5760.0 + result * x;
                result = 19135.0/9216.0 + result * x;
                result = -20131.0/4608.0 + result * x;
                result = 522103.0/92160.0 + result * x;
                result = -3300791.0/645120.0 + result * x;
                result = 6818531.0/2064384.0 + result * x;
                return result;
            } else if (x < 9.0/2.0)
            {
                double result = -1.0/6720.0;
                result = 13.0/2520.0 + result * x;
                result = -223.0/2880.0 + result * x;
                result = 943.0/1440.0 + result * x;
                result = -15643.0/4608.0 + result * x;
                result = 63073.0/5760.0 + result * x;
                result = -974263.0/46080.0 + result * x;
                result = 3498403.0/161280.0 + result * x;
                result = -43484083.0/5160960.0 + result * x;
                return result;
            } else
            {
                double result = 1.0/40320.0;
                result = -11.0/10080.0 + result * x;
                result = 121.0/5760.0 + result * x;
                result = -1331.0/5760.0 + result * x;
                result = 14641.0/9216.0 + result * x;
                result = -161051.0/23040.0 + result * x;
                result = 1771561.0/92160.0 + result * x;
                result = -19487171.0/645120.0 + result * x;
                result = 214358881.0/10321920.0 + result * x;
                return result;
            }
            break;
        case 9:
            if ((x < 0.0) || (x >= 6.0))
            {
                return 0.0;
            } else if (x < 1.0)
            {
                double result = 1.0/17280.0;
                result = -1.0/6720.0 + result * x;
                result = -1.0/2520.0 + result * x;
                result *= x;
                result = 1.0/360.0 + result * x;
                result = 1.0/120.0 + result * x;
                result = 7.0/540.0 + result * x;
                result = 1.0/84.0 + result * x;
                result = -5009.0/5040.0 + result * x;
                result = 1441.0/720.0 + result * x;
                return result;
            } else if (x < 2.0)
            {
                double result = -1.0/10368.0;
                result = 5.0/4032.0 + result * x;
                result = -1.0/168.0 + result * x;
                result = 7.0/540.0 + result * x;
                result = -1.0/60.0 + result * x;
                result = 1.0/36.0 + result * x;
                result *= x;
                result = 11.0/630.0 + result * x;
                result = -209.0/210.0 + result * x;
                result = 1297.0/648.0 + result * x;
                return result;
            } else if (x < 3.0)
            {
                double result = 1.0/10368.0;
                result = -1.0/448.0 + result * x;
                result = 11.0/504.0 + result * x;
                result = -7.0/60.0 + result * x;
                result = 67.0/180.0 + result * x;
                result = -3.0/4.0 + result * x;
                result = 28.0/27.0 + result * x;
                result = -61.0/70.0 + result * x;
                result = -347.0/630.0 + result * x;
                result = 137.0/72.0 + result * x;
                return result;
            } else if (x < 4.0)
            {
                double result = -1.0/17280.0;
                result = 13.0/6720.0 + result * x;
                result = -71.0/2520.0 + result * x;
                result = 7.0/30.0 + result * x;
                result = -433.0/360.0 + result * x;
                result = 159.0/40.0 + result * x;
                result = -4543.0/540.0 + result * x;
                result = 1579.0/140.0 + result * x;
                result = -48703.0/5040.0 + result * x;
                result = 3557.0/720.0 + result * x;
                return result;
            } else if (x < 5.0)
            {
                double result = 1.0/51840.0;
                result = -17.0/20160.0 + result * x;
                result = 41.0/2520.0 + result * x;
                result = -49.0/270.0 + result * x;
                result = 463.0/360.0 + result * x;
                result = -2153.0/360.0 + result * x;
                result = 9793.0/540.0 + result * x;
                result = -43133.0/1260.0 + result * x;
                result = 180673.0/5040.0 + result * x;
                result = -99059.0/6480.0 + result * x;
                return result;
            } else
            {
                double result = -1.0/362880.0;
                result = 1.0/6720.0 + result * x;
                result = -1.0/280.0 + result * x;
                result = 1.0/20.0 + result * x;
                result = -9.0/20.0 + result * x;
                result = 27.0/10.0 + result * x;
                result = -54.0/5.0 + result * x;
                result = 972.0/35.0 + result * x;
                result = -1458.0/35.0 + result * x;
                result = 972.0/35.0 + result * x;
                return result;
            }
            break;
        default:
            double y = 0.0;
            double x2 = x + static_cast<double>(p+1)/2.0 - 1.0;
            
            if (x2 > static_cast<double>(p) + 1.0)
            {
                return 0.0;
            }
            
            // the upper summation bound is defined to be ceil((p + 1) / 2.0),
            // which is the same as (p + 2) / 2 written in C
            for (size_t k = 0; k <= (p + 2) / 2; k++)
            {
                // x2 is chosen such that it holds: x2 = x + (p+1)/2 + k - 1
                y += static_cast<double>(k + 1) * bspline_basis.uniformBSpline(x2, p);
                // the rounding errors induced by this method can be neglected
                x2 += 1.0;
            }
            
            return y;
        }
    }

    inline double modifiedBSplineDx(double x, size_t p) const
    {
        switch (p)
        {
        case 1:
            if ((x < 0.0) || (x >= 2.0))
            {
                return 0.0;
            } else
            {
                return -1.0;
            }
            break;
        case 2:
            if ((x < 0.0) || (x >= 5.0/2.0))
            {
                return 0.0;
            } else if (x < 3.0/2.0)
            {
                return -1.0;
            } else
            {
                return x - 5.0/2.0;
            }
            break;
        case 3:
            if ((x < 0.0) || (x >= 3.0))
            {
                return 0.0;
            } else if (x < 1.0)
            {
                return -1.0;
            } else if (x < 2.0)
            {
                return 1.0/2.0*x*x - x - 1.0/2.0;
            } else
            {
                return -1.0/2.0*x*x + 3.0*x - 9.0/2.0;
            }
            break;
        case 4:
            if ((x < 0.0) || (x >= 7.0/2.0))
            {
                return 0.0;
            } else if (x < 1.0/2.0)
            {
                return -1.0;
            } else if (x < 3.0/2.0)
            {
                return 1.0/6.0*x*x*x - 1.0/4.0*x*x + 1.0/8.0*x - 49.0/48.0;
            } else if (x < 5.0/2.0)
            {
                return -1.0/3.0*x*x*x + 2.0*x*x - 13.0/4.0*x + 2.0/3.0;
            } else
            {
                return 1.0/6.0*x*x*x - 7.0/4.0*x*x + 49.0/8.0*x - 343.0/48.0;
            }
            break;
        case 5:
            if ((x < 0.0) || (x >= 4.0))
            {
                return 0.0;
            } else if (x < 1.0)
            {
                double result = 1.0/24.0;
                result *= x;
                result *= x;
                result *= x;
                result = -1.0 + result * x;
                return result;
            } else if (x < 2.0)
            {
                double result = -1.0/8.0;
                result = 2.0/3.0 + result * x;
                result = -1.0 + result * x;
                result = 2.0/3.0 + result * x;
                result = -7.0/6.0 + result * x;
                return result;
            } else if (x < 3.0)
            {
                double result = 1.0/8.0;
                result = -4.0/3.0 + result * x;
                result = 5.0 + result * x;
                result = -22.0/3.0 + result * x;
                result = 17.0/6.0 + result * x;
                return result;
            } else
            {
                double result = -1.0/24.0;
                result = 2.0/3.0 + result * x;
                result = -4.0 + result * x;
                result = 32.0/3.0 + result * x;
                result = -32.0/3.0 + result * x;
                return result;
            }
            break;
        case 6:
            if ((x < 0.0) || (x >= 9.0/2.0))
            {
                return 0.0;
            } else if (x < 1.0/2.0)
            {
                double result = 1.0/120.0;
                result = 1.0/48.0 + result * x;
                result = 1.0/48.0 + result * x;
                result = 1.0/96.0 + result * x;
                result = 1.0/384.0 + result * x;
                result = -3839.0/3840.0 + result * x;
                return result;
            } else if (x < 3.0/2.0)
            {
                double result = -1.0/30.0;
                result = 1.0/8.0 + result * x;
                result = -1.0/12.0 + result * x;
                result = 1.0/16.0 + result * x;
                result = -1.0/96.0 + result * x;
                result = -639.0/640.0 + result * x;
                return result;
            } else if (x < 5.0/2.0)
            {
                double result = 1.0/20.0;
                result = -1.0/2.0 + result * x;
                result = 43.0/24.0 + result * x;
                result = -11.0/4.0 + result * x;
                result = 403.0/192.0 + result * x;
                result = -261.0/160.0 + result * x;
                return result;
            } else if (x < 7.0/2.0)
            {
                double result = -1.0/30.0;
                result = 13.0/24.0 + result * x;
                result = -41.0/12.0 + result * x;
                result = 493.0/48.0 + result * x;
                result = -1361.0/96.0 + result * x;
                result = 12493.0/1920.0 + result * x;
                return result;
            } else
            {
                double result = 1.0/120.0;
                result = -3.0/16.0 + result * x;
                result = 27.0/16.0 + result * x;
                result = -243.0/32.0 + result * x;
                result = 2187.0/128.0 + result * x;
                result = -19683.0/1280.0 + result * x;
                return result;
            }
            break;
        case 7:
            if ((x < 0.0) || (x >= 5.0))
            {
                return 0.0;
            } else if (x < 1.0)
            {
                double result = -1.0/144.0;
                result = 1.0/120.0 + result * x;
                result = 1.0/48.0 + result * x;
                result = 1.0/36.0 + result * x;
                result = 1.0/48.0 + result * x;
                result = 1.0/120.0 + result * x;
                result = -719.0/720.0 + result * x;
                return result;
            } else if (x < 2.0)
            {
                double result = 1.0/72.0;
                result = -7.0/60.0 + result * x;
                result = 1.0/3.0 + result * x;
                result = -7.0/18.0 + result * x;
                result = 1.0/3.0 + result * x;
                result = -7.0/60.0 + result * x;
                result = -44.0/45.0 + result * x;
                return result;
            } else if (x < 3.0)
            {
                double result = -1.0/72.0;
                result = 13.0/60.0 + result * x;
                result = -4.0/3.0 + result * x;
                result = 73.0/18.0 + result * x;
                result = -19.0/3.0 + result * x;
                result = 313.0/60.0 + result * x;
                result = -124.0/45.0 + result * x;
                return result;
            } else if (x < 4.0)
            {
                double result = 1.0/144.0;
                result = -19.0/120.0 + result * x;
                result = 71.0/48.0 + result * x;
                result = -259.0/36.0 + result * x;
                result = 911.0/48.0 + result * x;
                result = -3019.0/120.0 + result * x;
                result = 8951.0/720.0 + result * x;
                return result;
            } else
            {
                double result = -1.0/720.0;
                result = 1.0/24.0 + result * x;
                result = -25.0/48.0 + result * x;
                result = 125.0/36.0 + result * x;
                result = -625.0/48.0 + result * x;
                result = 625.0/24.0 + result * x;
                result = -3125.0/144.0 + result * x;
                return result;
            }
            break;
        case 8:
            if ((x < 0.0) || (x >= 11.0/2.0))
            {
                return 0.0;
            } else if (x < 1.0/2.0)
            {
                double result = -1.0/840.0;
                result = -1.0/360.0 + result * x;
                result = 1.0/480.0 + result * x;
                result = 5.0/288.0 + result * x;
                result = 37.0/1152.0 + result * x;
                result = 59.0/1920.0 + result * x;
                result = 361.0/23040.0 + result * x;
                result = -32147.0/32256.0 + result * x;
                return result;
            } else if (x < 3.0/2.0)
            {
                double result = 1.0/336.0;
                result = -5.0/288.0 + result * x;
                result = 23.0/960.0 + result * x;
                result = -1.0/1152.0 + result * x;
                result = 95.0/2304.0 + result * x;
                result = 43.0/1536.0 + result * x;
                result = 743.0/46080.0 + result * x;
                result = -642961.0/645120.0 + result * x;
                return result;
            } else if (x < 5.0/2.0)
            {
                double result = -1.0/252.0;
                result = 1.0/18.0 + result * x;
                result = -73.0/240.0 + result * x;
                result = 59.0/72.0 + result * x;
                result = -685.0/576.0 + result * x;
                result = 109.0/96.0 + result * x;
                result = -6193.0/11520.0 + result * x;
                result = -35401.0/40320.0 + result * x;
                return result;
            } else if (x < 7.0/2.0)
            {
                double result = 1.0/336.0;
                result = -19.0/288.0 + result * x;
                result = 583.0/960.0 + result * x;
                result = -3431.0/1152.0 + result * x;
                result = 19135.0/2304.0 + result * x;
                result = -20131.0/1536.0 + result * x;
                result = 522103.0/46080.0 + result * x;
                result = -3300791.0/645120.0 + result * x;
                return result;
            } else if (x < 9.0/2.0)
            {
                double result = -1.0/840.0;
                result = 13.0/360.0 + result * x;
                result = -223.0/480.0 + result * x;
                result = 943.0/288.0 + result * x;
                result = -15643.0/1152.0 + result * x;
                result = 63073.0/1920.0 + result * x;
                result = -974263.0/23040.0 + result * x;
                result = 3498403.0/161280.0 + result * x;
                return result;
            } else
            {
                double result = 1.0/5040.0;
                result = -11.0/1440.0 + result * x;
                result = 121.0/960.0 + result * x;
                result = -1331.0/1152.0 + result * x;
                result = 14641.0/2304.0 + result * x;
                result = -161051.0/7680.0 + result * x;
                result = 1771561.0/46080.0 + result * x;
                result = -19487171.0/645120.0 + result * x;
                return result;
            }
            break;
        case 9:
            if ((x < 0.0) || (x >= 6.0))
            {
                return 0.0;
            } else if (x < 1.0)
            {
                double result = 1.0/1920.0;
                result = -1.0/840.0 + result * x;
                result = -1.0/360.0 + result * x;
                result *= x;
                result = 1.0/72.0 + result * x;
                result = 1.0/30.0 + result * x;
                result = 7.0/180.0 + result * x;
                result = 1.0/42.0 + result * x;
                result = -5009.0/5040.0 + result * x;
                return result;
            } else if (x < 2.0)
            {
                double result = -1.0/1152.0;
                result = 5.0/504.0 + result * x;
                result = -1.0/24.0 + result * x;
                result = 7.0/90.0 + result * x;
                result = -1.0/12.0 + result * x;
                result = 1.0/9.0 + result * x;
                result *= x;
                result = 11.0/315.0 + result * x;
                result = -209.0/210.0 + result * x;
                return result;
            } else if (x < 3.0)
            {
                double result = 1.0/1152.0;
                result = -1.0/56.0 + result * x;
                result = 11.0/72.0 + result * x;
                result = -7.0/10.0 + result * x;
                result = 67.0/36.0 + result * x;
                result = -3.0 + result * x;
                result = 28.0/9.0 + result * x;
                result = -61.0/35.0 + result * x;
                result = -347.0/630.0 + result * x;
                return result;
            } else if (x < 4.0)
            {
                double result = -1.0/1920.0;
                result = 13.0/840.0 + result * x;
                result = -71.0/360.0 + result * x;
                result = 7.0/5.0 + result * x;
                result = -433.0/72.0 + result * x;
                result = 159.0/10.0 + result * x;
                result = -4543.0/180.0 + result * x;
                result = 1579.0/70.0 + result * x;
                result = -48703.0/5040.0 + result * x;
                return result;
            } else if (x < 5.0)
            {
                double result = 1.0/5760.0;
                result = -17.0/2520.0 + result * x;
                result = 41.0/360.0 + result * x;
                result = -49.0/45.0 + result * x;
                result = 463.0/72.0 + result * x;
                result = -2153.0/90.0 + result * x;
                result = 9793.0/180.0 + result * x;
                result = -43133.0/630.0 + result * x;
                result = 180673.0/5040.0 + result * x;
                return result;
            } else
            {
                double result = -1.0/40320.0;
                result = 1.0/840.0 + result * x;
                result = -1.0/40.0 + result * x;
                result = 3.0/10.0 + result * x;
                result = -9.0/4.0 + result * x;
                result = 54.0/5.0 + result * x;
                result = -162.0/5.0 + result * x;
                result = 1944.0/35.0 + result * x;
                result = -1458.0/35.0 + result * x;
                return result;
            }
            break;
        default:
            double y = 0.0;
            double x2 = x + static_cast<double>(p+1)/2.0 - 1.0;
            
            if (x2 > static_cast<double>(p) + 1.0)
            {
                return 0.0;
            }
            
            for (size_t k = 0; k <= (p + 2) / 2; k++)
            {
                y += static_cast<double>(k + 1) * bspline_basis.uniformBSplineDx(x2, p);
                x2 += 1.0;
            }
            
            return y;
        }
    }

    inline double modifiedBSplineDxDx(double x, size_t p) const
    {
        switch (p)
        {
        case 1:
            return 0.0;
            break;
        case 2:
            if ((x < 0.0) || (x >= 5.0/2.0))
            {
                return 0.0;
            } else if (x < 3.0/2.0)
            {
                return 0.0;
            } else
            {
                return 1.0;
            }
            break;
        case 3:
            if ((x < 0.0) || (x >= 3.0))
            {
                return 0.0;
            } else if (x < 1.0)
            {
                return 0.0;
            } else if (x < 2.0)
            {
                return x - 1.0;
            } else
            {
                return -x + 3.0;
            }
            break;
        case 4:
            if ((x < 0.0) || (x >= 7.0/2.0))
            {
                return 0.0;
            } else if (x < 1.0/2.0)
            {
                return 0.0;
            } else if (x < 3.0/2.0)
            {
                return 1.0/2.0*x*x - 1.0/2.0*x + 1.0/8.0;
            } else if (x < 5.0/2.0)
            {
                return -x*x + 4.0*x - 13.0/4.0;
            } else
            {
                return 1.0/2.0*x*x - 7.0/2.0*x + 49.0/8.0;
            }
            break;
        case 5:
            if ((x < 0.0) || (x >= 4.0))
            {
                return 0.0;
            } else if (x < 1.0)
            {
                return 1.0/6.0*x*x*x;
            } else if (x < 2.0)
            {
                return -1.0/2.0*x*x*x + 2.0*x*x - 2.0*x + 2.0/3.0;
            } else if (x < 3.0)
            {
                return 1.0/2.0*x*x*x - 4.0*x*x + 10.0*x - 22.0/3.0;
            } else
            {
                return -1.0/6.0*x*x*x + 2.0*x*x - 8.0*x + 32.0/3.0;
            }
            break;
        case 6:
            if ((x < 0.0) || (x >= 9.0/2.0))
            {
                return 0.0;
            } else if (x < 1.0/2.0)
            {
                double result = 1.0/24.0;
                result = 1.0/12.0 + result * x;
                result = 1.0/16.0 + result * x;
                result = 1.0/48.0 + result * x;
                result = 1.0/384.0 + result * x;
                return result;
            } else if (x < 3.0/2.0)
            {
                double result = -1.0/6.0;
                result = 1.0/2.0 + result * x;
                result = -1.0/4.0 + result * x;
                result = 1.0/8.0 + result * x;
                result = -1.0/96.0 + result * x;
                return result;
            } else if (x < 5.0/2.0)
            {
                double result = 1.0/4.0;
                result = -2.0 + result * x;
                result = 43.0/8.0 + result * x;
                result = -11.0/2.0 + result * x;
                result = 403.0/192.0 + result * x;
                return result;
            } else if (x < 7.0/2.0)
            {
                double result = -1.0/6.0;
                result = 13.0/6.0 + result * x;
                result = -41.0/4.0 + result * x;
                result = 493.0/24.0 + result * x;
                result = -1361.0/96.0 + result * x;
                return result;
            } else
            {
                double result = 1.0/24.0;
                result = -3.0/4.0 + result * x;
                result = 81.0/16.0 + result * x;
                result = -243.0/16.0 + result * x;
                result = 2187.0/128.0 + result * x;
                return result;
            }
            break;
        case 7:
            if ((x < 0.0) || (x >= 5.0))
            {
                return 0.0;
            } else if (x < 1.0)
            {
                double result = -1.0/24.0;
                result = 1.0/24.0 + result * x;
                result = 1.0/12.0 + result * x;
                result = 1.0/12.0 + result * x;
                result = 1.0/24.0 + result * x;
                result = 1.0/120.0 + result * x;
                return result;
            } else if (x < 2.0)
            {
                double result = 1.0/12.0;
                result = -7.0/12.0 + result * x;
                result = 4.0/3.0 + result * x;
                result = -7.0/6.0 + result * x;
                result = 2.0/3.0 + result * x;
                result = -7.0/60.0 + result * x;
                return result;
            } else if (x < 3.0)
            {
                double result = -1.0/12.0;
                result = 13.0/12.0 + result * x;
                result = -16.0/3.0 + result * x;
                result = 73.0/6.0 + result * x;
                result = -38.0/3.0 + result * x;
                result = 313.0/60.0 + result * x;
                return result;
            } else if (x < 4.0)
            {
                double result = 1.0/24.0;
                result = -19.0/24.0 + result * x;
                result = 71.0/12.0 + result * x;
                result = -259.0/12.0 + result * x;
                result = 911.0/24.0 + result * x;
                result = -3019.0/120.0 + result * x;
                return result;
            } else
            {
                double result = -1.0/120.0;
                result = 5.0/24.0 + result * x;
                result = -25.0/12.0 + result * x;
                result = 125.0/12.0 + result * x;
                result = -625.0/24.0 + result * x;
                result = 625.0/24.0 + result * x;
                return result;
            }
            break;
        case 8:
            if ((x < 0.0) || (x >= 11.0/2.0))
            {
                return 0.0;
            } else if (x < 1.0/2.0)
            {
                double result = -1.0/120.0;
                result = -1.0/60.0 + result * x;
                result = 1.0/96.0 + result * x;
                result = 5.0/72.0 + result * x;
                result = 37.0/384.0 + result * x;
                result = 59.0/960.0 + result * x;
                result = 361.0/23040.0 + result * x;
                return result;
            } else if (x < 3.0/2.0)
            {
                double result = 1.0/48.0;
                result = -5.0/48.0 + result * x;
                result = 23.0/192.0 + result * x;
                result = -1.0/288.0 + result * x;
                result = 95.0/768.0 + result * x;
                result = 43.0/768.0 + result * x;
                result = 743.0/46080.0 + result * x;
                return result;
            } else if (x < 5.0/2.0)
            {
                double result = -1.0/36.0;
                result = 1.0/3.0 + result * x;
                result = -73.0/48.0 + result * x;
                result = 59.0/18.0 + result * x;
                result = -685.0/192.0 + result * x;
                result = 109.0/48.0 + result * x;
                result = -6193.0/11520.0 + result * x;
                return result;
            } else if (x < 7.0/2.0)
            {
                double result = 1.0/48.0;
                result = -19.0/48.0 + result * x;
                result = 583.0/192.0 + result * x;
                result = -3431.0/288.0 + result * x;
                result = 19135.0/768.0 + result * x;
                result = -20131.0/768.0 + result * x;
                result = 522103.0/46080.0 + result * x;
                return result;
            } else if (x < 9.0/2.0)
            {
                double result = -1.0/120.0;
                result = 13.0/60.0 + result * x;
                result = -223.0/96.0 + result * x;
                result = 943.0/72.0 + result * x;
                result = -15643.0/384.0 + result * x;
                result = 63073.0/960.0 + result * x;
                result = -974263.0/23040.0 + result * x;
                return result;
            } else
            {
                double result = 1.0/720.0;
                result = -11.0/240.0 + result * x;
                result = 121.0/192.0 + result * x;
                result = -1331.0/288.0 + result * x;
                result = 14641.0/768.0 + result * x;
                result = -161051.0/3840.0 + result * x;
                result = 1771561.0/46080.0 + result * x;
                return result;
            }
            break;
        case 9:
            if ((x < 0.0) || (x >= 6.0))
            {
                return 0.0;
            } else if (x < 1.0)
            {
                double result = 1.0/240.0;
                result = -1.0/120.0 + result * x;
                result = -1.0/60.0 + result * x;
                result *= x;
                result = 1.0/18.0 + result * x;
                result = 1.0/10.0 + result * x;
                result = 7.0/90.0 + result * x;
                result = 1.0/42.0 + result * x;
                return result;
            } else if (x < 2.0)
            {
                double result = -1.0/144.0;
                result = 5.0/72.0 + result * x;
                result = -1.0/4.0 + result * x;
                result = 7.0/18.0 + result * x;
                result = -1.0/3.0 + result * x;
                result = 1.0/3.0 + result * x;
                result *= x;
                result = 11.0/315.0 + result * x;
                return result;
            } else if (x < 3.0)
            {
                double result = 1.0/144.0;
                result = -1.0/8.0 + result * x;
                result = 11.0/12.0 + result * x;
                result = -7.0/2.0 + result * x;
                result = 67.0/9.0 + result * x;
                result = -9.0 + result * x;
                result = 56.0/9.0 + result * x;
                result = -61.0/35.0 + result * x;
                return result;
            } else if (x < 4.0)
            {
                double result = -1.0/240.0;
                result = 13.0/120.0 + result * x;
                result = -71.0/60.0 + result * x;
                result = 7.0 + result * x;
                result = -433.0/18.0 + result * x;
                result = 477.0/10.0 + result * x;
                result = -4543.0/90.0 + result * x;
                result = 1579.0/70.0 + result * x;
                return result;
            } else if (x < 5.0)
            {
                double result = 1.0/720.0;
                result = -17.0/360.0 + result * x;
                result = 41.0/60.0 + result * x;
                result = -49.0/9.0 + result * x;
                result = 463.0/18.0 + result * x;
                result = -2153.0/30.0 + result * x;
                result = 9793.0/90.0 + result * x;
                result = -43133.0/630.0 + result * x;
                return result;
            } else
            {
                double result = -1.0/5040.0;
                result = 1.0/120.0 + result * x;
                result = -3.0/20.0 + result * x;
                result = 3.0/2.0 + result * x;
                result = -9.0 + result * x;
                result = 162.0/5.0 + result * x;
                result = -324.0/5.0 + result * x;
                result = 1944.0/35.0 + result * x;
                return result;
            }
            break;
        default:
            double y = 0.0;
            double x2 = x + static_cast<double>(p+1)/2.0 - 1.0;
            
            if (x2 > static_cast<double>(p) + 1.0)
            {
                return 0.0;
            }
            
            for (size_t k = 0; k <= (p + 2) / 2; k++)
            {
                y += static_cast<double>(k + 1) * bspline_basis.uniformBSplineDxDx(x2, p);
                x2 += 1.0;
            }
            
            return y;
        }
    }

    inline double eval(LT level, IT index, double x)
    {
        //std::cout << "Level " << level <<" Index "<<index<<" Point "<<p<<" BasisValue ";
        if (static_cast<int>(level) == 1)
        {
            //std::cout<<1<<std::endl;
            return 1.0;
        }
        
        double hinv = static_cast<double>(1 << level);
        
        if (static_cast<int>(index) == static_cast<int>(1 << level) - 1)
        {
            x = 1.0 - x;
            index = 1;
        }
        
        if (index == 1)
        {
            //std::cout<<BoundaryBSpline(p,this->degree,level)<<std::endl;
            //return modifiedBSpline(x, bspline_basis.getDegree(), level);
            return modifiedBSpline(x * hinv, bspline_basis.getDegree());
        } else
        {
            //std::cout<<UniformBSpline(p*(1<<level)+(this->degree+1)/2-index,this->degree)<<std::endl;
            return bspline_basis.uniformBSpline(
                    x * hinv + static_cast<double>(bspline_basis.getDegree() + 1) / 2 - index,
                    bspline_basis.getDegree());
        }
    }
    
    inline double evalDx(LT level, IT index, double x)
    {
        if (static_cast<int>(level) == 1)
        {
            return 0.0;
        }
        
        double hinv = static_cast<double>(1 << level);
        // inner derivative
        double dx_factor = hinv;
        
        if (static_cast<int>(index) == static_cast<int>(1 << level) - 1)
        {
            x = 1.0 - x;
            index = 1;
            dx_factor *= -1.0;
        }
        
        if (index == 1)
        {
            //return dx_factor * modifiedBSplineDx(x, bspline_basis.getDegree(), level);
            return dx_factor * modifiedBSplineDx(x * hinv, bspline_basis.getDegree());
        } else
        {
            return dx_factor * bspline_basis.uniformBSplineDx(
                    x * hinv + static_cast<double>(bspline_basis.getDegree() + 1) / 2 - index,
                    bspline_basis.getDegree());
        }
    }
    
    inline double evalDxDx(LT level, IT index, double x)
    {
        if (static_cast<int>(level) == 1)
        {
            return 0.0;
        }
        
        double hinv = static_cast<double>(1 << level);
        // inner derivative
        double dx_factor = hinv;
        
        if (static_cast<int>(index) == static_cast<int>(1 << level) - 1)
        {
            x = 1.0 - x;
            index = 1;
        }
        
        if (index == 1)
        {
            //return dx_factor*dx_factor * modifiedBSplineDxDx(x, bspline_basis.getDegree(), level);
            return dx_factor*dx_factor * modifiedBSplineDxDx(x * hinv, bspline_basis.getDegree());
        } else
        {
            return dx_factor*dx_factor * bspline_basis.uniformBSplineDxDx(
                    x * hinv + static_cast<double>(bspline_basis.getDegree() + 1) / 2 - index,
                    bspline_basis.getDegree());
        }
    }
    
    inline size_t getDegree() const { return bspline_basis.getDegree(); }
};

// default type-def (unsigned int for level and index)
typedef ModifiedBsplineBasis<unsigned int, unsigned int> SModBsplineBase;

}
}

#endif /* MODIFIED_BSPLINE_BASE_HPP */
