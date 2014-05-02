/* ****************************************************************************
* Copyright (C) 2009 Technische Universitaet Muenchen                         *
* This file is part of the SG++ project. For conditions of distribution and   *
* use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
**************************************************************************** */
// @author Julian Valentin (julian.valentin@stud.mathematik.uni-stuttgart.de)

#ifndef BSPLINE_BASE_HPP
#define BSPLINE_BASE_HPP

#include <cmath>

namespace sg
{
namespace base
{

template<class LT, class IT>
class BsplineBasis
{
protected:
    size_t degree;
    
public:
    BsplineBasis(): degree(0) {}
    
    BsplineBasis(size_t degree) : degree(degree)
    {
        if (degree < 1) {
            this->degree = 1;
        }
        
        if (degree % 2 == 0)
        {
            this->degree = degree - 1;
        }
    }
    
    inline double uniformBSpline(double x, size_t p) const
    {
        switch (p)
        {
        case 0:
            if ((x < 0.0) || (x >= 1.0))
            {
                return 0.0;
            } else
            {
                return 1.0;
            }
            break;
        case 1:
            if ((x < 0.0) || (x >= 2.0))
            {
                return 0.0;
            } else if (x < 1.0)
            {
                return x;
            } else
            {
                return -x + 2.0;
            }
            break;
        case 2:
            if ((x < 0.0) || (x >= 3.0))
            {
                return 0.0;
            } else if (x < 1.0)
            {
                return 1.0/2.0*x*x;
            } else if (x < 2.0)
            {
                return -x*x + 3.0*x - 3.0/2.0;
            } else
            {
                return 1.0/2.0*x*x - 3.0*x + 9.0/2.0;
            }
            break;
        case 3:
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
        case 4:
            if ((x < 0.0) || (x >= 5.0))
            {
                return 0.0;
            } else if (x < 1.0)
            {
                double result = 1.0/24.0;
                result *= x;
                result *= x;
                result *= x;
                result *= x;
                return result;
            } else if (x < 2.0)
            {
                double result = -1.0/6.0;
                result = 5.0/6.0 + result * x;
                result = -5.0/4.0 + result * x;
                result = 5.0/6.0 + result * x;
                result = -5.0/24.0 + result * x;
                return result;
            } else if (x < 3.0)
            {
                double result = 1.0/4.0;
                result = -5.0/2.0 + result * x;
                result = 35.0/4.0 + result * x;
                result = -25.0/2.0 + result * x;
                result = 155.0/24.0 + result * x;
                return result;
            } else if (x < 4.0)
            {
                double result = -1.0/6.0;
                result = 5.0/2.0 + result * x;
                result = -55.0/4.0 + result * x;
                result = 65.0/2.0 + result * x;
                result = -655.0/24.0 + result * x;
                return result;
            } else
            {
                double result = 1.0/24.0;
                result = -5.0/6.0 + result * x;
                result = 25.0/4.0 + result * x;
                result = -125.0/6.0 + result * x;
                result = 625.0/24.0 + result * x;
                return result;
            }
            break;
        case 5:
            if ((x < 0.0) || (x >= 6.0))
            {
                return 0.0;
            } else if (x < 1.0)
            {
                double result = 1.0/120.0;
                result *= x;
                result *= x;
                result *= x;
                result *= x;
                result *= x;
                return result;
            } else if (x < 2.0)
            {
                double result = -1.0/24.0;
                result = 1.0/4.0 + result * x;
                result = -1.0/2.0 + result * x;
                result = 1.0/2.0 + result * x;
                result = -1.0/4.0 + result * x;
                result = 1.0/20.0 + result * x;
                return result;
            } else if (x < 3.0)
            {
                double result = 1.0/12.0;
                result = -1.0 + result * x;
                result = 9.0/2.0 + result * x;
                result = -19.0/2.0 + result * x;
                result = 39.0/4.0 + result * x;
                result = -79.0/20.0 + result * x;
                return result;
            } else if (x < 4.0)
            {
                double result = -1.0/12.0;
                result = 3.0/2.0 + result * x;
                result = -21.0/2.0 + result * x;
                result = 71.0/2.0 + result * x;
                result = -231.0/4.0 + result * x;
                result = 731.0/20.0 + result * x;
                return result;
            } else if (x < 5.0)
            {
                double result = 1.0/24.0;
                result = -1.0 + result * x;
                result = 19.0/2.0 + result * x;
                result = -89.0/2.0 + result * x;
                result = 409.0/4.0 + result * x;
                result = -1829.0/20.0 + result * x;
                return result;
            } else
            {
                double result = -1.0/120.0;
                result = 1.0/4.0 + result * x;
                result = -3.0 + result * x;
                result = 18.0 + result * x;
                result = -54.0 + result * x;
                result = 324.0/5.0 + result * x;
                return result;
            }
            break;
        case 6:
            if ((x < 0.0) || (x >= 7.0))
            {
                return 0.0;
            } else if (x < 1.0)
            {
                double result = 1.0/720.0;
                result *= x;
                result *= x;
                result *= x;
                result *= x;
                result *= x;
                result *= x;
                return result;
            } else if (x < 2.0)
            {
                double result = -1.0/120.0;
                result = 7.0/120.0 + result * x;
                result = -7.0/48.0 + result * x;
                result = 7.0/36.0 + result * x;
                result = -7.0/48.0 + result * x;
                result = 7.0/120.0 + result * x;
                result = -7.0/720.0 + result * x;
                return result;
            } else if (x < 3.0)
            {
                double result = 1.0/48.0;
                result = -7.0/24.0 + result * x;
                result = 77.0/48.0 + result * x;
                result = -161.0/36.0 + result * x;
                result = 329.0/48.0 + result * x;
                result = -133.0/24.0 + result * x;
                result = 1337.0/720.0 + result * x;
                return result;
            } else if (x < 4.0)
            {
                double result = -1.0/36.0;
                result = 7.0/12.0 + result * x;
                result = -119.0/24.0 + result * x;
                result = 196.0/9.0 + result * x;
                result = -1253.0/24.0 + result * x;
                result = 196.0/3.0 + result * x;
                result = -12089.0/360.0 + result * x;
                return result;
            } else if (x < 5.0)
            {
                double result = 1.0/48.0;
                result = -7.0/12.0 + result * x;
                result = 161.0/24.0 + result * x;
                result = -364.0/9.0 + result * x;
                result = 3227.0/24.0 + result * x;
                result = -700.0/3.0 + result * x;
                result = 59591.0/360.0 + result * x;
                return result;
            } else if (x < 6.0)
            {
                double result = -1.0/120.0;
                result = 7.0/24.0 + result * x;
                result = -203.0/48.0 + result * x;
                result = 1169.0/36.0 + result * x;
                result = -6671.0/48.0 + result * x;
                result = 7525.0/24.0 + result * x;
                result = -208943.0/720.0 + result * x;
                return result;
            } else
            {
                double result = 1.0/720.0;
                result = -7.0/120.0 + result * x;
                result = 49.0/48.0 + result * x;
                result = -343.0/36.0 + result * x;
                result = 2401.0/48.0 + result * x;
                result = -16807.0/120.0 + result * x;
                result = 117649.0/720.0 + result * x;
                return result;
            }
            break;
        case 7:
            if ((x < 0.0) || (x >= 8.0))
            {
                return 0.0;
            } else if (x < 1.0)
            {
                double result = 1.0/5040.0;
                result *= x;
                result *= x;
                result *= x;
                result *= x;
                result *= x;
                result *= x;
                result *= x;
                return result;
            } else if (x < 2.0)
            {
                double result = -1.0/720.0;
                result = 1.0/90.0 + result * x;
                result = -1.0/30.0 + result * x;
                result = 1.0/18.0 + result * x;
                result = -1.0/18.0 + result * x;
                result = 1.0/30.0 + result * x;
                result = -1.0/90.0 + result * x;
                result = 1.0/630.0 + result * x;
                return result;
            } else if (x < 3.0)
            {
                double result = 1.0/240.0;
                result = -1.0/15.0 + result * x;
                result = 13.0/30.0 + result * x;
                result = -3.0/2.0 + result * x;
                result = 55.0/18.0 + result * x;
                result = -37.0/10.0 + result * x;
                result = 223.0/90.0 + result * x;
                result = -149.0/210.0 + result * x;
                return result;
            } else if (x < 4.0)
            {
                double result = -1.0/144.0;
                result = 1.0/6.0 + result * x;
                result = -5.0/3.0 + result * x;
                result = 9.0 + result * x;
                result = -256.0/9.0 + result * x;
                result = 53.0 + result * x;
                result = -488.0/9.0 + result * x;
                result = 2477.0/105.0 + result * x;
                return result;
            } else if (x < 5.0)
            {
                double result = 1.0/144.0;
                result = -2.0/9.0 + result * x;
                result = 3.0 + result * x;
                result = -199.0/9.0 + result * x;
                result = 96.0 + result * x;
                result = -737.0/3.0 + result * x;
                result = 344.0 + result * x;
                result = -64249.0/315.0 + result * x;
                return result;
            } else if (x < 6.0)
            {
                double result = -1.0/240.0;
                result = 1.0/6.0 + result * x;
                result = -17.0/6.0 + result * x;
                result = 53.0/2.0 + result * x;
                result = -2647.0/18.0 + result * x;
                result = 967.0/2.0 + result * x;
                result = -15683.0/18.0 + result * x;
                result = 139459.0/210.0 + result * x;
                return result;
            } else if (x < 7.0)
            {
                double result = 1.0/720.0;
                result = -1.0/15.0 + result * x;
                result = 41.0/30.0 + result * x;
                result = -31.0/2.0 + result * x;
                result = 1889.0/18.0 + result * x;
                result = -4237.0/10.0 + result * x;
                result = 84881.0/90.0 + result * x;
                result = -187133.0/210.0 + result * x;
                return result;
            } else
            {
                double result = -1.0/5040.0;
                result = 1.0/90.0 + result * x;
                result = -4.0/15.0 + result * x;
                result = 32.0/9.0 + result * x;
                result = -256.0/9.0 + result * x;
                result = 2048.0/15.0 + result * x;
                result = -16384.0/45.0 + result * x;
                result = 131072.0/315.0 + result * x;
                return result;
            }
            break;
        /*case 8:
            if ((x < 0.0) || (x >= 9.0))
            {
                return 0.0;
            } else if (x < 1.0)
            {
                double result = 1.0/40320.0;
                result *= x;
                result *= x;
                result *= x;
                result *= x;
                result *= x;
                result *= x;
                result *= x;
                result *= x;
                return result;
            } else if (x < 2.0)
            {
                double result = -1.0/5040.0;
                result = 1.0/560.0 + result * x;
                result = -1.0/160.0 + result * x;
                result = 1.0/80.0 + result * x;
                result = -1.0/64.0 + result * x;
                result = 1.0/80.0 + result * x;
                result = -1.0/160.0 + result * x;
                result = 1.0/560.0 + result * x;
                result = -1.0/4480.0 + result * x;
                return result;
            } else if (x < 3.0)
            {
                double result = 1.0/1440.0;
                result = -1.0/80.0 + result * x;
                result = 3.0/32.0 + result * x;
                result = -31.0/80.0 + result * x;
                result = 63.0/64.0 + result * x;
                result = -127.0/80.0 + result * x;
                result = 51.0/32.0 + result * x;
                result = -73.0/80.0 + result * x;
                result = 1023.0/4480.0 + result * x;
                return result;
            } else if (x < 4.0)
            {
                double result = -1.0/720.0;
                result = 3.0/80.0 + result * x;
                result = -69.0/160.0 + result * x;
                result = 221.0/80.0 + result * x;
                result = -693.0/64.0 + result * x;
                result = 2141.0/80.0 + result * x;
                result = -6549.0/160.0 + result * x;
                result = 2843.0/80.0 + result * x;
                result = -60213.0/4480.0 + result * x;
                return result;
            } else if (x < 5.0)
            {
                double result = 1.0/576.0;
                result = -1.0/16.0 + result * x;
                result = 31.0/32.0 + result * x;
                result = -135.0/16.0 + result * x;
                result = 2891.0/64.0 + result * x;
                result = -2439.0/16.0 + result * x;
                result = 10159.0/32.0 + result * x;
                result = -5985.0/16.0 + result * x;
                result = 857291.0/4480.0 + result * x;
                return result;
            } else if (x < 6.0)
            {
                double result = -1.0/720.0;
                result = 1.0/16.0 + result * x;
                result = -39.0/32.0 + result * x;
                result = 215.0/16.0 + result * x;
                result = -5859.0/64.0 + result * x;
                result = 6311.0/16.0 + result * x;
                result = -33591.0/32.0 + result * x;
                result = 25265.0/16.0 + result * x;
                result = -4611459.0/4480.0 + result * x;
                return result;
            } else if (x < 7.0)
            {
                double result = 1.0/1440.0;
                result = -3.0/80.0 + result * x;
                result = 141.0/160.0 + result * x;
                result = -941.0/80.0 + result * x;
                result = 6237.0/64.0 + result * x;
                result = -41021.0/80.0 + result * x;
                result = 267501.0/160.0 + result * x;
                result = -246923.0/80.0 + result * x;
                result = 11064957.0/4480.0 + result * x;
                return result;
            } else if (x < 8.0)
            {
                double result = -1.0/5040.0;
                result = 1.0/80.0 + result * x;
                result = -11.0/32.0 + result * x;
                result = 431.0/80.0 + result * x;
                result = -3367.0/64.0 + result * x;
                result = 26207.0/80.0 + result * x;
                result = -40619.0/32.0 + result * x;
                result = 223673.0/80.0 + result * x;
                result = -11994247.0/4480.0 + result * x;
                return result;
            } else
            {
                double result = 1.0/40320.0;
                result = -1.0/560.0 + result * x;
                result = 9.0/160.0 + result * x;
                result = -81.0/80.0 + result * x;
                result = 729.0/64.0 + result * x;
                result = -6561.0/80.0 + result * x;
                result = 59049.0/160.0 + result * x;
                result = -531441.0/560.0 + result * x;
                result = 4782969.0/4480.0 + result * x;
                return result;
            }
            break;
        case 9:
            if ((x < 0.0) || (x >= 10.0))
            {
                return 0.0;
            } else if (x < 1.0)
            {
                double result = 1.0/362880.0;
                result *= x;
                result *= x;
                result *= x;
                result *= x;
                result *= x;
                result *= x;
                result *= x;
                result *= x;
                result *= x;
                return result;
            } else if (x < 2.0)
            {
                double result = -1.0/40320.0;
                result = 1.0/4032.0 + result * x;
                result = -1.0/1008.0 + result * x;
                result = 1.0/432.0 + result * x;
                result = -1.0/288.0 + result * x;
                result = 1.0/288.0 + result * x;
                result = -1.0/432.0 + result * x;
                result = 1.0/1008.0 + result * x;
                result = -1.0/4032.0 + result * x;
                result = 1.0/36288.0 + result * x;
                return result;
            } else if (x < 3.0)
            {
                double result = 1.0/10080.0;
                result = -1.0/504.0 + result * x;
                result = 17.0/1008.0 + result * x;
                result = -35.0/432.0 + result * x;
                result = 71.0/288.0 + result * x;
                result = -143.0/288.0 + result * x;
                result = 287.0/432.0 + result * x;
                result = -575.0/1008.0 + result * x;
                result = 1151.0/4032.0 + result * x;
                result = -329.0/5184.0 + result * x;
                return result;
            } else if (x < 4.0)
            {
                double result = -1.0/4320.0;
                result = 1.0/144.0 + result * x;
                result = -13.0/144.0 + result * x;
                result = 289.0/432.0 + result * x;
                result = -901.0/288.0 + result * x;
                result = 2773.0/288.0 + result * x;
                result = -8461.0/432.0 + result * x;
                result = 3667.0/144.0 + result * x;
                result = -11083.0/576.0 + result * x;
                result = 233893.0/36288.0 + result * x;
                return result;
            } else if (x < 5.0)
            {
                double result = 1.0/2880.0;
                result = -1.0/72.0 + result * x;
                result = 35.0/144.0 + result * x;
                result = -1055.0/432.0 + result * x;
                result = 4475.0/288.0 + result * x;
                result = -18731.0/288.0 + result * x;
                result = 77555.0/432.0 + result * x;
                result = -45485.0/144.0 + result * x;
                result = 185525.0/576.0 + result * x;
                result = -5271131.0/36288.0 + result * x;
                return result;
            } else if (x < 6.0)
            {
                double result = -1.0/2880.0;
                result = 5.0/288.0 + result * x;
                result = -55.0/144.0 + result * x;
                result = 2095.0/432.0 + result * x;
                result = -11275.0/288.0 + result * x;
                result = 60019.0/288.0 + result * x;
                result = -316195.0/432.0 + result * x;
                result = 235765.0/144.0 + result * x;
                result = -1220725.0/576.0 + result * x;
                result = 43947619.0/36288.0 + result * x;
                return result;
            } else if (x < 7.0)
            {
                double result = 1.0/4320.0;
                result = -1.0/72.0 + result * x;
                result = 53.0/144.0 + result * x;
                result = -2441.0/432.0 + result * x;
                result = 15941.0/288.0 + result * x;
                result = -103277.0/288.0 + result * x;
                result = 663581.0/432.0 + result * x;
                result = -604043.0/144.0 + result * x;
                result = 3818123.0/576.0 + result * x;
                result = -167683997.0/36288.0 + result * x;
                return result;
            } else if (x < 8.0)
            {
                double result = -1.0/10080.0;
                result = 1.0/144.0 + result * x;
                result = -31.0/144.0 + result * x;
                result = 1675.0/432.0 + result * x;
                result = -12871.0/288.0 + result * x;
                result = 98407.0/288.0 + result * x;
                result = -748207.0/432.0 + result * x;
                result = 807745.0/144.0 + result * x;
                result = -6064393.0/576.0 + result * x;
                result = 316559287.0/36288.0 + result * x;
                return result;
            } else if (x < 9.0)
            {
                double result = 1.0/40320.0;
                result = -1.0/504.0 + result * x;
                result = 71.0/1008.0 + result * x;
                result = -629.0/432.0 + result * x;
                result = 5561.0/288.0 + result * x;
                result = -49049.0/288.0 + result * x;
                result = 431441.0/432.0 + result * x;
                result = -3782969.0/1008.0 + result * x;
                result = 33046721.0/4032.0 + result * x;
                result = -287420489.0/36288.0 + result * x;
                return result;
            } else
            {
                double result = -1.0/362880.0;
                result = 1.0/4032.0 + result * x;
                result = -5.0/504.0 + result * x;
                result = 25.0/108.0 + result * x;
                result = -125.0/36.0 + result * x;
                result = 625.0/18.0 + result * x;
                result = -6250.0/27.0 + result * x;
                result = 62500.0/63.0 + result * x;
                result = -156250.0/63.0 + result * x;
                result = 1562500.0/567.0 + result * x;
                return result;
            }
            break;*/
        default:
            double p_double = static_cast<double>(p);
            
            if ((x < 0.0) || (x >= p_double+1.0))
            {
                return 0.0;
            } else
            {
                return (x / p_double) * uniformBSpline(x, p-1)
                        + ((p_double + 1.0 - x) / p_double) * uniformBSpline(x - 1.0, p-1);
            }
        }
    }
    
    inline double uniformBSplineDx(double x, size_t p) const
    {
        switch (p)
        {
        case 0:
            return 0.0;
            break;
        case 1:
            if ((x < 0.0) || (x >= 2.0))
            {
                return 0.0;
            } else if (x < 1.0)
            {
                return 1.0;
            } else
            {
                return -1.0;
            }
            break;
        case 2:
            if ((x < 0.0) || (x >= 3.0))
            {
                return 0.0;
            } else if (x < 1.0)
            {
                return x;
            } else if (x < 2.0)
            {
                return -2.0*x + 3.0;
            } else
            {
                return x - 3.0;
            }
            break;
        case 3:
            if ((x < 0.0) || (x >= 4.0))
            {
                return 0.0;
            } else if (x < 1.0)
            {
                return 1.0/2.0*x*x;
            } else if (x < 2.0)
            {
                return -3.0/2.0*x*x + 4.0*x - 2.0;
            } else if (x < 3.0)
            {
                return 3.0/2.0*x*x - 8.0*x + 10.0;
            } else
            {
                return -1.0/2.0*x*x + 4.0*x - 8.0;
            }
            break;
        case 4:
            if ((x < 0.0) || (x >= 5.0))
            {
                return 0.0;
            } else if (x < 1.0)
            {
                return 1.0/6.0*x*x*x;
            } else if (x < 2.0)
            {
                return -2.0/3.0*x*x*x + 5.0/2.0*x*x - 5.0/2.0*x + 5.0/6.0;
            } else if (x < 3.0)
            {
                return x*x*x - 15.0/2.0*x*x + 35.0/2.0*x - 25.0/2.0;
            } else if (x < 4.0)
            {
                return -2.0/3.0*x*x*x + 15.0/2.0*x*x - 55.0/2.0*x + 65.0/2.0;
            } else
            {
                return 1.0/6.0*x*x*x - 5.0/2.0*x*x + 25.0/2.0*x - 125.0/6.0;
            }
            break;
        case 5:
            if ((x < 0.0) || (x >= 6.0))
            {
                return 0.0;
            } else if (x < 1.0)
            {
                double result = 1.0/24.0;
                result *= x;
                result *= x;
                result *= x;
                result *= x;
                return result;
            } else if (x < 2.0)
            {
                double result = -5.0/24.0;
                result = 1.0 + result * x;
                result = -3.0/2.0 + result * x;
                result = 1.0 + result * x;
                result = -1.0/4.0 + result * x;
                return result;
            } else if (x < 3.0)
            {
                double result = 5.0/12.0;
                result = -4.0 + result * x;
                result = 27.0/2.0 + result * x;
                result = -19.0 + result * x;
                result = 39.0/4.0 + result * x;
                return result;
            } else if (x < 4.0)
            {
                double result = -5.0/12.0;
                result = 6.0 + result * x;
                result = -63.0/2.0 + result * x;
                result = 71.0 + result * x;
                result = -231.0/4.0 + result * x;
                return result;
            } else if (x < 5.0)
            {
                double result = 5.0/24.0;
                result = -4.0 + result * x;
                result = 57.0/2.0 + result * x;
                result = -89.0 + result * x;
                result = 409.0/4.0 + result * x;
                return result;
            } else
            {
                double result = -1.0/24.0;
                result = 1.0 + result * x;
                result = -9.0 + result * x;
                result = 36.0 + result * x;
                result = -54.0 + result * x;
                return result;
            }
            break;
        case 6:
            if ((x < 0.0) || (x >= 7.0))
            {
                return 0.0;
            } else if (x < 1.0)
            {
                double result = 1.0/120.0;
                result *= x;
                result *= x;
                result *= x;
                result *= x;
                result *= x;
                return result;
            } else if (x < 2.0)
            {
                double result = -1.0/20.0;
                result = 7.0/24.0 + result * x;
                result = -7.0/12.0 + result * x;
                result = 7.0/12.0 + result * x;
                result = -7.0/24.0 + result * x;
                result = 7.0/120.0 + result * x;
                return result;
            } else if (x < 3.0)
            {
                double result = 1.0/8.0;
                result = -35.0/24.0 + result * x;
                result = 77.0/12.0 + result * x;
                result = -161.0/12.0 + result * x;
                result = 329.0/24.0 + result * x;
                result = -133.0/24.0 + result * x;
                return result;
            } else if (x < 4.0)
            {
                double result = -1.0/6.0;
                result = 35.0/12.0 + result * x;
                result = -119.0/6.0 + result * x;
                result = 196.0/3.0 + result * x;
                result = -1253.0/12.0 + result * x;
                result = 196.0/3.0 + result * x;
                return result;
            } else if (x < 5.0)
            {
                double result = 1.0/8.0;
                result = -35.0/12.0 + result * x;
                result = 161.0/6.0 + result * x;
                result = -364.0/3.0 + result * x;
                result = 3227.0/12.0 + result * x;
                result = -700.0/3.0 + result * x;
                return result;
            } else if (x < 6.0)
            {
                double result = -1.0/20.0;
                result = 35.0/24.0 + result * x;
                result = -203.0/12.0 + result * x;
                result = 1169.0/12.0 + result * x;
                result = -6671.0/24.0 + result * x;
                result = 7525.0/24.0 + result * x;
                return result;
            } else
            {
                double result = 1.0/120.0;
                result = -7.0/24.0 + result * x;
                result = 49.0/12.0 + result * x;
                result = -343.0/12.0 + result * x;
                result = 2401.0/24.0 + result * x;
                result = -16807.0/120.0 + result * x;
                return result;
            }
            break;
        case 7:
            if ((x < 0.0) || (x >= 8.0))
            {
                return 0.0;
            } else if (x < 1.0)
            {
                double result = 1.0/720.0;
                result *= x;
                result *= x;
                result *= x;
                result *= x;
                result *= x;
                result *= x;
                return result;
            } else if (x < 2.0)
            {
                double result = -7.0/720.0;
                result = 1.0/15.0 + result * x;
                result = -1.0/6.0 + result * x;
                result = 2.0/9.0 + result * x;
                result = -1.0/6.0 + result * x;
                result = 1.0/15.0 + result * x;
                result = -1.0/90.0 + result * x;
                return result;
            } else if (x < 3.0)
            {
                double result = 7.0/240.0;
                result = -2.0/5.0 + result * x;
                result = 13.0/6.0 + result * x;
                result = -6.0 + result * x;
                result = 55.0/6.0 + result * x;
                result = -37.0/5.0 + result * x;
                result = 223.0/90.0 + result * x;
                return result;
            } else if (x < 4.0)
            {
                double result = -7.0/144.0;
                result = 1.0 + result * x;
                result = -25.0/3.0 + result * x;
                result = 36.0 + result * x;
                result = -256.0/3.0 + result * x;
                result = 106.0 + result * x;
                result = -488.0/9.0 + result * x;
                return result;
            } else if (x < 5.0)
            {
                double result = 7.0/144.0;
                result = -4.0/3.0 + result * x;
                result = 15.0 + result * x;
                result = -796.0/9.0 + result * x;
                result = 288.0 + result * x;
                result = -1474.0/3.0 + result * x;
                result = 344.0 + result * x;
                return result;
            } else if (x < 6.0)
            {
                double result = -7.0/240.0;
                result = 1.0 + result * x;
                result = -85.0/6.0 + result * x;
                result = 106.0 + result * x;
                result = -2647.0/6.0 + result * x;
                result = 967.0 + result * x;
                result = -15683.0/18.0 + result * x;
                return result;
            } else if (x < 7.0)
            {
                double result = 7.0/720.0;
                result = -2.0/5.0 + result * x;
                result = 41.0/6.0 + result * x;
                result = -62.0 + result * x;
                result = 1889.0/6.0 + result * x;
                result = -4237.0/5.0 + result * x;
                result = 84881.0/90.0 + result * x;
                return result;
            } else
            {
                double result = -1.0/720.0;
                result = 1.0/15.0 + result * x;
                result = -4.0/3.0 + result * x;
                result = 128.0/9.0 + result * x;
                result = -256.0/3.0 + result * x;
                result = 4096.0/15.0 + result * x;
                result = -16384.0/45.0 + result * x;
                return result;
            }
            break;
        /*case 8:
            if ((x < 0.0) || (x >= 9.0))
            {
                return 0.0;
            } else if (x < 1.0)
            {
                double result = 1.0/5040.0;
                result *= x;
                result *= x;
                result *= x;
                result *= x;
                result *= x;
                result *= x;
                result *= x;
                return result;
            } else if (x < 2.0)
            {
                double result = -1.0/630.0;
                result = 1.0/80.0 + result * x;
                result = -3.0/80.0 + result * x;
                result = 1.0/16.0 + result * x;
                result = -1.0/16.0 + result * x;
                result = 3.0/80.0 + result * x;
                result = -1.0/80.0 + result * x;
                result = 1.0/560.0 + result * x;
                return result;
            } else if (x < 3.0)
            {
                double result = 1.0/180.0;
                result = -7.0/80.0 + result * x;
                result = 9.0/16.0 + result * x;
                result = -31.0/16.0 + result * x;
                result = 63.0/16.0 + result * x;
                result = -381.0/80.0 + result * x;
                result = 51.0/16.0 + result * x;
                result = -73.0/80.0 + result * x;
                return result;
            } else if (x < 4.0)
            {
                double result = -1.0/90.0;
                result = 21.0/80.0 + result * x;
                result = -207.0/80.0 + result * x;
                result = 221.0/16.0 + result * x;
                result = -693.0/16.0 + result * x;
                result = 6423.0/80.0 + result * x;
                result = -6549.0/80.0 + result * x;
                result = 2843.0/80.0 + result * x;
                return result;
            } else if (x < 5.0)
            {
                double result = 1.0/72.0;
                result = -7.0/16.0 + result * x;
                result = 93.0/16.0 + result * x;
                result = -675.0/16.0 + result * x;
                result = 2891.0/16.0 + result * x;
                result = -7317.0/16.0 + result * x;
                result = 10159.0/16.0 + result * x;
                result = -5985.0/16.0 + result * x;
                return result;
            } else if (x < 6.0)
            {
                double result = -1.0/90.0;
                result = 7.0/16.0 + result * x;
                result = -117.0/16.0 + result * x;
                result = 1075.0/16.0 + result * x;
                result = -5859.0/16.0 + result * x;
                result = 18933.0/16.0 + result * x;
                result = -33591.0/16.0 + result * x;
                result = 25265.0/16.0 + result * x;
                return result;
            } else if (x < 7.0)
            {
                double result = 1.0/180.0;
                result = -21.0/80.0 + result * x;
                result = 423.0/80.0 + result * x;
                result = -941.0/16.0 + result * x;
                result = 6237.0/16.0 + result * x;
                result = -123063.0/80.0 + result * x;
                result = 267501.0/80.0 + result * x;
                result = -246923.0/80.0 + result * x;
                return result;
            } else if (x < 8.0)
            {
                double result = -1.0/630.0;
                result = 7.0/80.0 + result * x;
                result = -33.0/16.0 + result * x;
                result = 431.0/16.0 + result * x;
                result = -3367.0/16.0 + result * x;
                result = 78621.0/80.0 + result * x;
                result = -40619.0/16.0 + result * x;
                result = 223673.0/80.0 + result * x;
                return result;
            } else
            {
                double result = 1.0/5040.0;
                result = -1.0/80.0 + result * x;
                result = 27.0/80.0 + result * x;
                result = -81.0/16.0 + result * x;
                result = 729.0/16.0 + result * x;
                result = -19683.0/80.0 + result * x;
                result = 59049.0/80.0 + result * x;
                result = -531441.0/560.0 + result * x;
                return result;
            }
            break;
        case 9:
            if ((x < 0.0) || (x >= 10.0))
            {
                return 0.0;
            } else if (x < 1.0)
            {
                double result = 1.0/40320.0;
                result *= x;
                result *= x;
                result *= x;
                result *= x;
                result *= x;
                result *= x;
                result *= x;
                result *= x;
                return result;
            } else if (x < 2.0)
            {
                double result = -1.0/4480.0;
                result = 1.0/504.0 + result * x;
                result = -1.0/144.0 + result * x;
                result = 1.0/72.0 + result * x;
                result = -5.0/288.0 + result * x;
                result = 1.0/72.0 + result * x;
                result = -1.0/144.0 + result * x;
                result = 1.0/504.0 + result * x;
                result = -1.0/4032.0 + result * x;
                return result;
            } else if (x < 3.0)
            {
                double result = 1.0/1120.0;
                result = -1.0/63.0 + result * x;
                result = 17.0/144.0 + result * x;
                result = -35.0/72.0 + result * x;
                result = 355.0/288.0 + result * x;
                result = -143.0/72.0 + result * x;
                result = 287.0/144.0 + result * x;
                result = -575.0/504.0 + result * x;
                result = 1151.0/4032.0 + result * x;
                return result;
            } else if (x < 4.0)
            {
                double result = -1.0/480.0;
                result = 1.0/18.0 + result * x;
                result = -91.0/144.0 + result * x;
                result = 289.0/72.0 + result * x;
                result = -4505.0/288.0 + result * x;
                result = 2773.0/72.0 + result * x;
                result = -8461.0/144.0 + result * x;
                result = 3667.0/72.0 + result * x;
                result = -11083.0/576.0 + result * x;
                return result;
            } else if (x < 5.0)
            {
                double result = 1.0/320.0;
                result = -1.0/9.0 + result * x;
                result = 245.0/144.0 + result * x;
                result = -1055.0/72.0 + result * x;
                result = 22375.0/288.0 + result * x;
                result = -18731.0/72.0 + result * x;
                result = 77555.0/144.0 + result * x;
                result = -45485.0/72.0 + result * x;
                result = 185525.0/576.0 + result * x;
                return result;
            } else if (x < 6.0)
            {
                double result = -1.0/320.0;
                result = 5.0/36.0 + result * x;
                result = -385.0/144.0 + result * x;
                result = 2095.0/72.0 + result * x;
                result = -56375.0/288.0 + result * x;
                result = 60019.0/72.0 + result * x;
                result = -316195.0/144.0 + result * x;
                result = 235765.0/72.0 + result * x;
                result = -1220725.0/576.0 + result * x;
                return result;
            } else if (x < 7.0)
            {
                double result = 1.0/480.0;
                result = -1.0/9.0 + result * x;
                result = 371.0/144.0 + result * x;
                result = -2441.0/72.0 + result * x;
                result = 79705.0/288.0 + result * x;
                result = -103277.0/72.0 + result * x;
                result = 663581.0/144.0 + result * x;
                result = -604043.0/72.0 + result * x;
                result = 3818123.0/576.0 + result * x;
                return result;
            } else if (x < 8.0)
            {
                double result = -1.0/1120.0;
                result = 1.0/18.0 + result * x;
                result = -217.0/144.0 + result * x;
                result = 1675.0/72.0 + result * x;
                result = -64355.0/288.0 + result * x;
                result = 98407.0/72.0 + result * x;
                result = -748207.0/144.0 + result * x;
                result = 807745.0/72.0 + result * x;
                result = -6064393.0/576.0 + result * x;
                return result;
            } else if (x < 9.0)
            {
                double result = 1.0/4480.0;
                result = -1.0/63.0 + result * x;
                result = 71.0/144.0 + result * x;
                result = -629.0/72.0 + result * x;
                result = 27805.0/288.0 + result * x;
                result = -49049.0/72.0 + result * x;
                result = 431441.0/144.0 + result * x;
                result = -3782969.0/504.0 + result * x;
                result = 33046721.0/4032.0 + result * x;
                return result;
            } else
            {
                double result = -1.0/40320.0;
                result = 1.0/504.0 + result * x;
                result = -5.0/72.0 + result * x;
                result = 25.0/18.0 + result * x;
                result = -625.0/36.0 + result * x;
                result = 1250.0/9.0 + result * x;
                result = -6250.0/9.0 + result * x;
                result = 125000.0/63.0 + result * x;
                result = -156250.0/63.0 + result * x;
                return result;
            }
            break;*/
        default:
            double p_double = static_cast<double>(p);
            
            if ((x < 0.0) || (x >= p_double+1.0))
            {
                return 0.0;
            } else
            {
                return uniformBSpline(x, p-1) - uniformBSpline(x - 1.0, p-1);
            }
        }
    }
    
    inline double uniformBSplineDxDx(double x, size_t p) const
    {
        switch (p)
        {
        case 0:
        case 1:
            return 0.0;
            break;
        case 2:
            if ((x < 0.0) || (x >= 3.0))
            {
                return 0.0;
            } else if (x < 1.0)
            {
                return 1.0;
            } else if (x < 2.0)
            {
                return -2.0;
            } else
            {
                return 1.0;
            }
            break;
        case 3:
            if ((x < 0.0) || (x >= 4.0))
            {
                return 0.0;
            } else if (x < 1.0)
            {
                return x;
            } else if (x < 2.0)
            {
                return -3.0*x + 4.0;
            } else if (x < 3.0)
            {
                return 3.0*x - 8.0;
            } else
            {
                return -x + 4.0;
            }
            break;
        case 4:
            if ((x < 0.0) || (x >= 5.0))
            {
                return 0.0;
            } else if (x < 1.0)
            {
                return 1.0/2.0*x*x;
            } else if (x < 2.0)
            {
                return -2.0*x*x + 5.0*x - 5.0/2.0;
            } else if (x < 3.0)
            {
                return 3.0*x*x - 15.0*x + 35.0/2.0;
            } else if (x < 4.0)
            {
                return -2.0*x*x + 15.0*x - 55.0/2.0;
            } else
            {
                return 1.0/2.0*x*x - 5.0*x + 25.0/2.0;
            }
            break;
        case 5:
            if ((x < 0.0) || (x >= 6.0))
            {
                return 0.0;
            } else if (x < 1.0)
            {
                return 1.0/6.0*x*x*x;
            } else if (x < 2.0)
            {
                return -5.0/6.0*x*x*x + 3.0*x*x - 3.0*x + 1.0;
            } else if (x < 3.0)
            {
                return 5.0/3.0*x*x*x - 12.0*x*x + 27.0*x - 19.0;
            } else if (x < 4.0)
            {
                return -5.0/3.0*x*x*x + 18.0*x*x - 63.0*x + 71.0;
            } else if (x < 5.0)
            {
                return 5.0/6.0*x*x*x - 12.0*x*x + 57.0*x - 89.0;
            } else
            {
                return -1.0/6.0*x*x*x + 3.0*x*x - 18.0*x + 36.0;
            }
            break;
        case 6:
            if ((x < 0.0) || (x >= 7.0))
            {
                return 0.0;
            } else if (x < 1.0)
            {
                double result = 1.0/24.0;
                result *= x;
                result *= x;
                result *= x;
                result *= x;
                return result;
            } else if (x < 2.0)
            {
                double result = -1.0/4.0;
                result = 7.0/6.0 + result * x;
                result = -7.0/4.0 + result * x;
                result = 7.0/6.0 + result * x;
                result = -7.0/24.0 + result * x;
                return result;
            } else if (x < 3.0)
            {
                double result = 5.0/8.0;
                result = -35.0/6.0 + result * x;
                result = 77.0/4.0 + result * x;
                result = -161.0/6.0 + result * x;
                result = 329.0/24.0 + result * x;
                return result;
            } else if (x < 4.0)
            {
                double result = -5.0/6.0;
                result = 35.0/3.0 + result * x;
                result = -119.0/2.0 + result * x;
                result = 392.0/3.0 + result * x;
                result = -1253.0/12.0 + result * x;
                return result;
            } else if (x < 5.0)
            {
                double result = 5.0/8.0;
                result = -35.0/3.0 + result * x;
                result = 161.0/2.0 + result * x;
                result = -728.0/3.0 + result * x;
                result = 3227.0/12.0 + result * x;
                return result;
            } else if (x < 6.0)
            {
                double result = -1.0/4.0;
                result = 35.0/6.0 + result * x;
                result = -203.0/4.0 + result * x;
                result = 1169.0/6.0 + result * x;
                result = -6671.0/24.0 + result * x;
                return result;
            } else
            {
                double result = 1.0/24.0;
                result = -7.0/6.0 + result * x;
                result = 49.0/4.0 + result * x;
                result = -343.0/6.0 + result * x;
                result = 2401.0/24.0 + result * x;
                return result;
            }
            break;
        case 7:
            if ((x < 0.0) || (x >= 8.0))
            {
                return 0.0;
            } else if (x < 1.0)
            {
                double result = 1.0/120.0;
                result *= x;
                result *= x;
                result *= x;
                result *= x;
                result *= x;
                return result;
            } else if (x < 2.0)
            {
                double result = -7.0/120.0;
                result = 1.0/3.0 + result * x;
                result = -2.0/3.0 + result * x;
                result = 2.0/3.0 + result * x;
                result = -1.0/3.0 + result * x;
                result = 1.0/15.0 + result * x;
                return result;
            } else if (x < 3.0)
            {
                double result = 7.0/40.0;
                result = -2.0 + result * x;
                result = 26.0/3.0 + result * x;
                result = -18.0 + result * x;
                result = 55.0/3.0 + result * x;
                result = -37.0/5.0 + result * x;
                return result;
            } else if (x < 4.0)
            {
                double result = -7.0/24.0;
                result = 5.0 + result * x;
                result = -100.0/3.0 + result * x;
                result = 108.0 + result * x;
                result = -512.0/3.0 + result * x;
                result = 106.0 + result * x;
                return result;
            } else if (x < 5.0)
            {
                double result = 7.0/24.0;
                result = -20.0/3.0 + result * x;
                result = 60.0 + result * x;
                result = -796.0/3.0 + result * x;
                result = 576.0 + result * x;
                result = -1474.0/3.0 + result * x;
                return result;
            } else if (x < 6.0)
            {
                double result = -7.0/40.0;
                result = 5.0 + result * x;
                result = -170.0/3.0 + result * x;
                result = 318.0 + result * x;
                result = -2647.0/3.0 + result * x;
                result = 967.0 + result * x;
                return result;
            } else if (x < 7.0)
            {
                double result = 7.0/120.0;
                result = -2.0 + result * x;
                result = 82.0/3.0 + result * x;
                result = -186.0 + result * x;
                result = 1889.0/3.0 + result * x;
                result = -4237.0/5.0 + result * x;
                return result;
            } else
            {
                double result = -1.0/120.0;
                result = 1.0/3.0 + result * x;
                result = -16.0/3.0 + result * x;
                result = 128.0/3.0 + result * x;
                result = -512.0/3.0 + result * x;
                result = 4096.0/15.0 + result * x;
                return result;
            }
            break;
        /*case 8:
            if ((x < 0.0) || (x >= 9.0))
            {
                return 0.0;
            } else if (x < 1.0)
            {
                double result = 1.0/720.0;
                result *= x;
                result *= x;
                result *= x;
                result *= x;
                result *= x;
                result *= x;
                return result;
            } else if (x < 2.0)
            {
                double result = -1.0/90.0;
                result = 3.0/40.0 + result * x;
                result = -3.0/16.0 + result * x;
                result = 1.0/4.0 + result * x;
                result = -3.0/16.0 + result * x;
                result = 3.0/40.0 + result * x;
                result = -1.0/80.0 + result * x;
                return result;
            } else if (x < 3.0)
            {
                double result = 7.0/180.0;
                result = -21.0/40.0 + result * x;
                result = 45.0/16.0 + result * x;
                result = -31.0/4.0 + result * x;
                result = 189.0/16.0 + result * x;
                result = -381.0/40.0 + result * x;
                result = 51.0/16.0 + result * x;
                return result;
            } else if (x < 4.0)
            {
                double result = -7.0/90.0;
                result = 63.0/40.0 + result * x;
                result = -207.0/16.0 + result * x;
                result = 221.0/4.0 + result * x;
                result = -2079.0/16.0 + result * x;
                result = 6423.0/40.0 + result * x;
                result = -6549.0/80.0 + result * x;
                return result;
            } else if (x < 5.0)
            {
                double result = 7.0/72.0;
                result = -21.0/8.0 + result * x;
                result = 465.0/16.0 + result * x;
                result = -675.0/4.0 + result * x;
                result = 8673.0/16.0 + result * x;
                result = -7317.0/8.0 + result * x;
                result = 10159.0/16.0 + result * x;
                return result;
            } else if (x < 6.0)
            {
                double result = -7.0/90.0;
                result = 21.0/8.0 + result * x;
                result = -585.0/16.0 + result * x;
                result = 1075.0/4.0 + result * x;
                result = -17577.0/16.0 + result * x;
                result = 18933.0/8.0 + result * x;
                result = -33591.0/16.0 + result * x;
                return result;
            } else if (x < 7.0)
            {
                double result = 7.0/180.0;
                result = -63.0/40.0 + result * x;
                result = 423.0/16.0 + result * x;
                result = -941.0/4.0 + result * x;
                result = 18711.0/16.0 + result * x;
                result = -123063.0/40.0 + result * x;
                result = 267501.0/80.0 + result * x;
                return result;
            } else if (x < 8.0)
            {
                double result = -1.0/90.0;
                result = 21.0/40.0 + result * x;
                result = -165.0/16.0 + result * x;
                result = 431.0/4.0 + result * x;
                result = -10101.0/16.0 + result * x;
                result = 78621.0/40.0 + result * x;
                result = -40619.0/16.0 + result * x;
                return result;
            } else
            {
                double result = 1.0/720.0;
                result = -3.0/40.0 + result * x;
                result = 27.0/16.0 + result * x;
                result = -81.0/4.0 + result * x;
                result = 2187.0/16.0 + result * x;
                result = -19683.0/40.0 + result * x;
                result = 59049.0/80.0 + result * x;
                return result;
            }
            break;
        case 9:
            if ((x < 0.0) || (x >= 10.0))
            {
                return 0.0;
            } else if (x < 1.0)
            {
                double result = 1.0/5040.0;
                result *= x;
                result *= x;
                result *= x;
                result *= x;
                result *= x;
                result *= x;
                result *= x;
                return result;
            } else if (x < 2.0)
            {
                double result = -1.0/560.0;
                result = 1.0/72.0 + result * x;
                result = -1.0/24.0 + result * x;
                result = 5.0/72.0 + result * x;
                result = -5.0/72.0 + result * x;
                result = 1.0/24.0 + result * x;
                result = -1.0/72.0 + result * x;
                result = 1.0/504.0 + result * x;
                return result;
            } else if (x < 3.0)
            {
                double result = 1.0/140.0;
                result = -1.0/9.0 + result * x;
                result = 17.0/24.0 + result * x;
                result = -175.0/72.0 + result * x;
                result = 355.0/72.0 + result * x;
                result = -143.0/24.0 + result * x;
                result = 287.0/72.0 + result * x;
                result = -575.0/504.0 + result * x;
                return result;
            } else if (x < 4.0)
            {
                double result = -1.0/60.0;
                result = 7.0/18.0 + result * x;
                result = -91.0/24.0 + result * x;
                result = 1445.0/72.0 + result * x;
                result = -4505.0/72.0 + result * x;
                result = 2773.0/24.0 + result * x;
                result = -8461.0/72.0 + result * x;
                result = 3667.0/72.0 + result * x;
                return result;
            } else if (x < 5.0)
            {
                double result = 1.0/40.0;
                result = -7.0/9.0 + result * x;
                result = 245.0/24.0 + result * x;
                result = -5275.0/72.0 + result * x;
                result = 22375.0/72.0 + result * x;
                result = -18731.0/24.0 + result * x;
                result = 77555.0/72.0 + result * x;
                result = -45485.0/72.0 + result * x;
                return result;
            } else if (x < 6.0)
            {
                double result = -1.0/40.0;
                result = 35.0/36.0 + result * x;
                result = -385.0/24.0 + result * x;
                result = 10475.0/72.0 + result * x;
                result = -56375.0/72.0 + result * x;
                result = 60019.0/24.0 + result * x;
                result = -316195.0/72.0 + result * x;
                result = 235765.0/72.0 + result * x;
                return result;
            } else if (x < 7.0)
            {
                double result = 1.0/60.0;
                result = -7.0/9.0 + result * x;
                result = 371.0/24.0 + result * x;
                result = -12205.0/72.0 + result * x;
                result = 79705.0/72.0 + result * x;
                result = -103277.0/24.0 + result * x;
                result = 663581.0/72.0 + result * x;
                result = -604043.0/72.0 + result * x;
                return result;
            } else if (x < 8.0)
            {
                double result = -1.0/140.0;
                result = 7.0/18.0 + result * x;
                result = -217.0/24.0 + result * x;
                result = 8375.0/72.0 + result * x;
                result = -64355.0/72.0 + result * x;
                result = 98407.0/24.0 + result * x;
                result = -748207.0/72.0 + result * x;
                result = 807745.0/72.0 + result * x;
                return result;
            } else if (x < 9.0)
            {
                double result = 1.0/560.0;
                result = -1.0/9.0 + result * x;
                result = 71.0/24.0 + result * x;
                result = -3145.0/72.0 + result * x;
                result = 27805.0/72.0 + result * x;
                result = -49049.0/24.0 + result * x;
                result = 431441.0/72.0 + result * x;
                result = -3782969.0/504.0 + result * x;
                return result;
            } else
            {
                double result = -1.0/5040.0;
                result = 1.0/72.0 + result * x;
                result = -5.0/12.0 + result * x;
                result = 125.0/18.0 + result * x;
                result = -625.0/9.0 + result * x;
                result = 1250.0/3.0 + result * x;
                result = -12500.0/9.0 + result * x;
                result = 125000.0/63.0 + result * x;
                return result;
            }
            break;*/
        default:
            double p_double = static_cast<double>(p);
            
            if ((x < 0.0) || (x >= p_double+1.0))
            {
                return 0.0;
            } else
            {
                return uniformBSpline(x, p-2) - 2.0 * uniformBSpline(x - 1.0, p-2)
                                              + uniformBSpline(x - 2.0, p-2);
            }
        }
    }

    /*inline double UniformBSpline(double x, size_t p) const
    {
        if (p == 0)
        {
            return ((0.0 <= x) && (x < 1.0)) ? 1.0 : 0.0;
        } else if (p == 1)
        {
            if ((x < 0.0) || (x >= 2.0))
            {
                return 0.0;
            } else if (x < 1.0)
            {
                return x;
            } else
            {
                return 2.0 - x;
            }
        } else if (p == 2)
        {
            if ((x < 0.0) || (x >= 3.0))
            {
                return 0.0;
            } else if (x < 1.0)
            {
                return x*x / 2.0;
            } else if (x < 2.0)
            {
                return (-2.0*x*x + 6.0*x - 3.0) / 2.0;
            } else
            {
                return (x-3.0)*(x-3.0) / 2.0;
            }
        } else if (p == 3)
        {
            if ((x < 0.0) || (x >= 4.0))
            {
                return 0.0;
            } else if (x < 1.0)
            {
                return x*x*x / 6.0;
            } else if (x < 2.0)
            {
                return 2.0/3.0 - (x-2.0)*(x-2.0) - (x-2.0)*(x-2.0)*(x-2.0) / 2.0;
            } else if (x < 3.0)
            {
                return 2.0/3.0 - (x-2.0)*(x-2.0) + (x-2.0)*(x-2.0)*(x-2.0) / 2.0;
            } else
            {
                return -(x-4.0)*(x-4.0)*(x-4.0) / 6.0;
            }
        } else
        {
            double p_double = static_cast<double>(p);
            
            if ((x < 0.0) || (x >= p_double+1.0))
            {
                return 0.0;
            } else
            {
                return (x / p_double) * UniformBSpline(x, p - 1)
                       + ((p_double + 1.0 - x) / p_double) * UniformBSpline(x - 1.0, p - 1);
            }
        }
    }
    
    inline double UniformBSplineDx(double x, size_t p) const
    {
        if (p == 0)
        {
            return 0.0;
        } else if (p == 1)
        {
            if ((x < 0.0) || (x >= 2.0))
            {
                return 0.0;
            } else if (x < 1.0)
            {
                return 1.0;
            } else
            {
                return -1.0;
            }
        } else if (p == 2)
        {
            if ((x < 0.0) || (x >= 3.0))
            {
                return 0.0;
            } else if (x < 1.0)
            {
                return x;
            } else if (x < 2.0)
            {
                return 3.0 - 2.0 * x;
            } else
            {
                return x - 3.0;
            }
        } else if (p == 3)
        {
            if ((x < 0.0) || (x >= 4.0))
            {
                return 0.0;
            } else if (x < 1.0)
            {
                return x*x / 2.0;
            } else if (x < 2.0)
            {
                return -2.0 * (x-2.0) - (x-2.0)*(x-2.0) * 3.0 / 2.0;
            } else if (x < 3.0)
            {
                return -2.0 * (x-2.0) + (x-2.0)*(x-2.0) * 3.0 / 2.0;
            } else
            {
                return -(x-4.0)*(x-4.0) / 2.0;
            }
        } else
        {
            double p_double = static_cast<double>(p);
            
            if ((x < 0.0) || (x >= p_double+1.0))
            {
                return 0.0;
            } else
            {
                return UniformBSpline(x, p - 1) - UniformBSpline(x - 1.0, p - 1);
            }
        }
    }
    
    inline double UniformBSplineDxDx(double x, size_t p) const
    {
        if (p <= 1)
        {
            return 0.0;
        } else if (p == 2)
        {
            if ((x < 0.0) || (x >= 3.0))
            {
                return 0.0;
            } else if (x < 1.0)
            {
                return 1.0;
            } else if (x < 2.0)
            {
                return 2.0;
            } else
            {
                return 1.0;
            }
        } else if (p == 3)
        {
            if ((x < 0.0) || (x >= 4.0))
            {
                return 0.0;
            } else if (x < 1.0)
            {
                return x;
            } else if (x < 2.0)
            {
                return 4.0 - 3.0 * x;
            } else if (x < 3.0)
            {
                return 3.0 * x - 8.0;
            } else
            {
                return -(x-4.0);
            }
        } else
        {
            double p_double = static_cast<double>(p);
            
            if ((x < 0.0) || (x >= p_double+1.0))
            {
                return 0.0;
            } else
            {
                return UniformBSpline(x, p - 2) - 2.0 * UniformBSpline(x - 1.0, p - 2)
                                                + UniformBSpline(x - 2.0, p - 2);
            }
        }
    }*/
    
    inline double eval(LT level, IT index, double x) const
    {
        double hinv = static_cast<double>(1 << level);
        
        return uniformBSpline(
                x * hinv + static_cast<double>(this->degree + 1) / 2 - index,
                this->degree);
    }
    
    inline double evalDx(LT level, IT index, double x) const
    {
        double hinv = static_cast<double>(1 << level);
        
        return hinv * uniformBSplineDx(
                x * hinv + static_cast<double>(this->degree + 1) / 2 - index,
                this->degree);
    }
    
    inline double evalDxDx(LT level, IT index, double x) const
    {
        double hinv = static_cast<double>(1 << level);
        
        return hinv*hinv * uniformBSplineDxDx(
                x * hinv + static_cast<double>(this->degree + 1) / 2 - index,
                this->degree);
    }
    
    inline size_t getDegree() const { return degree; }
};

// default type-def (unsigned int for level and index)
typedef BsplineBasis<unsigned int, unsigned int> SBsplineBase;

}
}

#endif /* BSPLINE_BASE_HPP */
