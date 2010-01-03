/*****************************************************************************/
/* This file is part of sgpp, a program package making use of spatially      */
/* adaptive sparse grids to solve numerical problems                         */
/*                                                                           */
/* Copyright (C) 2008 Jörg Blank (blankj@in.tum.de)                          */
/* Copyright (C) 2009-2010 Alexander Heinecke (Alexander.Heinecke@mytum.de)  */
/* Copyright (C) 2009-2010 Dirk Pflueger (pflueged@in.tum.de)                */
/*                                                                           */
/* sgpp is free software; you can redistribute it and/or modify              */
/* it under the terms of the GNU Lesser General Public License as published  */
/* by the Free Software Foundation; either version 3 of the License, or      */
/* (at your option) any later version.                                       */
/*                                                                           */
/* sgpp is distributed in the hope that it will be useful,                   */
/* but WITHOUT ANY WARRANTY; without even the implied warranty of            */
/* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the             */
/* GNU Lesser General Public License for more details.                       */
/*                                                                           */
/* You should have received a copy of the GNU Lesser General Public License  */
/* along with sgpp; if not, write to the Free Software                       */
/* Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA */
/* or see <http://www.gnu.org/licenses/>.                                    */
/*****************************************************************************/

#ifndef MODIFIED_WAVELET_BASE_HPP
#define MODIFIED_WAVELET_BASE_HPP

#include <cmath>

namespace sg
{

/**
 * Modified Wavelet basis functions.
 *
 */
template<class LT, class IT>
class modified_wavelet_base
{
    public:
        double eval(LT level, IT index, double p)
        {   //std::cout << "Level " << level <<" Index "<<index<<" Point "<<p<<" BasisValue ";
            if (level==1)
            {
                //std::cout<<1<<std::endl;
                return 1;
            }
            double sup= 1<<level;
            sup=1/sup;
            if (index==((1<<level)-1) && p > (1-1.5602*sup) )
            {   //std::cout<<0.5014+1.3803*(0.5602+p/sup-index)<<std::endl;
                return 0.5014+1.3803*(0.5602+p/sup-index);
            }
            if ( index==1 && p < 1.5602*sup )
            {   //std::cout<<0.5014+1.3803*(0.5602-p/sup+1)<<std::endl;
                return 0.5014+1.3803*(0.5602-p/sup+1);
            }
            double t=pow(p/sup-index,2);
            //std::cout<<(1-t)*exp(-0.9986*t)<<std::endl;
            return (1-t)*exp(-0.9986*t);

        }

};

}

#endif /* MODIFIED_WAVELET_BASE_HPP */
