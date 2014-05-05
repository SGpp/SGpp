/* ****************************************************************************
* Copyright (C) 2009 Technische Universitaet Muenchen                         *
* This file is part of the SG++ project. For conditions of distribution and   *
* use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
**************************************************************************** */
// @author Julian Valentin (julian.valentin@stud.mathematik.uni-stuttgart.de)

#ifndef BASIS_HPP
#define BASIS_HPP

namespace sg
{
namespace base
{

template<class LT, class IT>
class Basis
{
public:
    virtual ~Basis() {};
    
    virtual double eval(LT level, IT index, double x) = 0;
};

typedef Basis<unsigned int, unsigned int> SBase;

}
}

#endif /* BASIS_HPP */
