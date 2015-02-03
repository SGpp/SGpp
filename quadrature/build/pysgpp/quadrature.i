/******************************************************************************
* Copyright (C) 2009 Technische Universitaet Muenchen                         *
* This file is part of the SG++ project. For conditions of distribution and   *
* use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
******************************************************************************/

// The Good, i.e. without any modifications
%include "quadrature/src/sgpp/quadrature/Random.hpp"
%include "quadrature/src/sgpp/quadrature/sample/SampleGenerator.hpp"
%include "quadrature/src/sgpp/quadrature/sample/NaiveSampleGenerator.hpp"
%include "quadrature/src/sgpp/quadrature/sample/LatinHypercubeSampleGenerator.hpp"
//TODO franzefn
//%apply (int DIM1, long long int* IN_ARRAY1) {(size_t dimensions, long long int* strataPerDimension)};
%include "quadrature/src/sgpp/quadrature/sample/StratifiedSampleGenerator.hpp"
%include "quadrature/src/sgpp/quadrature/sample/HaltonSampleGenerator.hpp"

%apply std::string *INPUT { std::string& istr };

%apply unsigned int *OUTPUT { unsigned int& l, unsigned int& i };

%apply std::vector<std::pair<size_t, double> > *OUTPUT { std::vector<std::pair<size_t, double> >& result };
%apply std::vector<double> *INPUT { std::vector<double>& point };
