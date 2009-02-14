# This file is part of sg++, a program package making use of spatially adaptive sparse grids to solve numerical problems.
#
# Copyright (C) 2007  Joerg Blank (blankj@in.tum.de)
# 
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU Lesser General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
# 
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# GNU Lesser General Public License for more details.
# 
# You should have received a copy of the GNU Lesser General Public License
# along with pyclass. If not, see <http://www.gnu.org/licenses/>.
#


%module pysgpp

%include "stl.i"
%include "std_vector.i"
%include "std_pair.i"

%include "typemaps.i"

%include "exception.i"

%exception {
  try {
    $action
  } catch (const std::exception& e) {
    SWIG_exception(SWIG_RuntimeError, e.what());
  }
}


namespace std {
	%template(DoubleVector) vector<double>;
	%template(IndexValPair) pair<size_t, double>;
	%template(IndexValVector) vector<pair<size_t, double> >;
}

// This should include all necessary header files
%{
#include "sgpp.hpp"
%}

// The Good,
%include "sg++/GridStorage.hpp"
%include "Operations.i"

%include "GridFactory.i"

// the Bad

%include "DataVector.i"

// and the rest

%include "sg++/sgpp.hpp"
%include "sg++/GridIndex.hpp"
%include "sg++/algorithms.hpp"
%include "sg++/generation.hpp"
%include "sg++/base.hpp"

%apply std::string *INPUT { std::string& istr };

%apply unsigned int *OUTPUT { unsigned int& l, unsigned int& i };

%template(GridIndex) sg::HashGridIndex<unsigned int, unsigned int>;
%template(GridStorage) sg::HashGridStorage<sg::GridIndex>;

%template(SLinearBase) sg::linear_base<unsigned int, unsigned int>;
%template(SModLinearBase) sg::modified_linear_base<unsigned int, unsigned int>;
%template(SPolyBase) sg::poly_base<unsigned int, unsigned int>;
%template(SModPolyBase) sg::modified_poly_base<unsigned int, unsigned int>;

%apply std::vector<std::pair<size_t, double> > *OUTPUT { std::vector<std::pair<size_t, double> >& result };
%apply std::vector<double> *INPUT { std::vector<double>& point }; 
%template(SGetAffectedBasisFunctions) sg::GetAffectedBasisFunctions<sg::SLinearBase>;

%template(test_dataset_linear) sg::test_dataset<sg::SLinearBase>;
%template(test_dataset_modlin) sg::test_dataset<sg::SModLinearBase>;
