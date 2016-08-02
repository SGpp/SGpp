// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <sgpp/combigrid/SingleFunction.hpp>

namespace sgpp {
namespace combigrid {

SingleFunction::SingleFunction(double (*ptr)(double)) : func(ptr) {}

double SingleFunction::operator()(double param) { return func(param); }

double SingleFunction::call(double param) { return func(param); }

} /* namespace combigrid */
} /* namespace sgpp*/
