// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <sgpp/combigrid/MultiFunction.hpp>

namespace sgpp {
namespace combigrid {

MultiFunction::MultiFunction(double (*ptr)(const base::DataVector&)) : func(ptr) {}

double MultiFunction::operator()(const base::DataVector& vec) { return func(vec); }

double MultiFunction::call(const base::DataVector& vec) { return func(vec); }

} /* namespace combigrid */
} /* namespace sgpp*/
