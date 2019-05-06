// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#ifndef SAMPLERTYPES_HPP_
#define SAMPLERTYPES_HPP_

#include <sgpp/globaldef.hpp>

namespace sgpp {
namespace quadrature {

enum class SamplerTypes { Naive, Stratified, LatinHypercube, Halton };

}  // namespace quadrature
}  // namespace sgpp

#endif /* SAMPLERTYPES_HPP_ */
