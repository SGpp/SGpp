// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <sgpp/base/grid/type/LinearGridStencil.hpp>
#include <sgpp/base/operation/hash/common/basis/LinearBasis.hpp>

#include <sgpp/base/exception/factory_exception.hpp>

#include <sgpp/globaldef.hpp>

#include <iostream>
#include <exception>

namespace sgpp {
namespace base {

LinearGridStencil::LinearGridStencil(std::istream& istr) :
    GridStencil(istr),
    generator(storage) {
}

LinearGridStencil::LinearGridStencil(size_t dim) :
    GridStencil(dim),
    generator(storage) {
}

LinearGridStencil::LinearGridStencil(BoundingBox& BB) :
    GridStencil(BB),
    generator(storage) {
}

LinearGridStencil::~LinearGridStencil() {
}

sgpp::base::GridType LinearGridStencil::getType() {
  return sgpp::base::GridType::LinearStencil;
}

SBasis& LinearGridStencil::getBasis() {
  throw factory_exception("Not implemented");
}

Grid* LinearGridStencil::unserialize(std::istream& istr) {
  return new LinearGridStencil(istr);
}

/**
 * Creates new GridGenerator
 * This must be changed if we add other storage types
 */
GridGenerator& LinearGridStencil::getGenerator() {
  return generator;
}


}  // namespace base
}  // namespace sgpp
