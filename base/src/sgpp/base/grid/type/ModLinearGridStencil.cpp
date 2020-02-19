// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <sgpp/base/grid/type/ModLinearGridStencil.hpp>
#include <sgpp/base/exception/factory_exception.hpp>
#include <sgpp/base/operation/hash/common/basis/LinearModifiedBasis.hpp>

#include <sgpp/globaldef.hpp>


namespace sgpp {
namespace base {

ModLinearGridStencil::ModLinearGridStencil(std::istream& istr) :
    GridStencil(istr),
    generator(storage) {
}

ModLinearGridStencil::ModLinearGridStencil(size_t dim) :
    GridStencil(dim),
    generator(storage) {
}

ModLinearGridStencil::ModLinearGridStencil(BoundingBox& BB) :
    GridStencil(BB),
    generator(storage) {
}

ModLinearGridStencil::~ModLinearGridStencil() {
}

sgpp::base::GridType ModLinearGridStencil::getType() {
  return sgpp::base::GridType::ModLinearStencil;
}

SBasis& ModLinearGridStencil::getBasis() {
  throw factory_exception("Not implemented");
}

Grid* ModLinearGridStencil::unserialize(std::istream& istr) {
  return new ModLinearGridStencil(istr);
}

/**
 * Creates new GridGenerator
 * This must be changed if we add other storage types
 */
GridGenerator& ModLinearGridStencil::getGenerator() {
  return generator;
}


}  // namespace base
}  // namespace sgpp
