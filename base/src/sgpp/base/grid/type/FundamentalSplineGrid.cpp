// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <sgpp/base/grid/Grid.hpp>
#include <sgpp/base/grid/GridStorage.hpp>
#include <sgpp/base/grid/type/FundamentalSplineGrid.hpp>

#include <sgpp/base/grid/generation/StandardGridGenerator.hpp>

#include <sgpp/base/exception/factory_exception.hpp>


#include <sgpp/globaldef.hpp>


namespace SGPP {
namespace base {

FundamentalSplineGrid::FundamentalSplineGrid(std::istream& istr) :
  Grid(istr),
  degree(1 << 16),
  basis_(NULL) {
  istr >> degree;
}


FundamentalSplineGrid::FundamentalSplineGrid(size_t dim, size_t degree) :
  Grid(dim),
  degree(degree),
  basis_(NULL) {
}

FundamentalSplineGrid::~FundamentalSplineGrid() {
  if (basis_ != NULL) {
    delete basis_;
  }
}

SGPP::base::GridType FundamentalSplineGrid::getType() {
  return SGPP::base::GridType::FundamentalSpline;
}

const SBasis& FundamentalSplineGrid::getBasis() {
  if (basis_ == NULL) {
    basis_ = new SFundamentalSplineBase(degree);
  }

  return *basis_;
}

size_t FundamentalSplineGrid::getDegree() {
  return this->degree;
}

std::unique_ptr<Grid> FundamentalSplineGrid::unserialize(std::istream& istr) {
  return std::unique_ptr<Grid>(new FundamentalSplineGrid(istr));
}

void FundamentalSplineGrid::serialize(std::ostream& ostr) {
  this->Grid::serialize(ostr);
  ostr << degree << std::endl;
}

/**
 * Creates new GridGenerator
 * This must be changed if we add other storage types
 */
std::unique_ptr<GridGenerator> FundamentalSplineGrid::createGridGenerator() {
  return std::unique_ptr<GridGenerator>(new StandardGridGenerator(storage));
}

}  // namespace base
}  // namespace SGPP
