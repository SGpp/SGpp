// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <sgpp/base/grid/Grid.hpp>

#include <sgpp/base/exception/factory_exception.hpp>

#include <sgpp/base/grid/type/LinearGridStencil.hpp>
#include <sgpp/base/grid/type/ModLinearGridStencil.hpp>
#include <sgpp/base/grid/type/LinearGrid.hpp>
#include <sgpp/base/grid/type/LinearStretchedGrid.hpp>
#include <sgpp/base/grid/type/LinearL0BoundaryGrid.hpp>
#include <sgpp/base/grid/type/LinearClenshawCurtisGrid.hpp>
#include <sgpp/base/grid/type/ModLinearGrid.hpp>
#include <sgpp/base/grid/type/PolyGrid.hpp>
#include <sgpp/base/grid/type/PolyBoundaryGrid.hpp>
#include <sgpp/base/grid/type/ModPolyGrid.hpp>
#include <sgpp/base/grid/type/BsplineGrid.hpp>
#include <sgpp/base/grid/type/BsplineBoundaryGrid.hpp>
#include <sgpp/base/grid/type/BsplineClenshawCurtisGrid.hpp>
#include <sgpp/base/grid/type/ModBsplineGrid.hpp>
#include <sgpp/base/grid/type/ModBsplineClenshawCurtisGrid.hpp>
#include <sgpp/base/grid/type/FundamentalSplineGrid.hpp>
#include <sgpp/base/grid/type/ModFundamentalSplineGrid.hpp>
#include <sgpp/base/grid/type/WaveletGrid.hpp>
#include <sgpp/base/grid/type/WaveletBoundaryGrid.hpp>
#include <sgpp/base/grid/type/ModWaveletGrid.hpp>
#include <sgpp/base/grid/type/SquareRootGrid.hpp>
#include <sgpp/base/grid/type/PrewaveletGrid.hpp>
#include <sgpp/base/grid/type/PeriodicGrid.hpp>

#include <sgpp/base/grid/generation/functors/SurplusRefinementFunctor.hpp>

#include <sgpp/base/operation/BaseOpFactory.hpp>

#include <sgpp/globaldef.hpp>
#include <sgpp/base/grid/type/LinearTruncatedBoundaryGrid.hpp>
#include <sgpp/base/grid/type/LinearStretchedBoundaryGrid.hpp>
#include <sgpp/base/grid/type/LinearBoundaryGrid.hpp>

#include <utility>
#include <map>
#include <string>
#include <vector>

namespace sgpp {
namespace base {

Grid* Grid::createLinearGridStencil(size_t dim) {
  return new LinearGridStencil(dim);
}

Grid* Grid::createModLinearGridStencil(size_t dim) {
  return new ModLinearGridStencil(dim);
}

Grid* Grid::createLinearGrid(size_t dim) {
  return new LinearGrid(dim);
}

Grid* Grid::createLinearStretchedGrid(size_t dim) {
  return new LinearStretchedGrid(dim);
}

Grid* Grid::createLinearBoundaryGrid(size_t dim, level_t boundaryLevel) {
  if (boundaryLevel == 0) {
    return new LinearL0BoundaryGrid(dim);
  } else {
    return new LinearBoundaryGrid(dim, boundaryLevel);
  }
}

Grid* Grid::createLinearStretchedBoundaryGrid(size_t dim) {
  return new LinearStretchedBoundaryGrid(dim);
}

Grid* Grid::createLinearClenshawCurtisGrid(size_t dim) {
  return new LinearClenshawCurtisGrid(dim);
}

Grid* Grid::createModLinearGrid(size_t dim) {
  return new ModLinearGrid(dim);
}

Grid* Grid::createPolyGrid(size_t dim, size_t degree) {
  return new PolyGrid(dim, degree);
}

Grid* Grid::createPolyBoundaryGrid(size_t dim, size_t degree) {
  return new PolyBoundaryGrid(dim, degree);
}

Grid* Grid::createWaveletGrid(size_t dim) {
  return new WaveletGrid(dim);
}

Grid* Grid::createWaveletBoundaryGrid(size_t dim) {
  return new WaveletBoundaryGrid(dim);
}

Grid* Grid::createModWaveletGrid(size_t dim) {
  return new ModWaveletGrid(dim);
}

Grid* Grid::createBsplineGrid(size_t dim, size_t degree) {
  return new BsplineGrid(dim, degree);
}

Grid* Grid::createBsplineBoundaryGrid(size_t dim, size_t degree) {
  return new BsplineBoundaryGrid(dim, degree);
}

Grid* Grid::createBsplineClenshawCurtisGrid(size_t dim, size_t degree) {
  return new BsplineClenshawCurtisGrid(dim, degree);
}

Grid* Grid::createModBsplineGrid(size_t dim, size_t degree) {
  return new ModBsplineGrid(dim, degree);
}

Grid* Grid::createModBsplineClenshawCurtisGrid(size_t dim, size_t degree) {
  return new ModBsplineClenshawCurtisGrid(dim, degree);
}

Grid* Grid::createFundamentalSplineGrid(size_t dim, size_t degree) {
  return new FundamentalSplineGrid(dim, degree);
}

Grid* Grid::createModFundamentalSplineGrid(size_t dim, size_t degree) {
  return new ModFundamentalSplineGrid(dim, degree);
}

Grid* Grid::createSquareRootGrid(size_t dim) {
  return new SquareRootGrid(dim);
}

Grid* Grid::createPrewaveletGrid(size_t dim) {
  return new PrewaveletGrid(dim);
}

Grid* Grid::createLinearTruncatedBoundaryGrid(size_t dim) {
  return new LinearTruncatedBoundaryGrid(dim);
}

Grid* Grid::createModPolyGrid(size_t dim, size_t degree) {
  return new ModPolyGrid(dim, degree);
}

Grid* Grid::createPeriodicGrid(size_t dim) {
  return new PeriodicGrid(dim);
}

Grid* Grid::unserialize(const std::string& istr) {
  std::istringstream istream;
  istream.str(istr);

  return Grid::unserialize(istream);
}

Grid* Grid::unserialize(std::istream& istr) {
  std::string gridtype;
  istr >> gridtype;

  if (typeMap().count(gridtype) > 0) {
    // it calls a function pointer out of a map.
    return typeMap()[gridtype](istr);
  } else {
    // compose error message and throw factory_exception
    std::string errMsg = "factory_exception unserialize: unknown gridtype.\n";
    errMsg += "Got: '" + gridtype + "'; possible choices: ";
    factoryMap::iterator it;

    for (it = typeMap().begin(); it != typeMap().end(); it++) {
      errMsg += "'" + (*it).first + "' ";
    }

    throw factory_exception(errMsg.c_str());
  }

  return nullptr;
}

std::map<std::string, Grid::Factory>& Grid::typeMap() {
  // This is only executed once!
  static factoryMap* tMap = new factoryMap();

  if (tMap->size() == 0) {
/*
 * Insert factories here. This methods may NOT read the grid type.
 * This map takes a string, function pointer pair.
 */
#ifdef _WIN32
    tMap->insert(std::pair<std::string, Grid::Factory>("NULL", Grid::nullFactory));
    tMap->insert(std::pair<std::string, Grid::Factory>("linear", LinearGrid::unserialize));
    tMap->insert(
        std::pair<std::string, Grid::Factory>("linearStretched", LinearStretchedGrid::unserialize));
    tMap->insert(std::pair<std::string, Grid::Factory>("linearL0Boundary",
                                                       LinearL0BoundaryGrid::unserialize));
    tMap->insert(
        std::pair<std::string, Grid::Factory>("linearstencil", LinearGridStencil::unserialize));
    tMap->insert(std::pair<std::string, Grid::Factory>("modlinearstencil",
                                                       ModLinearGridStencil::unserialize));
    tMap->insert(
        std::pair<std::string, Grid::Factory>("linearBoundary", LinearBoundaryGrid::unserialize));
    tMap->insert(std::pair<std::string, Grid::Factory>("linearStretchedBoundary",
                                                       LinearStretchedBoundaryGrid::unserialize));
    tMap->insert(std::pair<std::string, Grid::Factory>("linearClenshawCurtis",
                                                       LinearClenshawCurtisGrid::unserialize));
    tMap->insert(std::pair<std::string, Grid::Factory>("modlinear", ModLinearGrid::unserialize));
    tMap->insert(std::pair<std::string, Grid::Factory>("poly", PolyGrid::unserialize));
    tMap->insert(
        std::pair<std::string, Grid::Factory>("polyBoundary", PolyBoundaryGrid::unserialize));
    tMap->insert(std::pair<std::string, Grid::Factory>("modpoly", ModPolyGrid::unserialize));
    tMap->insert(std::pair<std::string, Grid::Factory>("wavelet", WaveletGrid::unserialize));
    tMap->insert(
        std::pair<std::string, Grid::Factory>("waveletBoundary", WaveletBoundaryGrid::unserialize));
    tMap->insert(std::pair<std::string, Grid::Factory>("modWavelet", ModWaveletGrid::unserialize));
    tMap->insert(std::pair<std::string, Grid::Factory>("bspline", BsplineGrid::unserialize));
    tMap->insert(
        std::pair<std::string, Grid::Factory>("bsplineBoundary", BsplineBoundaryGrid::unserialize));
    tMap->insert(std::pair<std::string, Grid::Factory>("bsplineClenshawCurtis",
                                                       BsplineClenshawCurtisGrid::unserialize));
    tMap->insert(std::pair<std::string, Grid::Factory>("modBspline", ModBsplineGrid::unserialize));
    tMap->insert(std::pair<std::string, Grid::Factory>("fundamentalSpline",
                                                       FundamentalSplineGrid::unserialize));
    tMap->insert(std::pair<std::string, Grid::Factory>("modFundamentalSpline",
                                                       ModFundamentalSplineGrid::unserialize));
    tMap->insert(std::pair<std::string, Grid::Factory>("modBsplineClenshawCurtis",
                                                       ModBsplineClenshawCurtisGrid::unserialize));
    tMap->insert(std::pair<std::string, Grid::Factory>("prewavelet", PrewaveletGrid::unserialize));
    tMap->insert(std::pair<std::string, Grid::Factory>("periodic", PeriodicGrid::unserialize));
    tMap->insert(std::pair<std::string, Grid::Factory>("linearTruncatedBoundary",
                                                       LinearTruncatedBoundaryGrid::unserialize));
#else
    tMap->insert(std::make_pair("NULL", Grid::nullFactory));
    tMap->insert(std::make_pair("linear", LinearGrid::unserialize));
    tMap->insert(std::make_pair("linearStretched", LinearStretchedGrid::unserialize));
    tMap->insert(std::make_pair("linearL0Boundary", LinearL0BoundaryGrid::unserialize));
    tMap->insert(std::make_pair("linearstencil", LinearGridStencil::unserialize));
    tMap->insert(std::make_pair("modlinearstencil", ModLinearGridStencil::unserialize));
    tMap->insert(std::make_pair("linearBoundary", LinearBoundaryGrid::unserialize));
    tMap->insert(
        std::make_pair("linearStretchedBoundary", LinearStretchedBoundaryGrid::unserialize));
    tMap->insert(std::make_pair("linearClenshawCurtis", LinearClenshawCurtisGrid::unserialize));
    tMap->insert(std::make_pair("modlinear", ModLinearGrid::unserialize));
    tMap->insert(std::make_pair("poly", PolyGrid::unserialize));
    tMap->insert(std::make_pair("polyBoundary", PolyBoundaryGrid::unserialize));
    tMap->insert(std::make_pair("modpoly", ModPolyGrid::unserialize));
    tMap->insert(std::make_pair("wavelet", WaveletGrid::unserialize));
    tMap->insert(std::make_pair("waveletBoundary", WaveletBoundaryGrid::unserialize));
    tMap->insert(std::make_pair("modWavelet", ModWaveletGrid::unserialize));
    tMap->insert(std::make_pair("bspline", BsplineGrid::unserialize));
    tMap->insert(std::make_pair("bsplineBoundary", BsplineBoundaryGrid::unserialize));
    tMap->insert(std::make_pair("bsplineClenshawCurtis", BsplineClenshawCurtisGrid::unserialize));
    tMap->insert(std::make_pair("modBspline", ModBsplineGrid::unserialize));
    tMap->insert(std::make_pair("fundamentalSpline", FundamentalSplineGrid::unserialize));
    tMap->insert(std::make_pair("modFundamentalSpline", ModFundamentalSplineGrid::unserialize));
    tMap->insert(
        std::make_pair("modBsplineClenshawCurtis", ModBsplineClenshawCurtisGrid::unserialize));
    tMap->insert(std::make_pair("prewavelet", PrewaveletGrid::unserialize));
    tMap->insert(std::make_pair("periodic", PeriodicGrid::unserialize));
    tMap->insert(
        std::make_pair("linearTruncatedBoundary", LinearTruncatedBoundaryGrid::unserialize));
#endif
  }

  return *tMap;
}

std::map<sgpp::base::GridType, std::string>& Grid::typeVerboseMap() {
  // This is only executed once!
  static gridTypeVerboseMap* verboseMap = new gridTypeVerboseMap();

  if (verboseMap->size() == 0) {
/*
 * Insert strings here.
 */
#ifdef _WIN32
    verboseMap->insert(std::pair<sgpp::base::GridType, std::string>(GridType::Linear, "linear"));
    verboseMap->insert(
        std::pair<sgpp::base::GridType, std::string>(GridType::LinearStretched, "linearStretched"));
    verboseMap->insert(std::pair<sgpp::base::GridType, std::string>(GridType::LinearL0Boundary,
                                                                    "linearL0Boundary"));
    verboseMap->insert(
        std::pair<sgpp::base::GridType, std::string>(GridType::LinearStencil, "linearstencil"));
    verboseMap->insert(std::pair<sgpp::base::GridType, std::string>(GridType::ModLinearStencil,
                                                                    "modlinearstencil"));
    verboseMap->insert(
        std::pair<sgpp::base::GridType, std::string>(GridType::LinearBoundary, "linearBoundary"));
    verboseMap->insert(std::pair<sgpp::base::GridType, std::string>(
        GridType::LinearStretchedBoundary, "linearStretchedBoundary"));
    verboseMap->insert(std::pair<sgpp::base::GridType, std::string>(GridType::LinearClenshawCurtis,
                                                                    "linearClenshawCurtis"));
    verboseMap->insert(
        std::pair<sgpp::base::GridType, std::string>(GridType::ModLinear, "modlinear"));
    verboseMap->insert(std::pair<sgpp::base::GridType, std::string>(GridType::Poly, "poly"));
    verboseMap->insert(
        std::pair<sgpp::base::GridType, std::string>(GridType::PolyBoundary, "polyBoundary"));
    verboseMap->insert(std::pair<sgpp::base::GridType, std::string>(GridType::ModPoly, "modpoly"));
    verboseMap->insert(std::pair<sgpp::base::GridType, std::string>(GridType::Wavelet, "wavelet"));
    verboseMap->insert(
        std::pair<sgpp::base::GridType, std::string>(GridType::WaveletBoundary, "waveletBoundary"));
    verboseMap->insert(
        std::pair<sgpp::base::GridType, std::string>(GridType::ModWavelet, "modWavelet"));
    verboseMap->insert(std::pair<sgpp::base::GridType, std::string>(GridType::Bspline, "bspline"));
    verboseMap->insert(
        std::pair<sgpp::base::GridType, std::string>(GridType::BsplineBoundary, "bsplineBoundary"));
    verboseMap->insert(std::pair<sgpp::base::GridType, std::string>(GridType::BsplineClenshawCurtis,
                                                                    "bsplineClenshawCurtis"));
    verboseMap->insert(
        std::pair<sgpp::base::GridType, std::string>(GridType::ModBspline, "modBspline"));
    verboseMap->insert(std::pair<sgpp::base::GridType, std::string>(GridType::FundamentalSpline,
                                                                    "fundamentalSpline"));
    verboseMap->insert(std::pair<sgpp::base::GridType, std::string>(GridType::ModFundamentalSpline,
                                                                    "modFundamentalSpline"));
    verboseMap->insert(std::pair<sgpp::base::GridType, std::string>(
        GridType::ModBsplineClenshawCurtis, "modBsplineClenshawCurtis"));
    verboseMap->insert(
        std::pair<sgpp::base::GridType, std::string>(GridType::Prewavelet, "prewavelet"));
    verboseMap->insert(
        std::pair<sgpp::base::GridType, std::string>(GridType::Periodic, "periodic"));
    verboseMap->insert(std::pair<sgpp::base::GridType, std::string>(
        GridType::LinearTruncatedBoundary, "linearTruncatedBoundary"));
#else
    verboseMap->insert(std::make_pair(GridType::Linear, "linear"));
    verboseMap->insert(std::make_pair(GridType::LinearStretched, "linearStretched"));
    verboseMap->insert(std::make_pair(GridType::LinearL0Boundary, "linearL0Boundary"));
    verboseMap->insert(std::make_pair(GridType::LinearStencil, "linearstencil"));
    verboseMap->insert(std::make_pair(GridType::ModLinearStencil, "modlinearstencil"));
    verboseMap->insert(std::make_pair(GridType::LinearBoundary, "linearBoundary"));
    verboseMap->insert(
        std::make_pair(GridType::LinearStretchedBoundary, "linearStretchedBoundary"));
    verboseMap->insert(std::make_pair(GridType::LinearClenshawCurtis, "linearClenshawCurtis"));
    verboseMap->insert(std::make_pair(GridType::ModLinear, "modlinear"));
    verboseMap->insert(std::make_pair(GridType::Poly, "poly"));
    verboseMap->insert(std::make_pair(GridType::PolyBoundary, "polyBoundary"));
    verboseMap->insert(std::make_pair(GridType::ModPoly, "modpoly"));
    verboseMap->insert(std::make_pair(GridType::Wavelet, "wavelet"));
    verboseMap->insert(std::make_pair(GridType::WaveletBoundary, "waveletBoundary"));
    verboseMap->insert(std::make_pair(GridType::ModWavelet, "modWavelet"));
    verboseMap->insert(std::make_pair(GridType::Bspline, "bspline"));
    verboseMap->insert(std::make_pair(GridType::BsplineBoundary, "bsplineBoundary"));
    verboseMap->insert(std::make_pair(GridType::BsplineClenshawCurtis, "bsplineClenshawCurtis"));
    verboseMap->insert(std::make_pair(GridType::ModBspline, "modBspline"));
    verboseMap->insert(std::make_pair(GridType::FundamentalSpline, "fundamentalSpline"));
    verboseMap->insert(std::make_pair(GridType::ModFundamentalSpline, "modFundamentalSpline"));
    verboseMap->insert(
        std::make_pair(GridType::ModBsplineClenshawCurtis, "modBsplineClenshawCurtis"));
    verboseMap->insert(std::make_pair(GridType::Prewavelet, "prewavelet"));
    verboseMap->insert(std::make_pair(GridType::Periodic, "periodic"));
    verboseMap->insert(
        std::make_pair(GridType::LinearTruncatedBoundary, "linearTruncatedBoundary"));
#endif
  }

  return *verboseMap;
}

/**
 * Factory for everything we don't know.
 */
Grid* Grid::nullFactory(std::istream&) {
  throw factory_exception("factory_exeception unserialize: unsupported gridtype");
  return nullptr;
}

Grid::Grid(std::istream& istr) : storage(istr) {}

Grid::Grid(size_t dim) : storage(dim) {}

Grid::Grid(BoundingBox& boundingBox) : storage(boundingBox) {}

Grid::Grid(Stretching& stretching) : storage(stretching) {}

Grid::~Grid() {}

GridStorage& Grid::getStorage() { return storage; }

BoundingBox& Grid::getBoundingBox() { return *storage.getBoundingBox(); }

Stretching& Grid::getStretching() { return *storage.getStretching(); }

void Grid::setBoundingBox(BoundingBox& boundingBox) { storage.setBoundingBox(boundingBox); }

void Grid::setStretching(Stretching& stretching) { storage.setStretching(stretching); }

void Grid::serialize(std::string& ostr, int version) {
  std::ostringstream ostream;
  this->serialize(ostream, version);

  ostr = ostream.str();
}

std::string Grid::serialize(int version) {
  std::ostringstream ostream;
  this->serialize(ostream, version);

  return ostream.str();
}

void Grid::serialize(std::ostream& ostr, int version) {
  ostr << typeVerboseMap()[this->getType()] << std::endl;
  storage.serialize(ostr, version);
}

void Grid::refine(DataVector& vector, int numOfPoints) {
  SurplusRefinementFunctor functor(vector, numOfPoints);
  getGenerator().refine(functor);
}

void Grid::insertPoint(size_t dim, unsigned int levels[], unsigned int indices[], bool isLeaf) {
  // create HashGridPoint object for the point
  GridPoint pointIndex(dim);

  for (unsigned int i = 0; i < dim - 1; i++) {
    pointIndex.push(i, levels[i], indices[i]);
  }

  // insert last level/index and hash
  pointIndex.set(dim - 1, levels[dim - 1], indices[dim - 1], isLeaf);
  // insert point to the GridStorage
  storage.insert(pointIndex);
}

size_t Grid::getDimension() { return storage.getDimension(); }

size_t Grid::getSize() { return storage.getSize(); }

std::vector<size_t> Grid::getAlgorithmicDimensions() { return storage.getAlgorithmicDimensions(); }

void Grid::setAlgorithmicDimensions(std::vector<size_t> newAlgoDims) {
  this->storage.setAlgorithmicDimensions(newAlgoDims);
}

}  // namespace base
}  // namespace sgpp
