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
#include <sgpp/base/grid/type/LinearBoundaryGrid.hpp>
#include <sgpp/base/grid/type/LinearTrapezoidBoundaryGrid.hpp>
#include <sgpp/base/grid/type/LinearStretchedTrapezoidBoundaryGrid.hpp>
#include <sgpp/base/grid/type/ModLinearGrid.hpp>
#include <sgpp/base/grid/type/PolyGrid.hpp>
#include <sgpp/base/grid/type/ModPolyGrid.hpp>
#include <sgpp/base/grid/type/ModBsplineGrid.hpp>
#include <sgpp/base/grid/type/ModWaveletGrid.hpp>
#include <sgpp/base/grid/type/SquareRootGrid.hpp>
#include <sgpp/base/grid/type/TruncatedTrapezoidGrid.hpp>
#include <sgpp/base/grid/type/PrewaveletGrid.hpp>
#include <sgpp/base/grid/type/PeriodicGrid.hpp>

#include <sgpp/base/grid/generation/functors/SurplusRefinementFunctor.hpp>

#include <sgpp/base/operation/BaseOpFactory.hpp>


#include <sgpp/globaldef.hpp>


namespace SGPP {
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

    Grid* Grid::createLinearBoundaryGrid(size_t dim) {
      return new LinearBoundaryGrid(dim);
    }

    Grid* Grid::createLinearTrapezoidBoundaryGrid(size_t dim) {
      return new LinearTrapezoidBoundaryGrid(dim);
    }

    Grid* Grid::createLinearStretchedTrapezoidBoundaryGrid(size_t dim) {
      return new LinearStretchedTrapezoidBoundaryGrid(dim);
    }

    Grid* Grid::createModLinearGrid(size_t dim) {
      return new ModLinearGrid(dim);
    }

    Grid* Grid::createPolyGrid(size_t dim, size_t degree) {
      return new PolyGrid(dim, degree);
    }

    Grid* Grid::createModWaveletGrid(size_t dim) {
      return new ModWaveletGrid(dim);
    }

    Grid* Grid::createModBsplineGrid(size_t dim, size_t degree) {
      return new ModBsplineGrid(dim, degree);
    }

    Grid* Grid::createSquareRootGrid(size_t dim) {
      return new SquareRootGrid(dim);
    }

    Grid* Grid::createPrewaveletGrid(size_t dim) {
      return new PrewaveletGrid(dim);
    }

    Grid* Grid::createTruncatedTrapezoidGrid(size_t dim) {
      return new TruncatedTrapezoidGrid(dim);
    }

    //OperationMatrix* Grid::createOperationIdentity()
    //{
    //  return new OperationIdentity();
    //}

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
        std::string errMsg = "factory_exeception unserialize: unkown gridtype.\n";
        errMsg += "Got: '" + gridtype + "'; possible choices: ";
        factoryMap::iterator it;

        for (it = typeMap().begin(); it != typeMap().end(); it++) {
          errMsg += "'" + (*it).first + "' ";
        }

        throw factory_exception(errMsg.c_str());
      }

      return NULL;
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
        tMap->insert(std::pair<std::string, Grid::Factory>("linearStretched", LinearStretchedGrid::unserialize));
        tMap->insert(std::pair<std::string, Grid::Factory>("linearBoundary", LinearBoundaryGrid::unserialize));
        tMap->insert(std::pair<std::string, Grid::Factory>("linearstencil", LinearGridStencil::unserialize));
        tMap->insert(std::pair<std::string, Grid::Factory>("modlinearstencil", ModLinearGridStencil::unserialize));
        tMap->insert(std::pair<std::string, Grid::Factory>("linearTrapezoidBoundary", LinearTrapezoidBoundaryGrid::unserialize));
        tMap->insert(std::pair<std::string, Grid::Factory>("linearStretchedTrapezoidBoundary", LinearStretchedTrapezoidBoundaryGrid::unserialize));
        tMap->insert(std::pair<std::string, Grid::Factory>("modlinear", ModLinearGrid::unserialize));
        tMap->insert(std::pair<std::string, Grid::Factory>("poly", PolyGrid::unserialize));
        tMap->insert(std::pair<std::string, Grid::Factory>("modpoly", ModPolyGrid::unserialize));
        tMap->insert(std::pair<std::string, Grid::Factory>("modWavelet", ModWaveletGrid::unserialize));
        tMap->insert(std::pair<std::string, Grid::Factory>("modBspline", ModBsplineGrid::unserialize));
        tMap->insert(std::pair<std::string, Grid::Factory>("prewavelet", PrewaveletGrid::unserialize));
        tMap->insert(std::pair<std::string, Grid::Factory>("periodic", PeriodicGrid::unserialize));
#else
        tMap->insert(std::make_pair("NULL", Grid::nullFactory));
        tMap->insert(std::make_pair("linear", LinearGrid::unserialize));
        tMap->insert(std::make_pair("linearStretched", LinearStretchedGrid::unserialize));
        tMap->insert(std::make_pair("linearBoundary", LinearBoundaryGrid::unserialize));
        tMap->insert(std::make_pair("linearstencil", LinearGridStencil::unserialize));
        tMap->insert(std::make_pair("modlinearstencil", ModLinearGridStencil::unserialize));
        tMap->insert(std::make_pair("linearTrapezoidBoundary", LinearTrapezoidBoundaryGrid::unserialize));
        tMap->insert(std::make_pair("linearStretchedTrapezoidBoundary", LinearStretchedTrapezoidBoundaryGrid::unserialize));
        tMap->insert(std::make_pair("modlinear", ModLinearGrid::unserialize));
        tMap->insert(std::make_pair("poly", PolyGrid::unserialize));
        tMap->insert(std::make_pair("modpoly", ModPolyGrid::unserialize));
        tMap->insert(std::make_pair("modWavelet", ModWaveletGrid::unserialize));
        tMap->insert(std::make_pair("modBspline", ModBsplineGrid::unserialize));
        tMap->insert(std::make_pair("prewavelet", PrewaveletGrid::unserialize));
        tMap->insert(std::make_pair("periodic", PeriodicGrid::unserialize));
#endif
      }

      return *tMap;
    }

    /**
     * Factory for everything we don't know.
     */
    Grid* Grid::nullFactory(std::istream&) {
      throw factory_exception("factory_exeception unserialize: unsupported gridtype");
      return NULL;
    }

    Grid::Grid(std::istream& istr) : storage(NULL) {
      int hasStorage;
      istr >> hasStorage;

      if (hasStorage == 1) {
        storage = new GridStorage(istr);
      }
    }

    Grid::Grid() : storage(NULL) {
    }

    Grid::~Grid() {
      if (storage != NULL) {
        delete storage;
      }

      if (evalOp != NULL) {
        delete evalOp;
        evalOp = NULL;
      }

    }

    GridStorage* Grid::getStorage() {
      return this->storage;
    }

    BoundingBox* Grid::getBoundingBox() {
      return this->storage->getBoundingBox();
    }

    Stretching* Grid::getStretching() {
      return this->storage->getStretching();
    }

    void Grid::setBoundingBox(BoundingBox& bb) {
      this->storage->setBoundingBox(bb);
    }

    void Grid::setStretching(Stretching& bb) {
      this->storage->setStretching(bb);
    }

    void Grid::serialize(std::string& ostr) {
      std::ostringstream ostream;
      this->serialize(ostream);

      ostr = ostream.str();
    }

    std::string Grid::serialize() {
      std::ostringstream ostream;
      this->serialize(ostream);

      return ostream.str();
    }

    void Grid::serialize(std::ostream& ostr) {
      ostr << this->getType() << std::endl;

      if (storage != NULL) {
        ostr << "1" << std::endl;
        storage->serialize(ostr);
      } else {
        ostr << "0" << std::endl;
      }
    }

    void Grid::refine(DataVector* vector, int numOfPoints) {
      // @todo (khakhutv) (low) different refinemente Functors
      this->createGridGenerator()->refine(new SurplusRefinementFunctor(vector, numOfPoints));
    }

    OperationEval* Grid::evalOp(NULL);

    double Grid::eval(DataVector& alpha, DataVector& point) {
      if (this->evalOp == NULL) this->evalOp = SGPP::op_factory::createOperationEval(*this);

      return this->evalOp->eval(alpha, point);
    }

    void Grid::insertPoint(size_t dim, unsigned int levels[], unsigned int indices[], bool isLeaf) {
      //create HashGridIndex object for the point
      GridIndex pointIndex = new GridIndex(dim);

      for (unsigned int i = 0; i < dim - 1; i++) {
        pointIndex.push(i, levels[i], indices[i]);
      }

      //insert last level/index and hash
      pointIndex.set(dim - 1, levels[dim - 1], indices[dim - 1], isLeaf);
      //insert point to the GridStorage
      storage->insert(pointIndex);
    }

    size_t Grid::getSize() {
      return this->storage->size();
    }

    std::vector<size_t> Grid::getAlgorithmicDimensions() {
      return this->storage->getAlgorithmicDimensions();
    }

    void Grid::setAlgorithmicDimensions(std::vector<size_t> newAlgoDims) {
      this->storage->setAlgorithmicDimensions(newAlgoDims);
    }

  }
}