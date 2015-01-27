// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at 
// sgpp.sparsegrids.org

#include <sgpp/base/grid/Grid.hpp>
#include <sgpp/base/grid/type/TruncatedTrapezoidGrid.hpp>

#include <sgpp/base/grid/generation/TruncatedTrapezoidGridGenerator.hpp>
//#include <sgpp/base/operation/hash/OperationEvalLinearBoundary.hpp>
//#include <sgpp/base/operation/hash/OperationHierarchisationLinearBoundary.hpp>

#include <sgpp/base/basis/linear/boundary/LinearBoundaryBasis.hpp>


#include <iostream>

#include <sgpp/globaldef.hpp>


namespace SGPP {
  namespace base {

    TruncatedTrapezoidGrid::TruncatedTrapezoidGrid(std::istream& istr) : Grid(istr) {

    }

    TruncatedTrapezoidGrid::TruncatedTrapezoidGrid(size_t dim) {
      this->storage = new GridStorage(dim);
    }

    TruncatedTrapezoidGrid::TruncatedTrapezoidGrid(BoundingBox& BB) {
      this->storage = new GridStorage(BB);
    }

    TruncatedTrapezoidGrid::~TruncatedTrapezoidGrid() {
    }

    const char* TruncatedTrapezoidGrid::getType() {
      return "TruncatedTrapezoid";
    }

    const SBasis& TruncatedTrapezoidGrid::getBasis(){
		static SLinearBoundaryBase basis;
		return basis;
	}

    Grid* TruncatedTrapezoidGrid::unserialize(std::istream& istr) {
      return new TruncatedTrapezoidGrid(istr);
    }
    /**
     * Creates new GridGenerator
     * This must be changed if we add other storage types
     */
    GridGenerator* TruncatedTrapezoidGrid::createGridGenerator() {
      return new TruncatedTrapezoidGridGenerator(this->storage);
    }
    //OperationHierarchisation* TruncatedTrapezoidGrid::createOperationHierarchisation()
    //{
    //  return new OperationHierarchisationLinearBoundary(this->storage);
    //}
    //OperationEval* TruncatedTrapezoidGrid::createOperationEval()
    //{
    //  return new OperationEvalLinearBoundary(this->storage);
    //}

    //OperationConvert* TruncatedTrapezoidGrid::createOperationConvert()
    //{
    //  throw factory_exception("Unsupported operation");
    //}

  }
}
