#pragma once

#include <vector>

#include <base/datatypes/DataVector.hpp>

using namespace std;
using namespace sg::base;

namespace sg {
  namespace parallel {
    class X86VectorizedSubspaceNode {
    public:
      vector<uint32_t> level;
      vector<uint32_t> hInverse;
      uint32_t actualGridPointsOnLevel;
      vector<uint32_t> indices; //for stream computations
      
      X86VectorizedSubspaceNode(DataVector &level, DataVector &hInverse, DataVector &index) {
	size_t dim = level.getSize();
	for (size_t i = 0; i < dim; i++) {
	  this->level.push_back(level.get(i));
	  this->hInverse.push_back(hInverse.get(i));
	  this->indices.push_back(index.get(i));
	}
	this->actualGridPointsOnLevel = 1;
      }
      
      //increases number of grid points on the subspace
      void addGridPoint(DataVector &index) {
	size_t dim = index.getSize();
	for (size_t i = 0; i < dim; i++) {
	  this->indices.push_back(index.get(i));
	}
	this->actualGridPointsOnLevel += 1;
      }
    };

  }
}
