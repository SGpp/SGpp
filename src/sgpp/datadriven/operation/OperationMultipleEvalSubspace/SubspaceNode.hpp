#ifndef SUBSPACENODE_HPP
#define SUBSPACENODE_HPP

#include <vector>

#include <base/datatypes/DataVector.hpp>

using namespace std;
using namespace sg::base;

namespace sg {
namespace datadriven {
class SubspaceNode {
public:
    vector<size_t> level;
    vector<size_t> hInverse;
    size_t actualGridPointsOnLevel;
    vector<size_t> indices; //for stream computations
    size_t gridPointsOnLevel;

    SubspaceNode(DataVector &level, DataVector &hInverse, DataVector &index) {
        size_t dim = level.getSize();
        for (size_t i = 0; i < dim; i++) {
            this->level.push_back(static_cast<size_t>(level.get(i))); //TODO 3 ugly conversion
            this->hInverse.push_back(static_cast<size_t>(hInverse.get(i)));
            this->indices.push_back(static_cast<size_t>(index.get(i)));
        }
        this->actualGridPointsOnLevel = 1;
        gridPointsOnLevel = 1;
        for (size_t j = 0; j < dim; j++) {
            //cout << hInversePtr[j] << endl;
            //int dimTemp = hInverse[j];
            int dimTemp = static_cast<int>(this->hInverse[j]); //TODO ugly conversion
            dimTemp >>= 1; //skip even indices
            gridPointsOnLevel *= dimTemp;
        }
    }

    //increases number of grid points on the subspace
    void addGridPoint(DataVector &index) {
        size_t dim = index.getSize();
        for (size_t i = 0; i < dim; i++) {
            this->indices.push_back(static_cast<size_t>(index.get(i))); //TODO ugly conversion
        }
        this->actualGridPointsOnLevel += 1;
    }
};

}
}

#endif /* SUBSPACENODE_HPP */
