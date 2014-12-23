#pragma once

#include <vector>
#include <limits>
#include <omp.h>

#include <base/datatypes/DataVector.hpp>

#include "../../OperationMultipleEvalSubspace/perfect/ListSubspaceTuple.hpp"

using namespace std;
using namespace sg::base;

namespace sg {
namespace parallel {
class X86PerfectSubspaceNode {
public:
    enum SubspaceType {
        NOT_SET, ARRAY, LIST
    };

    vector<uint32_t> level;
    vector<uint32_t> hInverse;
    uint32_t gridPointsOnLevel;
    uint32_t existingGridPointsOnLevel;
    SubspaceType type;
    vector<uint32_t> indices; //for list representation (and future streaming subspaces
    vector<ListSubspaceTuple> indexFlatSurplusPairs;
    double *subspaceArray;
    omp_lock_t subspaceLock;

    uint32_t jumpTargetIndex;
    uint32_t flatLevel;
    //if the next node is reached, which dims are to be recomputed
    uint32_t nextDiff;
    //if a jump is performed, which components are to be recomputed
    uint32_t jumpDiff;

    X86PerfectSubspaceNode(DataVector &level, size_t flatLevel, DataVector &hInverse, DataVector &index) {
        size_t dim = level.getSize();
        for (size_t i = 0; i < dim; i++) {
            this->level.push_back(level.get(i));
            this->hInverse.push_back(hInverse.get(i));
            this->indices.push_back(index.get(i));
        }
        // exactly one gp is added in the loop above
        this->existingGridPointsOnLevel = 1;
        this->flatLevel = flatLevel;
        this->type = NOT_SET;

        this->gridPointsOnLevel = 1;
        for (size_t j = 0; j < dim; j++) {
            int dimTemp = hInverse[j];
            dimTemp >>= 1; //skip even indices
            this->gridPointsOnLevel *= dimTemp;
        }
        this->subspaceArray = nullptr;

        //initalize other member variables with dummies
        this->jumpTargetIndex = 9999;
        this->nextDiff = 9999;
        this->jumpDiff = 9999;

        //initialize the lock for this subspace
        omp_init_lock(&this->subspaceLock);
    }

    X86PerfectSubspaceNode(size_t dim, size_t index) {
        for (size_t i = 0; i < dim; i++) {
            this->level.push_back(1);
            this->hInverse.push_back(2);
        }
        this->gridPointsOnLevel = 0;
        this->existingGridPointsOnLevel = 0;
        this->flatLevel = 0;
        this->type = NOT_SET;

        //initalize other member variables with dummies
        this->jumpTargetIndex = index;
        this->nextDiff = dim;
        this->jumpDiff = dim;
        this->subspaceArray = nullptr;
    }

    ~X86PerfectSubspaceNode() {
        if (this->type == ARRAY) {
            delete[] this->subspaceArray;
        }
    }

    void lockSubspace() {
        omp_set_lock(&this->subspaceLock);
    }

    void unlockSubspace() {
        omp_unset_lock(&this->subspaceLock);
    }

    //increases number of grid points on the subspace
    void addGridPoint(DataVector &index) {
        size_t dim = index.getSize();
        for (size_t i = 0; i < dim; i++) {
            this->indices.push_back(index.get(i));
        }
        this->existingGridPointsOnLevel += 1;
    }

    void printLevel() {
        for (size_t i = 0; i < level.size(); i++) {
            if (i > 0) {
                cout << ", ";
            }
            cout << level[i];
        }
        cout << endl;
    }

    // unpack has to be called when the subspace is set up (except for surplus valus)
    // this method will decide how to best represent the subspace (list or array type)
    // and prepare the subspace for its representation
    void unpack() {
        double usageRatio = (double) this->existingGridPointsOnLevel / (double) this->gridPointsOnLevel;

        if (usageRatio < X86PERFECT_LIST_RATIO && this->existingGridPointsOnLevel < X86PERFECT_STREAMING_THRESHOLD) {
            this->type = LIST;
        } else {
            this->type = ARRAY;
            if (this->subspaceArray == nullptr) {
                // space for the coefficient/isLeaf-tuples
                this->subspaceArray = new double[this->gridPointsOnLevel * 2];
                for (size_t i = 0; i < this->gridPointsOnLevel * 2; i++) {
                    this->subspaceArray[i] = std::numeric_limits<double>::quiet_NaN();
                }
            }
        }
    }

    // the first call initializes the array for ARRAY type subspaces
    //
    void setSurplus(size_t indexFlat, double surplus, bool isLeaf) {
        if (this->type == ARRAY) {
            this->subspaceArray[indexFlat] = surplus;

            if (isLeaf) {
                this->subspaceArray[indexFlat + 1] = 1.0;
            } else {
                this->subspaceArray[indexFlat + 1] = 0.0;
            }

        } else if (this->type == LIST) {
            //TODO: only for debugging
            //this->subspaceArray[indexFlat] = surplus;
            bool found = false;

            //for (pair<uint32_t, double> &tuple : this->indexFlatSurplusPairs) {
            for (ListSubspaceTuple &tuple : this->indexFlatSurplusPairs) {
                //if (tuple.first == indexFlat) {
                if (tuple.indexFlat == indexFlat) {
                    //tuple.second = surplus;
                    tuple.surplus = surplus;
                    tuple.isLeaf = isLeaf;
                    found = true;
                    break;
                }
            }
            if (!found) {
                //this->indexFlatSurplusPairs.emplace_back(make_pair(indexFlat, surplus));
                this->indexFlatSurplusPairs.emplace_back(ListSubspaceTuple{indexFlat, surplus, isLeaf});
            }
        }
    }

    // the first call initializes the array for ARRAY type subspaces
    //
    double getSurplus(size_t indexFlat) {

        if (this->type == ARRAY) {
            return this->subspaceArray[indexFlat];
        } else if (this->type == LIST) {
            //for (tuple<uint32_t, double, double> myTuple : this->indexFlatSurplusPairs) {
            for (ListSubspaceTuple tuple : this->indexFlatSurplusPairs) {
                //if (tuple.first == indexFlat) {
                if (tuple.indexFlat == indexFlat) {
                    //return tuple.second;
                    return tuple.surplus;
                }
            }
        }
        throw;
    }

    static size_t compareLexicographically(X86PerfectSubspaceNode &current, X86PerfectSubspaceNode &last) {
        for (size_t i = 0; i < current.level.size(); i++) {
            if (current.level[i] != last.level[i]) {
                return i;
            }
        }
        throw "illegal input";
    }

    static bool subspaceCompare(X86PerfectSubspaceNode left, X86PerfectSubspaceNode right) {
        for (size_t i = 0; i < left.level.size(); i++) {
            if (left.level[i] >= right.level[i]) {
                if (left.level[i] > right.level[i]) {
                    return 0;
                }
            } else {
                return 1;
            }
        }
        return 1;
    }

};

}
}
