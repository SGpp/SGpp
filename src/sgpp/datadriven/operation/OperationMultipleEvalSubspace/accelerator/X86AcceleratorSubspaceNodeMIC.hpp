#pragma once

#include <vector>
#include <limits>

#include <omp.h>

#include <base/datatypes/DataVector.hpp>

using namespace std;
using namespace sg::base;

namespace sg {
namespace parallel {
class X86AcceleratorSubspaceNodeMIC {
public:

    //vector<uint32_t> level;
    //vector<uint32_t> hInverse;
    uint32_t hInverseOffset;

    uint32_t gridPointsOnLevel;
    uint32_t existingGridPointsOnLevel;


    uint32_t jumpTargetIndex;
    uint32_t flatLevel;

    // every node that reaches this subspace has to calculate this diff
    uint32_t arriveDiff;

    uint32_t surplusOffset;

};

}
}
