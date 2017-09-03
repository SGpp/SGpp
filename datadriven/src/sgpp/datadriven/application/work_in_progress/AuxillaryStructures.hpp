//
// Created by Vincent_Bode on 03.09.2017.
//

#ifndef SGPP_AUXILLARYSTRUCTURES_H
#define SGPP_AUXILLARYSTRUCTURES_H

#include <sgpp_datadriven.hpp>

namespace sgpp {
    namespace datadriven {
        typedef ClassDensityConntainer ClassDensityContainer;

        struct LevelIndexPair {
            unsigned long level;
            unsigned long index;
        };

        typedef std::vector<LevelIndexPair> LevelIndexVector;

        struct RefinementResult {
            std::list<LevelIndexVector> addedGridPoints;
            std::list<size_t> deletedGridPointsIndexes;
        };

    }
}

#endif //SGPP_AUXILLARYSTRUCTURES_H
