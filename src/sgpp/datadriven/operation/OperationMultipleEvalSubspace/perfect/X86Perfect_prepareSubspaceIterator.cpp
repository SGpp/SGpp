#include "../../OperationMultipleEvalSubspace/perfect/X86Perfect.hpp"

/*
 levels needs to have size of the number of subspaces, stores a level tuple for each subspace
 level level part of a level index tuple, needs to have size of grid
 index index part of a level index tuple, needs to have size of grid

 Important: In contrast to other eval preparation method, this creates the actual level, not 2**level
 */
void X86Perfect::prepareSubspaceIterator() {

    /////////////////////////////////////////////////////
    // extract the subspace and grid points
    // and put them in a map with flatLevel -> subspace
    /////////////////////////////////////////////////////

    DataVector level(this->dim);
    DataVector index(this->dim);
    DataVector maxIndex(this->dim);

    unsigned int curLevel;
    unsigned int curIndex;

    this->allSubspaceNodes.clear();
    this->subspaceCount = 0;
    this->maxLevel = 0;

    //calculate the maxLevel first - required for level vector flattening
    for (size_t gridIndex = 0; gridIndex < this->storage->size(); gridIndex++) {
        sg::base::GridIndex *point = this->storage->get(gridIndex);
        for (size_t d = 0; d < this->dim; d++) {
            point->get(d, curLevel, curIndex);
            if (curLevel > this->maxLevel) {
                this->maxLevel = curLevel;
            }
        }
    }

    //create a list of subspaces and setup the grid points - now we know which subspaces actually exist in the grid
    //create map of flatLevel -> subspaceIndex (for convenience operations, not for time-critical code)
    //also find out how many grid points the largest subspace contains
    this->maxGridPointsOnLevel = 0;
    for (size_t gridIndex = 0; gridIndex < this->storage->size(); gridIndex++) {
        sg::base::GridIndex *point = this->storage->get(gridIndex);

        for (size_t d = 0; d < this->dim; d++) {
            point->get(d, curLevel, curIndex);
            level.set(d, curLevel);
            index.set(d, curIndex);
            maxIndex.set(d, 1 << curLevel);
        }

        uint32_t flatLevel = X86Perfect::flattenLevel(this->dim, maxLevel, level);

        map<uint32_t, uint32_t>::iterator it = this->allLevelsIndexMap.find(flatLevel);
        if (it == this->allLevelsIndexMap.end()) {
            this->allLevelsIndexMap.insert(make_pair(flatLevel, this->subspaceCount));

            this->allSubspaceNodes.emplace_back(X86PerfectSubspaceNode(level, flatLevel, maxIndex, index));

            X86PerfectSubspaceNode &subspace = this->allSubspaceNodes[this->subspaceCount];
            if (subspace.gridPointsOnLevel > this->maxGridPointsOnLevel) {
                this->maxGridPointsOnLevel = subspace.gridPointsOnLevel;
            }
            this->subspaceCount += 1;
        } else {
            // add the current grid point to its subspace
            uint32_t subspaceIndex = it->second;
            X86PerfectSubspaceNode &subspace = this->allSubspaceNodes[subspaceIndex];
            subspace.addGridPoint(index);
        }
    }

    // sort the subspaces lexicographically to get efficient curve though the subspaces
    std::sort(this->allSubspaceNodes.begin(), this->allSubspaceNodes.end(), X86PerfectSubspaceNode::subspaceCompare);

    // - after sorting the allLevelsIndexMap indices have to be recomputed
    // - the grid points are "unpacked", they choose their representation depending on their configuration
    // - also collect some statistical information
    this->allLevelsIndexMap.clear();
    size_t totalRegularGridPoints = 0;
    size_t actualGridPoints = 0;
    size_t largestArraySubspace = 0;
    size_t largestListSubspace = 0;
    size_t numberOfListSubspaces = 0;
    size_t nonVirtualGridPoints = 0;
    for (size_t subspaceIndex = 0; subspaceIndex < this->subspaceCount; subspaceIndex++) {
        X86PerfectSubspaceNode &subspace = this->allSubspaceNodes[subspaceIndex];
        //select representation
        subspace.unpack();

        this->allLevelsIndexMap.insert(make_pair(subspace.flatLevel, subspaceIndex));

        //collect statistics
        totalRegularGridPoints += subspace.gridPointsOnLevel;
        nonVirtualGridPoints += subspace.existingGridPointsOnLevel;
        if (subspace.type == X86PerfectSubspaceNode::SubspaceType::LIST) {
            actualGridPoints += subspace.existingGridPointsOnLevel;
            numberOfListSubspaces += 1;
            if (subspace.existingGridPointsOnLevel > largestListSubspace) {
                largestListSubspace = subspace.existingGridPointsOnLevel;
            }
        } else {
            actualGridPoints += subspace.gridPointsOnLevel;
            if (subspace.gridPointsOnLevel > largestArraySubspace) {
                largestArraySubspace = subspace.gridPointsOnLevel;
            }
        }

    }

    //cout << "maxGridPointsOnLevel: " << this->maxGridPointsOnLevel << endl;
    //cout << "no. of subspaces: " << subspaceCount << endl;
    //cout << "maxLevel: " << maxLevel << endl;

#ifdef X86PERFECT_WRITE_STATS
    // cout << "nonVirtualGridPoints: " << nonVirtualGridPoints << endl;
    // cout << "totalRegularGridPoints: " << totalRegularGridPoints << endl;
    // cout << "actualGridPoints: " << actualGridPoints << endl;
    // cout << "largestArraySubspace: " << largestArraySubspace << endl;
    // cout << "largestListSubspace: " << largestListSubspace << endl;
    // cout << "numberOfListSubspaces: " << numberOfListSubspaces << endl;
    // cout << "subspaceCount: " << subspaceCount << endl;
    // cout << "avr. points per subspace: " << ((double) nonVirtualGridPoints / (double) subspaceCount) << endl;

    this->statsFile << this->refinementStep << this->csvSep;
    this->statsFile << nonVirtualGridPoints << this->csvSep;
    this->statsFile << totalRegularGridPoints << this->csvSep;
    this->statsFile << actualGridPoints << this->csvSep;
    this->statsFile << largestArraySubspace << this->csvSep;
    this->statsFile << largestListSubspace << this->csvSep;
    this->statsFile << numberOfListSubspaces << this->csvSep;
    this->statsFile << subspaceCount << this->csvSep;
    this->statsFile << ((double) nonVirtualGridPoints / (double) subspaceCount) << this->csvSep;
    size_t numberOfThreads = omp_get_max_threads();
    this->statsFile << ((double) (this->maxGridPointsOnLevel * numberOfThreads + actualGridPoints) * 8.0) / (1024.0 * 1024.0) << this->csvSep;
    this->statsFile << (double) (this->maxGridPointsOnLevel * numberOfThreads + actualGridPoints) / nonVirtualGridPoints << endl;
    this->refinementStep += 1;
#endif

    //add padding subspace at the end
    this->allSubspaceNodes.emplace_back(X86PerfectSubspaceNode(this->dim, subspaceCount));

    //////////////////////////////////////
    // now the jump links can be created
    //////////////////////////////////////

    //marker for the padding subspace - traversal is finished when this subspace is reached
    uint32_t computationFinishedMarker = subspaceCount;

    //tracks at which position the next subspace can be found that has a higher value in at least the n-th component
    uint32_t *jumpIndexMap = new uint32_t[this->dim];
    for (size_t i = 0; i < this->dim; i++) {
        jumpIndexMap[i] = computationFinishedMarker;
    }

    for (size_t i = subspaceCount - 1; i > 0; i--) {
        X86PerfectSubspaceNode &currentNode = this->allSubspaceNodes[i];
        X86PerfectSubspaceNode &previousNode = this->allSubspaceNodes[i - 1];

        //number of dimension for which this subspace is responsible
        uint32_t recomputeComponentPreviousNode = X86PerfectSubspaceNode::compareLexicographically(currentNode,
                previousNode);

        // diff to the next element, if no jump is taken, required to get (near) O(1) index calculation
        uint32_t nextDiff = this->dim; //will result in doing nothing for the padding
        // diff to the jump destination (if the jump should be taken), required to get (near) O(1) index calculation
        uint32_t jumpDiff = this->dim;

        if (i != subspaceCount - 1) {
            X86PerfectSubspaceNode &nextNode = this->allSubspaceNodes[i + 1];
            nextDiff = X86PerfectSubspaceNode::compareLexicographically(currentNode, nextNode);
            uint32_t jumpDestination = jumpIndexMap[recomputeComponentPreviousNode];
            if (jumpDestination != computationFinishedMarker) {
                X86PerfectSubspaceNode &jumpNode = this->allSubspaceNodes[jumpDestination];
                jumpDiff = X86PerfectSubspaceNode::compareLexicographically(currentNode, jumpNode);
            }
        }

        //which dims are to be recomputed in the next processed (!) subspace (successor or jump)
        currentNode.nextDiff = nextDiff;
        currentNode.jumpDiff = jumpDiff;

        //where to jump to?
        currentNode.jumpTargetIndex = jumpIndexMap[recomputeComponentPreviousNode];

        //update jump map
        //current node is a valid jump target for all predeccessors that change in one of the higher dims
        // sssssssssss <change> uuuuuuuu (s = same component, u = different component)
        // for predecessors that changed in the "u" array, jump to the current node
        // they also update the jump map
        for (size_t j = recomputeComponentPreviousNode + 1; j < dim; j++) {
            jumpIndexMap[j] = i;
        }

    }
    delete jumpIndexMap;

    //make sure that the tensor products are calculated entirely at the first subspace
    X86PerfectSubspaceNode &firstNode = this->allSubspaceNodes[0];
    firstNode.jumpTargetIndex = computationFinishedMarker;
    firstNode.nextDiff = 0;
    firstNode.jumpDiff = 0;
    if (subspaceCount > 1) {
        X86PerfectSubspaceNode &secondNode = this->allSubspaceNodes[1];
        firstNode.nextDiff = X86PerfectSubspaceNode::compareLexicographically(firstNode, secondNode);
    }
}
