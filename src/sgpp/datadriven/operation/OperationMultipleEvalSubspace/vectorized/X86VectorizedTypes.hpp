#pragma once

#include <immintrin.h>

//calculates the array index of a subspace member variable
#define GET_SUBS(subspaceIndex, offset) subspaceIndex * this->subspaceSize + offset
#define GET_SUBS_ARRAY(subspaceIndex, offset, index) subspaceIndex * this->subspaceSize + offset + index
#define GET_NEXT(subspaceIndex) subspaceIndex * this->subspaceSize + this->nextOffset
#define GET_FLATLEVEL(subspaceIndex) subspaceIndex * this->subspaceSize + this->flatLevelOffset
#define GET_NEXTDIFF(subspaceIndex) subspaceIndex * this->subspaceSize + this->nextDiffOffset
#define GET_JUMPDIFF(subspaceIndex) subspaceIndex * this->subspaceSize + this->jumpDiffOffset
//#define GET_FLATLEVELPTR(subspaceIndex) subspaceIndex * this->subspaceSize + this->flatLevelPointerOffset

#define GET_FLATLEVELPTR_ON_PTR(subspaceIndex) reinterpret_cast<double **>(this->allSubspaces + subspaceIndex * this->subspaceSize + this->flatLevelPointerOffset)

#define GET_HINVERSE_LOOP(actualIndex) actualIndex + this->hInverseOffset
#define GET_NEXT_LOOP(actualIndex) actualIndex + this->nextOffset
#define GET_FLATLEVEL_LOOP(actualIndex) actualIndex + this->flatLevelOffset
#define GET_NEXTDIFF_LOOP(actualIndex) actualIndex + this->nextDiffOffset
#define GET_JUMPDIFF_LOOP(actualIndex) actualIndex + this->jumpDiffOffset

#define GET_LEVEL_PTR(subspaceIndex) this->allSubspaces + subspaceIndex * this->subspaceSize + this->levelOffset
#define GET_HINVERSE_PTR(subspaceIndex) this->allSubspaces + subspaceIndex * this->subspaceSize + this->hInverseOffset
#define GET_NEXT_PTR(subspaceIndex) this->allSubspaces + subspaceIndex * this->subspaceSize + this->nextOffset
#define GET_FLATLEVEL_PTR(subspaceIndex) this->allSubspaces + subspaceIndex * this->subspaceSize + this->flatLevelOffset
#define GET_NEXTDIFF_PTR(subspaceIndex) this->allSubspaces + subspaceIndex * this->subspaceSize + this->nextDiffOffset
#define GET_JUMPDIFF_PTR(subspaceIndex) this->allSubspaces + subspaceIndex * this->subspaceSize + this->jumpDiffOffset

