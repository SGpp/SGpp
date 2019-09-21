// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <sgpp/base/algorithm/AlgorithmDGEMV.hpp>

#include <sgpp/base/operation/hash/OperationMultipleEvalInterModLinear.hpp>
#include <sgpp/base/operation/hash/common/basis/LinearModifiedBasis.hpp>

#include <sgpp/globaldef.hpp>

#include <algorithm>
#include <set>
#include <vector>

namespace sgpp {
namespace base {

void OperationMultipleEvalInterModLinear::mult(DataVector& alpha, DataVector& result) {
  /*
  AlgorithmDGEMV<SLinearModifiedBase> op;
  LinearModifiedBasis<unsigned int, unsigned int> base;

  op.mult(storage, base, alpha, this->dataset, result);
  */
  result.setAll(0.0);


  #pragma omp parallel
  {
    DataVector line(dataset.getNcols());
    DataVector privateResult(result.getSize());
    privateResult.setAll(0.0);

    GridStorage::grid_iterator working(storage);
    LinearModifiedBasis<unsigned int, unsigned int> basis;
    size_t dimensions = storage.getDimension();

    #pragma omp for
    for (size_t j = 0; j < dataset.getNrows(); j++) {
      dataset.getRow(j, line);
      // iterate over all interactions
      for (const std::set<size_t>& in : interactions) {
        bool pointComputed = true;
        bool lvlComplete = false;
        size_t relLvl = 0;
        level_t* level = new level_t[in.size()];
        for (size_t i = 0; i < in.size(); i++) {
          level[i] = 0;
        }
        // only view affected subspaces
        do {
          pointComputed = (!lvlComplete) & pointComputed;
          lvlComplete = false;
          for (size_t d = 0; d < dimensions; d++) {
            working.push(d, 1, 1);
          }

          size_t index = 0;
          for (size_t i : in) {
            level_t lvl = static_cast<level_t> (2 + level[index]);
            index_t idx = static_cast<index_t>(std::min(1+2*
              floor(line[i] * static_cast<double>(1 << (lvl-1))),
              static_cast<double>(1 << lvl)-1.));
            // only hash for last element
            if (index++ == in.size()-1)
              working.set(i, lvl, idx);
            else
              working.push(i, lvl, idx);
          }

          // hash first dimension if there was no element in the interaction
          if (in.size() == 0) working.set(0, 1, 1);

          size_t seq = working.seq();

          if (!storage.isInvalidSequenceNumber(seq)) {
            double value = alpha[seq];
            index_t work_index;
            level_t work_level;
            for (size_t i : in) {
              working.get(i, work_level, work_index);
              value *= basis.eval(work_level, work_index, line[i]);
            }
            privateResult[j] += value;
            pointComputed = true;
          }

          if (in.size() == 0) break;
          // update subspace lvl
          if (level[in.size()-1] == relLvl) {
            relLvl++;
            lvlComplete = true;
            for (size_t i = 0; i < in.size(); i++) {
              level[i] = 0;
            }
          }
          size_t levelsum;
          // find next subspace for currentlvl
          do {
            levelsum = 0;
            bool carry = true;
            for (size_t i = 0; i < in.size(); i++) {
              if (carry) level[i]++;
              if (level[i] == relLvl+1) {
                level[i] = 0;
                carry = true;
              } else {
                carry = false;
              }
            }
            for (size_t i = 0; i < in.size(); i++) {
              levelsum += level[i];
            }
          } while (levelsum != relLvl);
        }while((pointComputed||!lvlComplete));
        delete[] level;
      }
    }
    #pragma omp critical
    {
      result.add(privateResult);
    }
  }
}

void OperationMultipleEvalInterModLinear::multTranspose(DataVector& source, DataVector& result) {
  result.setAll(0.0);

  #pragma omp parallel
  {
    DataVector line(dataset.getNcols());
    DataVector privateResult(result.getSize());
    privateResult.setAll(0.0);

    GridStorage::grid_iterator working(storage);
    LinearModifiedBasis<unsigned int, unsigned int> basis;

    size_t dimensions = storage.getDimension();

    #pragma omp for
    for (size_t j = 0; j < dataset.getNrows(); j++) {
      dataset.getRow(j, line);
      // iterate over all interactions
      for (const std::set<size_t>& in : interactions) {
        bool pointComputed = true;
        bool lvlComplete = false;
        size_t relLvl = 0;
        level_t* level = new level_t[in.size()];
        for (size_t i = 0; i < in.size(); i++) {
          level[i] = 0;
        }
        // only view affected subspaces
        do {
          pointComputed = (!lvlComplete) & pointComputed;
          lvlComplete = false;
          for (size_t d = 0; d < dimensions; d++) {
            working.push(d, 1, 1);
          }

          size_t index = 0;
          for (size_t i : in) {
            level_t lvl = static_cast<level_t> (2 + level[index]);
            index_t idx = static_cast<index_t>(std::min(1+2*
              floor(line[i] * static_cast<double>(1 << (lvl-1))),
              static_cast<double>(1 << lvl)-1.));
            // only hash for last element
            if (index++ == in.size()-1) working.set(i, lvl, idx);
            else
              working.push(i, lvl, idx);
          }

          // hash first dimension if there was no element in the interaction
          if (in.size() == 0) working.set(0, 1, 1);

          size_t seq = working.seq();

          if (!storage.isInvalidSequenceNumber(seq)) {
            double value = source[j];
            index_t work_index;
            level_t work_level;
            for (size_t i : in) {
              working.get(i, work_level, work_index);
              value *= basis.eval(work_level, work_index, line[i]);
            }
            privateResult[seq] += value;
            pointComputed = true;
          }

          if (in.size() == 0) break;
          // update subspace lvl
          if (level[in.size()-1] == relLvl) {
            relLvl++;
            lvlComplete = true;
            for (size_t i = 0; i < in.size(); i++) {
              level[i] = 0;
            }
          }
          size_t levelsum;
          // find next subspace for currentlvl
          do {
            levelsum = 0;
            bool carry = true;
            for (size_t i = 0; i < in.size(); i++) {
              if (carry) level[i]++;
              if (level[i] == relLvl+1) {
                level[i] = 0;
                carry = true;
              } else {
                carry = false;
              }
            }
            for (size_t i = 0; i < in.size(); i++) {
              levelsum += level[i];
            }
          } while (levelsum != relLvl);
        }while((pointComputed||!lvlComplete));
        delete[] level;
      }
    }
    #pragma omp critical
    {
      result.add(privateResult);
    }
  }
}


}  // namespace base
}  // namespace sgpp
