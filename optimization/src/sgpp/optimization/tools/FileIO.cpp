// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <sgpp/globaldef.hpp>

#include <sgpp/optimization/tools/FileIO.hpp>
#include <sgpp/base/grid/GridStorage.hpp>

#include <iostream>
#include <stdexcept>

namespace SGPP {
  namespace optimization {
    namespace file_io {

      void writeGrid(const std::string& filename,
                     const base::GridStorage& gridStorage) {
        const size_t N = gridStorage.size();
        const std::vector<float_t> functionValues(N, 0.0);
        writeGrid(filename, gridStorage, functionValues);
      }

      void writeGrid(const std::string& filename,
                     const base::GridStorage& gridStorage,
                     const std::vector<float_t>& functionValues) {
        const size_t N = gridStorage.size();
        const size_t d = gridStorage.dim();

        if (functionValues.size() != N) {
          throw std::invalid_argument("functionValues must have as many "
                                      "elements as there are grid "
                                      "points in gridStorage.");
        }

        std::ofstream f;
        f.exceptions(std::ofstream::failbit | std::ofstream::badbit);
        f.open(filename.c_str(), std::ios::out | std::ios::binary);

        // header (number of points and dimensions)
        f.write(reinterpret_cast<const char*>(&N), sizeof(N));
        f.write(reinterpret_cast<const char*>(&d), sizeof(d));

        for (size_t j = 0; j < N; j++) {
          base::GridIndex& gp = *gridStorage.get(j);

          for (size_t t = 0; t < d; t++) {
            float_t x = gp.getCoord(t);
            base::GridIndex::level_type l = gp.getLevel(t);
            base::GridIndex::index_type i = gp.getIndex(t);

            // coordinate, level and index of current grid point
            f.write(reinterpret_cast<const char*>(&x), sizeof(x));
            f.write(reinterpret_cast<const char*>(&l), sizeof(l));
            f.write(reinterpret_cast<const char*>(&i), sizeof(i));
          }

          // function value at the current grid point
          f.write(reinterpret_cast<const char*>(&functionValues[j]),
                  sizeof(float_t));
        }
      }

      void readGrid(const std::string& filename,
                    base::GridStorage& gridStorage) {
        std::vector<float_t> functionValues;
        readGrid(filename, gridStorage, functionValues);
      }

      void readGrid(const std::string& filename,
                    base::GridStorage& gridStorage,
                    std::vector<float_t>& functionValues) {
        std::ifstream f;
        f.exceptions(std::ifstream::failbit | std::ifstream::badbit);
        f.open(filename.c_str(), std::ios::in | std::ios::binary);

        size_t N, d;
        f.read(reinterpret_cast<char*>(&N), sizeof(N));
        f.read(reinterpret_cast<char*>(&d), sizeof(d));

        gridStorage.emptyStorage();
        functionValues.clear();

        base::GridIndex gp(d);
        float_t x;
        base::GridIndex::level_type l;
        base::GridIndex::index_type i;
        float_t functionValue;

        for (size_t j = 0; j < N; j++) {
          for (size_t t = 0; t < d; t++) {
            f.read(reinterpret_cast<char*>(&x), sizeof(x));
            f.read(reinterpret_cast<char*>(&l), sizeof(l));
            f.read(reinterpret_cast<char*>(&i), sizeof(i));
            gp.set(t, l, i);
          }

          f.read(reinterpret_cast<char*>(&functionValue), sizeof(functionValue));
          gridStorage.insert(gp);
          functionValues.push_back(functionValue);
        }
      }

      void writeMatrix(const std::string& filename,
                       base::DataMatrix& A) {
        // convert DataMatrix to std::vector
        const std::vector<float_t> AVector(A.getPointer(), A.getPointer() +
                                           A.getNrows() * A.getNcols());
        writeMatrix(filename, AVector, A.getNrows(), A.getNcols());
      }

      void readMatrix(const std::string& filename, base::DataMatrix& A) {
        std::vector<float_t> AVector;
        size_t m, n;
        readMatrix(filename, AVector, m, n);
        A.resize(m, n);
        A = base::DataMatrix(&AVector[0], m, n);
      }

      void writeVector(const std::string& filename,
                       base::DataVector& x) {
        // convert DataVector to std::vector
        const std::vector<float_t> xVector(x.getPointer(),
                                           x.getPointer() + x.getSize());
        writeMatrix(filename, xVector, 1, x.getSize());
      }

      void readVector(const std::string& filename, base::DataVector& x) {
        std::vector<float_t> xVector;
        size_t m, n;
        readMatrix(filename, xVector, m, n);
        x.resize(xVector.size());
        x = base::DataVector(&xVector[0], xVector.size());
      }

    }
  }
}
