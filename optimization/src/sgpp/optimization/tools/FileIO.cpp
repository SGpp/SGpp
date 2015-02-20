// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <sgpp/globaldef.hpp>

#include <sgpp/optimization/tools/FileIO.hpp>
#include <sgpp/base/grid/GridStorage.hpp>

#include <iostream>

namespace SGPP {
  namespace optimization {
    namespace file_io {

      void writeGridToFile(const std::string& filename,
                           const IterativeGridGenerator& gridGen) {
        base::GridStorage& grid_storage = *gridGen.getGrid().getStorage();
        const std::vector<float_t>& function_values =
          gridGen.getFunctionValues();
        size_t N = grid_storage.size();
        size_t d = grid_storage.dim();

        std::ofstream f(filename.c_str(), std::ios::out | std::ios::binary);

        // header (number of points and dimensions)
        f.write(reinterpret_cast<const char*>(&N), sizeof(N));
        f.write(reinterpret_cast<const char*>(&d), sizeof(d));

        for (size_t j = 0; j < N; j++) {
          base::GridIndex& gp = *grid_storage.get(j);

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
          f.write(reinterpret_cast<const char*>(&function_values[j]),
                  sizeof(float_t));
        }

        f.close();
      }

      void writeMatrixToFile(const std::string& filename,
                             base::DataMatrix& A) {
        // convert DataMatrix to std::vector
        std::vector<float_t> AVector(A.getPointer(), A.getPointer() +
                                     A.getNrows() * A.getNcols());
        writeMatrixToFile(filename, AVector, A.getNrows(), A.getNcols());
      }

      void writeVectorToFile(const std::string& filename,
                             base::DataVector& x) {
        // convert DataVector to std::vector
        std::vector<float_t> xVector(x.getPointer(),
                                     x.getPointer() + x.getSize());
        writeMatrixToFile(filename, xVector, x.getSize(), 1);
      }

    }
  }
}
