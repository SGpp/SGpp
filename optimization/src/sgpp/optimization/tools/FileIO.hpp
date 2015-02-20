// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#ifndef SGPP_OPTIMIZATION_TOOLS_FILEIO_HPP
#define SGPP_OPTIMIZATION_TOOLS_FILEIO_HPP

#include <sgpp/globaldef.hpp>

#include <fstream>

#include <sgpp/base/datatypes/DataVector.hpp>
#include <sgpp/base/datatypes/DataMatrix.hpp>
#include <sgpp/optimization/gridgen/IterativeGridGenerator.hpp>

namespace SGPP {
  namespace optimization {

    /**
     * Namespace with functions to write data (vectors, matrices, grids, ...)
     * to a file.
     */
    namespace file_io {
      /**
       * Write bytes of entry representation to output stream
       * (called by printMatrixToFile()).
       *
       * @param f         output stream to be written to
       * @param entry     entry to be written
       */
      template <class T>
      void writeEntryToFile(std::ofstream& f, const T& entry) {
        f.write(reinterpret_cast<const char*>(&entry), sizeof(T));
      }

      /**
       * Write string to output stream (called by printMatrixToFile()).
       *
       * @param f         output stream to be written to
       * @param entry     entry to be written
       */
      template <>
      void writeEntryToFile(std::ofstream& f, const std::string& entry) {
        // write string terminated by null character
        const char null_char[1] = {'\0'};
        f << entry;
        f.write(null_char, 1);
      }

      /**
       * @return  type string for unknown types (right-padded "other")
       */
      template <class T>
      const char* getTypeString(const std::vector<T>& A) {
        (void)A;
        return "other           ";
      }

      /**
       * @return  type string for uint8_t (right-padded "uint8")
       */
      template <>
      const char* getTypeString(const std::vector<uint8_t>& A) {
        (void)A;
        return "uint8           ";
      }

      /**
       * @return  type string for uint16_t (right-padded "uint16")
       */
      template <>
      const char* getTypeString(const std::vector<uint16_t>& A) {
        (void)A;
        return "uint16          ";
      }

      /**
       * @return  type string for uint32_t (right-padded "uint32")
       */
      template <>
      const char* getTypeString(const std::vector<uint32_t>& A) {
        (void)A;
        return "uint32          ";
      }

      /**
       * @return  type string for uint64_t (right-padded "uint64")
       */
      template <>
      const char* getTypeString(const std::vector<uint64_t>& A) {
        (void)A;
        return "uint64          ";
      }

      /**
       * @return  type string for float (right-padded "float")
       */
      template <>
      const char* getTypeString(const std::vector<float>& A) {
        (void)A;
        return "float           ";
      }

      /**
       * @return  type string for double (right-padded "double")
       */
      template <>
      const char* getTypeString(const std::vector<double>& A) {
        (void)A;
        return "double          ";
      }

      /**
       * @return  type string for string (right-padded "string")
       */
      template <>
      const char* getTypeString(const std::vector<std::string>& A) {
        (void)A;
        return "string          ";
      }

      /**
       * Write a grid (grid points and function values) to a file.
       *
       * The format is as follows:
       *
       * <pre>
       * size_t   N (number of grid points)
       * size_t   d (dimension)
       * for j = 0, ..., N-1
       *     for t = 0, ..., d-1
       *         float_t           grid_point[j].getCoord(t)
       *         unsigned int     grid_point[j].level(t)
       *         unsigned int     grid_point[j].index(t)
       *     end
       *     float_t   function_value[j]
       * end
       * </pre>
       *
       * @param filename      filename of the file to be written
       * @param gridGen       iterative grid generator containing the
       *                      grid points and function values
       */
      void writeGridToFile(const std::string& filename,
                           const IterativeGridGenerator& gridGen);

      /**
       * Write a base::DataMatrix to a file.
       *
       * @param filename      filename of the file to be written
       * @param A             matrix
       */
      void writeMatrixToFile(const std::string& filename, base::DataMatrix& A);

      /**
       * Write a matrix (stored row-wise in a std::vector) to a file.
       *
       * The format is as follows:
       *
       * <pre>
       * size_t       m
       * size_t       n
       * char[16]     type string (one of "uint8", "uint16", "uint32",
       *              "uint64", "float_t", "string", or "other",
       *              right-padded with spaces to 16 characters)
       * for i = 0, ..., m*n - 1
       *     T        A[i] (size depending on template parameter,
       *              strings are written null-terminatedly)
       * end
       * </pre>
       *
       * @param filename      filename of the file to be written
       * @param A             matrix
       * @param m             number of rows
       * @param n             number of columns
       */
      template <class T>
      void writeMatrixToFile(const std::string& filename,
                             const std::vector<T>& A,
                             size_t m, size_t n) {
        std::ofstream f;
        f.exceptions(std::ofstream::failbit | std::ofstream::badbit);
        const char* type = getTypeString(A);

        f.open(filename.c_str(), std::ios::binary);

        // header (size and type)
        f.write(reinterpret_cast<const char*>(&m), sizeof(m));
        f.write(reinterpret_cast<const char*>(&n), sizeof(n));
        f.write(type, 16);

        // entries
        for (size_t i = 0; i < m * n; i++) {
          writeEntryToFile(f, A[i]);
        }

        f.close();
      }

      /**
       * Write a matrix (stored as a vector of row vectors) to a file.
       * Copies the matrix and calls the other version with the
       * vectorized matrix.
       *
       * @param filename      filename of the file to be written
       * @param A             matrix
       */
      template <class T>
      void writeMatrixToFile(const std::string& filename,
                             const std::vector<std::vector<T>>& A) {
        size_t m = A.size();
        size_t n = (A.empty() ? 0 : A[0].size());
        std::vector<T> B(m * n);

        for (size_t i = 0; i < m; i++) {
          std::copy(A[i].begin(), A[i].end(), B.begin() + i * n);
        }

        writeMatrixToFile(filename, B, m, n);
      }

      /**
       * Write a base::DataVector to a file.
       * It's printMatrixToFile with the vector as one column.
       *
       * @param filename      filename of the file to be written
       * @param x             vector
       */
      void writeVectorToFile(const std::string& filename, base::DataVector& x);

      /**
       * Write a std::vector to a file.
       * It's printMatrixToFile with the vector as one column.
       *
       * @param filename      filename of the file to be written
       * @param x             vector
       */
      template <class T>
      void writeVectorToFile(const std::string& filename,
                             const std::vector<T>& x) {
        writeMatrixToFile(filename, x, x.size(), 1);
      }

    }
  }
}

#endif /* SGPP_OPTIMIZATION_TOOLS_FILEIO_HPP */
