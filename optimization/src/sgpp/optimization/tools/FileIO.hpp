// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#ifndef SGPP_OPTIMIZATION_TOOLS_FILEIO_HPP
#define SGPP_OPTIMIZATION_TOOLS_FILEIO_HPP

#include <sgpp/globaldef.hpp>

#include <sgpp/base/datatypes/DataVector.hpp>
#include <sgpp/base/datatypes/DataMatrix.hpp>
#include <sgpp/optimization/gridgen/IterativeGridGenerator.hpp>

#include <fstream>
#include <stdexcept>
#include <string>
#include <vector>

namespace sgpp {
namespace optimization {

/**
 * Namespace with functions to write data (vectors, matrices, grids, ...)
 * to a file.
 */
namespace file_io {
/**
 * Write bytes of entry representation to output stream
 * (called by writeMatrix()).
 *
 * @param f         output stream to be written to
 * @param entry     entry to be written
 */
template <class T>
void writeEntry(std::ofstream& f, const T& entry) {
  f.write(reinterpret_cast<const char*>(&entry), sizeof(T));
}

/**
 * Write string to output stream (called by writeMatrix()).
 *
 * @param f         output stream to be written to
 * @param entry     entry to be written
 */
template <>
void writeEntry(std::ofstream& f, const std::string& entry);

/**
 * Read bytes of entry representation from input stream
 * (called by readMatrix()).
 *
 * @param f           input stream to be read
 * @param[out] entry  entry to be read
 */
template <class T>
void readEntry(std::ifstream& f, T& entry) {
  f.read(reinterpret_cast<char*>(&entry), sizeof(T));
}

/**
 * Read string from input stream (called by readMatrix()).
 *
 * @param f           input stream to be read
 * @param[out] entry  entry to be read
 */
template <>
void readEntry(std::ifstream& f, std::string& entry);

/**
 * @param A ignored
 * @return  type string for unknown types (right-padded "other")
 */
template <class T>
const char* getTypeString(const std::vector<T>& A) {
  (void)A;
  return "other           ";
}

/**
 * @param A ignored
 * @return  type string for uint8_t (right-padded "uint8")
 */
template <>
const char* getTypeString(const std::vector<uint8_t>& A);

/**
 * @param A ignored
 * @return  type string for uint16_t (right-padded "uint16")
 */
template <>
const char* getTypeString(const std::vector<uint16_t>& A);

/**
 * @param A ignored
 * @return  type string for uint32_t (right-padded "uint32")
 */
template <>
const char* getTypeString(const std::vector<uint32_t>& A);

/**
 * @param A ignored
 * @return  type string for uint64_t (right-padded "uint64")
 */
template <>
const char* getTypeString(const std::vector<uint64_t>& A);

/**
 * @param A ignored
 * @return  type string for float (right-padded "float")
 */
template <>
const char* getTypeString(const std::vector<float>& A);

/**
 * @param A ignored
 * @return  type string for double (right-padded "double")
 */
template <>
const char* getTypeString(const std::vector<double>& A);

/**
 * @param A ignored
 * @return  type string for string (right-padded "string")
 */
template <>
const char* getTypeString(const std::vector<std::string>& A);

/**
 * Write a grid (only grid points) to a file.
 * The format is the same as the version with functions values
 * with all function values set to zero.
 *
 * @param filename        filename of the file to be written
 * @param gridStorage     grid storage containing the grid points
 */
void writeGrid(const std::string& filename, const base::GridStorage& gridStorage);

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
 *         double          grid_point[j].getCoord(t)
 *         unsigned int     grid_point[j].level(t)
 *         unsigned int     grid_point[j].index(t)
 *     end
 *     double   function_value[j]
 * end
 * </pre>
 *
 * @param filename        filename of the file to be written
 * @param gridStorage     grid storage containing the grid points
 * @param functionValues  vector of function values
 */
void writeGrid(const std::string& filename, const base::GridStorage& gridStorage,
               const base::DataVector& functionValues);

/**
 * Read a grid (only grid points) from a file.
 * The format is as in writeGrid (discarding function values).
 *
 * @param filename              filename of the file to be read
 * @param[out] gridStorage      grid storage containing the grid points
 */
void readGrid(const std::string& filename, base::GridStorage& gridStorage);

/**
 * Read a grid (grid points and function values) from a file.
 * The format is as in writeGrid.
 *
 * @param filename              filename of the file to be read
 * @param[out] gridStorage      grid storage containing the grid points
 * @param[out] functionValues   vector of function values
 */
void readGrid(const std::string& filename, base::GridStorage& gridStorage,
              base::DataVector& functionValues);

/**
 * Write a base::DataMatrix to a file.
 *
 * @param filename      filename of the file to be written
 * @param A             matrix
 */
void writeMatrix(const std::string& filename, base::DataMatrix& A);

/**
 * Write a matrix (stored row-wise in a std::vector) to a file.
 *
 * The format is as follows:
 *
 * <pre>
 * size_t       m
 * size_t       n
 * char[16]     type string (one of "uint8", "uint16", "uint32",
 *              "uint64", "double", "string", or "other",
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
void writeMatrix(const std::string& filename, const std::vector<T>& A, size_t m, size_t n) {
  std::ofstream f;
  f.exceptions(std::ofstream::failbit | std::ofstream::badbit);
  const char* type = getTypeString(A);

  f.open(filename.c_str(), std::ios::out | std::ios::binary);

  // header (size and type)
  f.write(reinterpret_cast<const char*>(&m), sizeof(m));
  f.write(reinterpret_cast<const char*>(&n), sizeof(n));
  f.write(type, 16);

  // entries
  for (size_t i = 0; i < m * n; i++) {
    writeEntry(f, A[i]);
  }
}

/**
 * Read a matrix from a file.
 * The format is as in writeMatrix.
 *
 * @param filename      filename of the file to be read
 * @param[out] A        matrix
 */
void readMatrix(const std::string& filename, base::DataMatrix& A);

/**
 * Read a matrix (stored row-wise in a std::vector) from a file.
 * The format is as in writeMatrix.
 *
 * @param filename      filename of the file to be written
 * @param[out] A        matrix
 * @param[out] m        number of rows
 * @param[out] n        number of columns
 */
template <class T>
void readMatrix(const std::string& filename, std::vector<T>& A, size_t& m, size_t& n) {
  std::ifstream f;
  f.exceptions(std::ofstream::failbit | std::ofstream::badbit);
  const char* type = getTypeString(A);

  f.open(filename.c_str(), std::ios::in | std::ios::binary);

  // header (size and type)
  f.read(reinterpret_cast<char*>(&m), sizeof(m));
  f.read(reinterpret_cast<char*>(&n), sizeof(n));

  char type2[17];
  f.read(type2, 16);
  type2[16] = 0;

  if (std::string(type) != std::string(type2)) {
    throw std::invalid_argument(
        "The type of the entries in the file "
        "differ from the type of the "
        "entries of A.");
  }

  A.empty();

  // entries
  for (size_t i = 0; i < m * n; i++) {
    T entry;
    readEntry(f, entry);
    A.push_back(entry);
  }
}

/**
 * Write a base::DataVector to a file.
 * It's writeMatrix with the vector as one row.
 *
 * @param filename      filename of the file to be written
 * @param x             vector
 */
void writeVector(const std::string& filename, base::DataVector& x);

/**
 * Write a std::vector to a file.
 * It's writeMatrix with the vector as one row.
 *
 * @param filename      filename of the file to be written
 * @param x             vector
 */
template <class T>
void writeVector(const std::string& filename, const std::vector<T>& x) {
  writeMatrix(filename, x, 1, x.size());
}

/**
 * Read a base::DataVector from a file.
 * It's readMatrix with the vector as one row.
 *
 * @param filename      filename of the file to be read
 * @param[out] x        vector
 */
void readVector(const std::string& filename, base::DataVector& x);

/**
 * Read a std::vector from a file.
 * It's readMatrix with the vector as one row.
 *
 * @param filename      filename of the file to be read
 * @param[out] x        vector
 */
template <class T>
void readVector(const std::string& filename, std::vector<T>& x) {
  size_t m, n;
  readMatrix(filename, x, m, n);
}
}  // namespace file_io
}  // namespace optimization
}  // namespace sgpp

#endif /* SGPP_OPTIMIZATION_TOOLS_FILEIO_HPP */
