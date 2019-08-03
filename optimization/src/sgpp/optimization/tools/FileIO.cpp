// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <sgpp/globaldef.hpp>

#include <sgpp/base/grid/GridStorage.hpp>
#include <sgpp/optimization/tools/FileIO.hpp>

#include <stdexcept>
#include <string>
#include <vector>  // namespace

namespace sgpp {
namespace optimization {
namespace file_io {

template <>
void writeEntry(std::ofstream& f, const std::string& entry) {
  // write string terminated by null character
  f << entry << '\0';
}

template <>
void readEntry(std::ifstream& f, std::string& entry) {
  char ch;
  entry = "";

  // read until null character
  while (f.get(ch)) {
    if (ch != '\0') {
      entry += ch;
    } else {
      break;
    }
  }
}

template <>
const char* getTypeString(const std::vector<uint8_t>& A) {
  (void)A;
  return "uint8           ";
}

template <>
const char* getTypeString(const std::vector<uint16_t>& A) {
  (void)A;
  return "uint16          ";
}

template <>
const char* getTypeString(const std::vector<uint32_t>& A) {
  (void)A;
  return "uint32          ";
}

template <>
const char* getTypeString(const std::vector<uint64_t>& A) {
  (void)A;
  return "uint64          ";
}

template <>
const char* getTypeString(const std::vector<float>& A) {
  (void)A;
  return "float           ";
}

template <>
const char* getTypeString(const std::vector<double>& A) {
  (void)A;
  return "double          ";
}

template <>
const char* getTypeString(const std::vector<std::string>& A) {
  (void)A;
  return "string          ";
}

void writeGrid(const std::string& filename, const base::GridStorage& gridStorage) {
  const size_t N = gridStorage.getSize();
  const base::DataVector functionValues(N, 0.0);
  writeGrid(filename, gridStorage, functionValues);
}

void writeGrid(const std::string& filename, const base::GridStorage& gridStorage,
               const base::DataVector& functionValues) {
  const size_t N = gridStorage.getSize();
  const size_t d = gridStorage.getDimension();

  if (functionValues.getSize() != N) {
    throw std::invalid_argument(
        "functionValues must have as many "
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
    const base::GridPoint& gp = gridStorage[j];

    for (size_t t = 0; t < d; t++) {
      const double x = gridStorage.getCoordinate(gp, t);
      base::GridPoint::level_type l = gp.getLevel(t);
      base::GridPoint::index_type i = gp.getIndex(t);

      // coordinate, level and index of current grid point
      f.write(reinterpret_cast<const char*>(&x), sizeof(x));
      f.write(reinterpret_cast<const char*>(&l), sizeof(l));
      f.write(reinterpret_cast<const char*>(&i), sizeof(i));
    }

    // function value at the current grid point
    const double fX = functionValues[j];
    f.write(reinterpret_cast<const char*>(&fX), sizeof(double));
  }
}

void readGrid(const std::string& filename, base::GridStorage& gridStorage) {
  base::DataVector functionValues(0);
  readGrid(filename, gridStorage, functionValues);
}

void readGrid(const std::string& filename, base::GridStorage& gridStorage,
              base::DataVector& functionValues) {
  std::ifstream f;
  f.exceptions(std::ifstream::failbit | std::ifstream::badbit);
  f.open(filename.c_str(), std::ios::in | std::ios::binary);

  size_t N, d;
  f.read(reinterpret_cast<char*>(&N), sizeof(N));
  f.read(reinterpret_cast<char*>(&d), sizeof(d));

  gridStorage.clear();
  functionValues.resize(N);
  base::GridPoint gp(d);

  for (size_t j = 0; j < N; j++) {
    for (size_t t = 0; t < d; t++) {
      double x;
      base::GridPoint::level_type l;
      base::GridPoint::index_type i;
      f.read(reinterpret_cast<char*>(&x), sizeof(x));
      f.read(reinterpret_cast<char*>(&l), sizeof(l));
      f.read(reinterpret_cast<char*>(&i), sizeof(i));
      gp.set(t, l, i);
    }

    double functionValue;
    f.read(reinterpret_cast<char*>(&functionValue), sizeof(functionValue));
    gridStorage.insert(gp);
    functionValues[j] = functionValue;
  }
}

void writeMatrix(const std::string& filename, base::DataMatrix& A) {
  // convert DataMatrix to std::vector
  const std::vector<double> AVector(A.getPointer(), A.getPointer() + A.getNrows() * A.getNcols());
  writeMatrix(filename, AVector, A.getNrows(), A.getNcols());
}

void readMatrix(const std::string& filename, base::DataMatrix& A) {
  std::vector<double> AVector;
  size_t m, n;
  readMatrix(filename, AVector, m, n);
  A.resize(m, n);
  A = base::DataMatrix(&AVector[0], m, n);
}

void writeVector(const std::string& filename, base::DataVector& x) {
  // convert DataVector to std::vector
  const std::vector<double> xVector(x.getPointer(), x.getPointer() + x.getSize());
  writeMatrix(filename, xVector, 1, x.getSize());
}

void readVector(const std::string& filename, base::DataVector& x) {
  std::vector<double> xVector;
  size_t m, n;
  readMatrix(filename, xVector, m, n);
  x.resize(xVector.size());
  x = base::DataVector(&xVector[0], xVector.size());
}
}  // namespace file_io
}  // namespace optimization
}  // namespace sgpp
