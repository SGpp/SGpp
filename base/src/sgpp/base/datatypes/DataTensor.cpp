// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <sgpp/base/datatypes/DataTensor.hpp>
#include <sgpp/base/exception/data_exception.hpp>

namespace sgpp {
namespace base {

DataTensor::DataTensor() : DataTensor(0, 0, 0) {}

DataTensor::DataTensor(size_t ndepth, size_t nrows, size_t ncols)
    : DataTensor(ndepth, nrows, ncols, 0.0) {}

DataTensor::DataTensor(size_t ndepth, size_t nrows, size_t ncols, double value)
    : ndepth(ndepth), nrows(nrows), ncols(ncols), matrixsize(nrows * ncols) {
  this->assign(ndepth * nrows * ncols, value);
}

void DataTensor::resize(size_t ndepth, size_t nrows, size_t ncols) {
  this->std::vector<double>::resize(ndepth * nrows * ncols);
}

double DataTensor::get(size_t depth, size_t row, size_t col) {
  if ((depth >= ndepth) || (row >= nrows) || (col >= ncols)) {
    throw sgpp::base::data_exception("DataTensor::get Entry does not exist");
  }
  return (*this)[depth * matrixsize + row * ncols + col];
}

void DataTensor::set(size_t depth, size_t row, size_t col, double value) {
  if ((depth >= ndepth) || (row >= nrows) || (col >= ncols)) {
    throw sgpp::base::data_exception("DataTensor::set Entry does not exist");
  }
  (*this)[depth * matrixsize + row * ncols + col] = value;
}

void DataTensor::getMatrix(size_t depth, sgpp::base::DataMatrix& matrix) {
  if (depth >= ndepth) {
    throw sgpp::base::data_exception("DataTensor::getMatrix Entry does not exist");
  }
  if (matrix.getNcols() != this->ncols) {
    matrix.resize(this->nrows, this->ncols);
  }
  if (matrix.getNrows() != this->nrows) {
    matrix.resize(this->nrows, this->ncols);
  }

  sgpp::base::DataVector row(this->ncols);
  for (size_t i = 0; i < this->nrows; i++) {
    this->getRow(depth, i, row);
    matrix.setRow(i, row);
  }
}

void DataTensor::getRow(size_t depth, size_t row, sgpp::base::DataVector& vec) {
  if (row >= nrows) {
    throw sgpp::base::data_exception("DataTensor::getRow Entry does not exist");
  }
  vec.clear();
  for (size_t i = 0; i < this->ncols; ++i) {
    vec.push_back((*this)[depth * matrixsize + row * ncols + i]);
  }
}

void DataTensor::getColumn(size_t depth, size_t col, sgpp::base::DataVector& vec) {
  if (col >= ncols) {
    throw sgpp::base::data_exception("DataTensor::getColumn Entry does not exist");
  }
  if (vec.getSize() != this->nrows) {
    vec.resize(nrows);
  }

  for (size_t j = 0; j < this->nrows; ++j) {
    vec[j] = (*this)[depth * matrixsize + j * ncols + col];
  }
}

void DataTensor::toString(std::string& text) const {
  std::stringstream str;

  str << std::scientific;
  str.precision(20);

  str << "{";
  for (size_t h = 0; h < ndepth; ++h) {
    str << "\n[";
    for (size_t i = 0; i < nrows; ++i) {
      str << "[";

      for (size_t j = 0; j < ncols; ++j) {
        if (j != 0) {
          str << ", ";
          // add linebreak for readability
          if (j % 20 == 0) {
            str << std::endl;
          }
        }

        str << (*this)[h * matrixsize + i * ncols + j];
      }

      if (i == nrows - 1) {
        str << "]";
      } else {
        str << "]," << std::endl;
      }
    }
    str << "]\n";
  }
  str << "}";
  text = str.str();
}

std::string DataTensor::toString() const {
  std::string str;
  toString(str);
  return str;
}

}  // namespace base
}  // namespace sgpp
