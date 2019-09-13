// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#ifndef DATAVECTORSP_HPP
#define DATAVECTORSP_HPP

#include <sgpp/globaldef.hpp>

#include <initializer_list>
#include <string>
#include <vector>

namespace sgpp {
namespace base {

/**
 * A class to store one-dimensional data.
 * Typically, an object of type DataVectorSP will contain an array
 * of (hierarchical) coefficients (or surplusses), or the coordinates
 * of a data point at which a sparse grid function should be
 * evaluated.
 *
 * This is an re-implementation of the standard DataVector
 * for single precision floating point numbers in order to
 * increase support for GPUs.
 */
class DataVectorSP : public std::vector<float> {
 public:
  /**
   * Create an empty DataVectorSP.
   */
  DataVectorSP();

  /**
   * Copy constructor.
   */
  DataVectorSP(const DataVectorSP&) = default;

  /**
   * Move constructor
   */
  DataVectorSP(DataVectorSP&&) = default;

  /**
   * Copy assignment operator
   */
  DataVectorSP& operator=(const DataVectorSP&) = default;

  /**
   * Move assignment operator
   */
  DataVectorSP& operator=(DataVectorSP&&) = default;

  /**
   * Destructor
   */
  ~DataVectorSP() = default;

  /**
   * Create a DataVectorSP with @em size elements (uninitialized values).
   *
   * @param size Number of elements
   */
  explicit DataVectorSP(size_t size);

  /**
   * Create a DataVectorSP with @em size elements and initializes
   * all elements with the same value.
   *
   * @param size Number of elements
   * @param value Value for all entries
   */
  DataVectorSP(size_t size, float value);

  /**
   * Create a new DataVectorSP from a float array with size elements.
   *
   * @param input float array that contains the data
   * @param size number of elements
   */
  DataVectorSP(float* input, size_t size);

  /**
   * Create a new DataVectorSP from a std::vector<float>.
   *
   * @param input std::vector<float> that contains the data
   */
  explicit DataVectorSP(std::vector<float> input);

  /**
   * Create a new DataVector from a std::initializer_list<float>.
   *
   * @param input std::initializer_list<float> that contains the data
   */
  explicit DataVectorSP(std::initializer_list<float> input);

  /**
   * Create a new DataVectorSP from a std::vector<int>.
   *
   * @param input std::vector<int> that contains the data
   */
  explicit DataVectorSP(std::vector<int> input);

  static DataVectorSP fromFile(const std::string& fileName);

  static DataVectorSP fromString(const std::string& serializedVector);

  /**
   * Resizes the DataVectorSP to size elements.
   * All new additional entries are set to zero.
   * If nrows is smaller than the current number of rows,
   * all superfluous entries are removed.
   *
   * @param nrows New number of rows of the DataVectorSP
   */
  void resizeZero(size_t nrows);

  /**
   * Resizes the DataVectorSP by removing entries. Throws an exception
   * if boundaries a violated.
   *
   * @param remainingIndex vector that contains the remaining indices of the DataVectorSP
   */
  void restructure(std::vector<size_t>& remainingIndex);

  /**
   * Removes indexes form the vector. Throws an exception if the boundaries are violated
   *
   * @param indexesToRemove a vector if indexes that will be removed from the vector
   */
  void remove(std::vector<size_t>& indexesToRemove);

  /**
   * Appends a new element and returns index of it.
   *
   * @return Index of new element
   */
  size_t append();

  /**
   * Appends a new element and returns index of new element.
   *
   * @param value Value of new element
   * @return Index of new element
   */
  size_t append(float value);

  /**
   * Inserts a new element at the given index.
   *
   * @param index Index of new element
   * @param value Value of new element
   */
  void insert(size_t index, float value);

  /**
   * Sets all values of DataVectorSP to value
   *
   * @param value New value for all entries
   */
  void setAll(float value);

  /**
   * Returns the i-th element.
   *
   * @param i position of the element
   * @return data[i]
   */
  inline float get(size_t i) const { return (*this)[i]; }

  /**
   * Sets the element at index i to value.
   *
   * @param i Index
   * @param value New value for element
   */
  void set(size_t i, float value);

  /**
   * Copies the data from another DataVectorSP vec.
   * Disregards the number of entries set for the two vectors,
   * i.e., just copies the data entry by entry.
   * If the size matches, the current DataVectorSP is an
   * exact copy of vec. If not, as many elements as possible are
   * copied, and everything else is left untouched.
   *
   * @param vec The source DataVectorSP containing the data
   */
  void copyFrom(const DataVectorSP& vec);

  /**
   * Adds the values from another DataVectorSP to the current values.
   * Modifies the current values.
   *
   * @param vec The DataVectorSP which is added to the current values
   */
  void add(const DataVectorSP& vec);

  /***
   * Accumulation (summation) of vectors using Kahan's summation formula
   * for better precision.
   *
   * All vectors need to be added to one for the summation formula to work. It
   * will not work in a tree-like summation pattern.
   *
   * @param vec The DataVectorSP that will be added to the current one
   */
  void accumulate(const DataVectorSP& vec);

  /**
   * Subtracts the values from another DataVectorSP of the current values.
   * Modifies the current values.
   *
   * @param vec The DataVectorSP which is subtracted from the current values
   */
  void sub(const DataVectorSP& vec);

  /**
   * Multiplies the current DataVectorSP component-wise with another DataVectorSP.
   * Modifies the current values.
   * Performs
   * @code
   * for i from 1 to this.getSize()
   *   this[i] *= vec[i]
   * @endcode
   *
   * @param vec the DataVectorSP which is multiplied to current DataVectorSP
   */
  void componentwise_mult(const DataVectorSP& vec);

  /**
   * Divides the current DataVectorSP component-wise by another DataVectorSP.
   * Modifies the current values.
   * Performs
   * @code
   * for i from 1 to this.getSize()
   *   this[i] /= vec[i]
   * @endcode
   * Note: <b>No check for division by zero!</b>
   *
   * @param vec the DataVectorSP which the current DataVectorSP is divided by
   */
  void componentwise_div(const DataVectorSP& vec);

  /**
   * Returns the dot product of the two vectors.
   *
   * @param vec Reference to another vector
   *
   * @return The dot-product
   */
  float dotProduct(const DataVectorSP& vec) const;

  /**
   * multiplies all elements by a constant factor
   *
   * @param scalar the constant
   */
  void mult(float scalar);

  /**
   * Squares all elements of the DataVectorSP
   */
  void sqr();

  /**
   * Takes the square root of all elements of the DataVectorSP
   */
  void sqrt();

  /**
   * Sets all elements to their absolute value.
   *
   */
  void abs();

  /**
   * Returns the sum of all elements
   *
   * @return The sum of all elements
   */
  float sum() const;

  /**
   * calculates the vector's max norm
   *
   * @return the vector's max norm
   */
  float maxNorm() const;

  /**
   * Returns the vector's root mean square (RMS)-norm, i.e.,
   * @f$\sqrt{ 1/N \sum_{i=1}^N x_i^2 }@f$. If the vector's entries
   * correspond to function values on a full grid, this is the
   * discrete @f$L^2@f$-norm of the corresponding function.
   *
   * @return The vector's root mean square-norm.
   */
  float RMSNorm() const;

  /**
   * Returns the vector's @f$l^2@f$-norm, i.e.,
   * @f$\sqrt{ \sum_i x_i^2 }@f$.
   *
   * @return The vector's @f$l^2@f$-norm.
   */
  float l2Norm() const;

  /**
   * Returns the minimum over all entries.
   *
   * @return Minimal value
   */
  float min() const;

  /**
   * Returns the maximum over all entries.
   *
   * @return global maximum
   */
  float max() const;

  /**
   * Determines minimum and maximum over all entries.
   *
   * @param min Reference variable for the minimum
   * @param max Reference variable for the maximum
   */
  void minmax(float* min, float* max) const;

  /**
   * Adds a*x to current vector.
   * BLAS Level 1 (elementary vector operations) operation: axpy.
   *
   * @param a A scalar
   * @param x Reference to the DataVectorSP
   */
  void axpy(float a, DataVectorSP& x);

  /**
   * gets a pointer to the data array
   *
   * @return pointer to the data array
   */
  float* getPointer();

  /**
   * gets a const pointer to the data array
   *
   * @return const pointer to the data array
   */
  const float* getPointer() const;

  /**
   * gets the elements stored in the vector
   * \deprecated in favour of the equivalent size() method
   *
   * @return elements stored in the vector
   */
  inline size_t getSize() const { return this->size(); }

  /**
   * Determines the number of non-zero elements in the vector.
   *
   * @return The number of non-zero elements
   */
  size_t getNumberNonZero() const;

  /**
   * Partitions vector into two classes using a choosen border.
   *
   * @param threshold value of the border
   */
  void partitionClasses(float threshold);

  /**
   * Normalizes vector entries to [0,1]
   *
   */
  void normalize();

  /**
   * Normalizes vector entries to [border, 1-border]
   *
   * @param border width of border
   */
  void normalize(float border);

  /**
   * Writes the data stored in the DataVectorSP into a string
   *
   * @param text string to which the data is written
   */
  void toString(std::string& text) const;

  /**
   * Returns a description of the DataVectorSP as a string.
   *
   * @returns string of the DataVectorSP
   */
  std::string toString() const;

  void toFile(const std::string& fileName) const;

 private:
  using std::vector<float>::insert;
  /// Corrections for Kahan's summation in accumulate()
  std::vector<float> correction;
};

}  // namespace base
}  // namespace sgpp

#endif /* DATAVECTORSP_HPP */
