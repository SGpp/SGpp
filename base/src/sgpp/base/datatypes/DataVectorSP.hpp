// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#ifndef DATAVECTORSP_HPP
#define DATAVECTORSP_HPP

#include <string>
#include <vector>

#include <sgpp/globaldef.hpp>


namespace SGPP {
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
    class DataVectorSP {
      public:
        /**
         * Create a DataVectorSP with @em size elements (uninitialized values).
         *
         * @param size Number of elements
         */
        DataVectorSP(size_t size);

        /**
         * Create a DataVectorSP with @em size elements and initializes
         * all elements with the same value.
         *
         * @param size Number of elements
         * @param value Value for all entries
         */
        DataVectorSP(size_t size, float value);

        /**
         * Create a new DataVectorSP that is a copy of vec.
         *
         * @param vec Reference to another instance of DataVectorSP
         */
        DataVectorSP(const DataVectorSP& vec);

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
        DataVectorSP(std::vector<float> input);

        /**
         * Create a new DataVectorSP from a std::vector<int>.
         *
         * @param input std::vector<int> that contains the data
         */
        DataVectorSP(std::vector<int> input);

        /**
         * Resizes the DataVectorSP to size elements.
         * All new additional entries are uninitialized.
         * If nrows is smaller than the current number of rows,
         * all superfluous entries are removed.
         *
         * @param size New number of elements of the DataVectorSP
         */
        void resize(size_t size);

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
         * Add add potentially new elements to the DataVectorSP. The size remains unchanged
         * Reserves memory for potentially inc_elems new elements;
         * the actual number of elements remains unchanged.
         * Corresponds to a resize to size+inc_elems new elements while leaving
         * the current vector's size unchanged.
         *
         * @param inc_elems Number of additional elements for which storage is to be reserved.
         */
        void addSize(size_t inc_elems);

        /**
         * Appends a new element and returns index of it.
         * If the new element does not fit into the reserved memory,
         * reserves memory for getInc() additional elements.
         *
         * @return Index of new element
         */
        size_t append();

        /**
         * Appends a new element and returns index of new element.
         * If the new element does not fit into the reserved memory,
         * reserves memory for getInc() additional elements.
         *
         * @param value Value of new element
         * @return Index of new element
         */
        size_t append(float value);

        /**
         * Inserts a new element at the given index.
         * If the new element does not fit into the reserved memory,
         * reserves memory for getInc() additional elements.
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
         * Copies the data from another DataVectorSP.
         * Dimensions have to match.
         *
         * @param vec the DataVectorSP containing the data
         * @return *this
         */
        DataVectorSP& operator=(const DataVectorSP& vec);

        /**
         * Returns a reference to the i-th element.
         *
         * @param i position of the element
         * @return data[i]
         */
        inline float& operator[](size_t i) {
          return data[i];
        }

        /**
         * Returns a constant reference to the i-th element.
         *
         * @param i position of the element
         * @return data[i]
         */
        inline const float& operator[](size_t i) const {
          return data[i];
        }

        /**
         * Returns the i-th element.
         *
         * @param i position of the element
         * @return data[i]
         */
        inline float get(size_t i) const {
          return data[i];
        }

        /**
         * Sets the element at index i to value.
         *
         * @param i Index
         * @param value New value for element
         */
        void set(size_t i, float value);

        /**
         * Adds the values from another DataVectorSP to the current values.
         * Modifies the current values.
         *
         * @param vec The DataVectorSP which is added to the current values
         */
        void add(DataVectorSP& vec);

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
        void componentwise_mult(DataVectorSP& vec);

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
        void componentwise_div(DataVectorSP& vec);

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
         * Returns the dot product of the two vectors.
         *
         * @param vec Reference to another vector
         *
         * @return The dot-product
         */
        float dotProduct(DataVectorSP& vec) const;

        /**
         * gets a pointer to the data array
         *
         * @return pointer to the data array
         */
        float* getPointer();

        /**
         * gets the elements stored in the vector
         *
         * @return elements stored in the vector
         */
        inline size_t getSize() const {
          return size;
        }

        /**
         * Returns the number of unused elements.
         *
         * @return number of unused elements
         */
        inline size_t getUnused() const {
          return unused;
        }

        /**
         * Determines the number of non-zero elements in the vector.
         *
         * @return The number of non-zero elements
         */
        size_t getNumberNonZero() const;

        /**
         * Get the current number of elements by which the DataVectorSP is extended,
         * if append() or insert() is called and no unused rows are left
         *
         * @return Increment
         */
        inline size_t getInc() const {
          return inc_elems;
        }

        /**
         * Sets the current number of elements by which the DataVectorSP is extended,
         * if append() or insert() is called and no unused elements are left.
         * Defaults to 100.
         *
         * @param inc_elems Increment
         */
        void setInc(size_t inc_elems) {
          this->inc_elems = inc_elems;
        }

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

        /**
         * Destructor
         */
        virtual ~DataVectorSP();

      private:
        /// Array to store the data
        float* data;
        /// Number of elements of the data vector
        size_t size;
        /// Number of additional rows for which memory has already been reserved
        size_t unused;
        /// Number of elements by which the reserved memory is increased, if adding a row would exceed the storage reserved so far.
        size_t inc_elems;
    };

  }
}

#endif /* DATAVECTORSP_HPP */
