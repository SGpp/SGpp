/* ****************************************************************************
* Copyright (C) 2010 Technische Universitaet Muenchen                         *
* This file is part of the SG++ project. For conditions of distribution and   *
* use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
**************************************************************************** */
// @author Alexander Heinecke (Alexander.Heinecke@mytum.de)

#ifndef DATAMATRIXSP_HPP
#define DATAMATRIXSP_HPP

#include <string>
#include <vector>
#include <sgpp/base/datatypes/DataVectorSP.hpp>

#include <sgpp/globaldef.hpp>


namespace SGPP {
  namespace base {

    /**
     * A class to store two-dimensional data.
     * Typically, DataMatrixSP would contain a set of (d-dimensional) data or evaluation points, i.e.,
     * the DataMatrixSP consists of d columns, and each row is one of the points.
     * Thus, typical functionality like obtaining the maximum for a certain dimension (or attribute),
     * or normalizing all data points to the unit interval for a certain dimension are
     * provided.
     *
     * This is an re-implementation of the standard DataMatrixSP
     * for single precision floating point numbers in order to
     * increase support for GPUs.
     */
    class DataMatrixSP {
      public:
        /**
         * Create a two-dimensional DataMatrixSP with @em nrows rows and
         * @em ncols columns.
         *
         * @param nrows Number of rows
         * @param ncols Number of columns
         */
        DataMatrixSP(size_t nrows, size_t ncols);

        /**
         * Create a new DataMatrixSP that is a copy of matr.
         *
         * @param matr Reference to another instance of DataMatrixSP
         */
        DataMatrixSP(const DataMatrixSP& matr);

        /**
         * Create a new DataMatrixSP from a float array.
         * The float array contains the entries row-wise:
         * x0_0,x0_1,...,x0_ncol-1,
         * x1_0,x1_1,...
         * ...
         * xnrow_0, xnrow_1,...,xnrow_ncol-1
         *
         * @param input float array that contains the data
         * @param nrows number of rows
         * @param ncols number of columns
         */
        DataMatrixSP(float* input, size_t nrows, size_t ncols);


        /**
         * Resizes the DataVectorSP to size elements.
         * All new additional entries are uninitialized.
         * If nrows is smaller than the current number of rows,
         * all superfluous entries are removed.
         *
         * @param nrows New number of rows of the DataMatrixSP
         */
        void resize(size_t nrows);

        /**
         * Resizes the DataMatrixSP to nrows rows.
         * All new additional entries are set to zero.
         * If nrows is smaller than the current number of rows,
         * all superfluous entries are removed.
         *
         * @param nrows New number of rows of the DataMatrixSP
         */
        void resizeZero(size_t nrows);

        /**
         * Reserves memory for potentially inc_nrows new rows;
         * the actual number of rows remains unchanged.
         * Corresponds to a resize to nrows+inc_nrows new rows while leaving
         * the current matrix' size unchanged.
         *
         * @param inc_nrows Number of additional rows for which storage is to be reserved.
         */
        void addSize(size_t inc_nrows);

        /**
         * Appends a new row and returns index of it.
         * If the new row does not fit into the reserved memory,
         * reserves memory for getIncRows() additional rows.
         *
         * @return Index of new row
         */
        size_t appendRow();

        /**
           * Appends a new row with data contained in DataVectorSP vec
           * and returns index of new row.
           * If the new row does not fit into the reserved memory,
           * reserves memory for getIncRows() additional rows.
           *
         * @param vec DataVectorSP (length has to match getNcols()) with data
         * @return Index of new row
         */
        size_t appendRow(DataVectorSP& vec);


        /**
        * Sets all entries of DataMatrixSP to value.
        *
        * @param value New value for all entries
        */
        void setAll(float value);

        /**
         * Copies the data from another DataMatrixSP matr.
         * Disregards the number of rows and columns set for the two matrices,
         * i.e., just copies the data entry by entry (and row by row).
         * If the dimensions match (nrows, ncols), the current DataMatrixSP is an
         * exact copy of matr. If not, as many elements as possible are
         * copied, and everything else is left untouched.
         *
         * @param matr The source DataMatrixSP containing the data
         */
        void copyFrom(const DataMatrixSP& matr);

        /**
         * Transposes this DataMatrixSP
         */
        void transpose();

        /**
         * Copies the data from another DataMatrixSP.
         * Dimensions have to match.
         *
         * @param matr the DataMatrixSP containing the data
         */
        DataMatrixSP& operator=(const DataMatrixSP& matr);

        /**
         * Returns the i-th element.
         * For the 5th element in the third row, i would be 2*getNcols()+4.
         *
         * @param i position of the element
         */
        inline float& operator[](size_t i) {
          return data[i];
        };

        /**
         * Returns the value of the element at position [row,col]
         *
         * @param row Row
         * @param col Column
         * @return Value of the element
         */
        inline float get(size_t row, size_t col) const {
          return data[row * ncols + col];
        };

        /**
         * Sets the element at position [row,col] to value.
         *
         * @param row Row
         * @param col Column
         * @param value New value for element
         */
        inline void set(size_t row, size_t col, float value) {
          data[row * ncols + col] = value;
        };

        /**
         * Copies the values of a row to the DataVectorSP vec.
         *
         * @param row The row
         * @param vec DataVectorSP into which the data is written
         */
        void getRow(size_t row, DataVectorSP& vec);

        /**
         * Copies the values of a row to the std::vector vec.
         *
         * @param row The row
         * @param vec std::vector into which the data is written
         */
        void getRow(size_t row, std::vector<float>& vec);

        /**
         * Sets a row of the DataMatrixSP to the values of a DataVectorSP vec.
         *
         * @param row The row which is to be overwritten
         * @param vec DataVectorSP containing the data of the row
         */
        void setRow(size_t row, DataVectorSP& vec);

        /**
           * Copies the values of a column to the DataVectorSP vec.
         *
         * @param col The column
         * @param vec DataVectorSP into which the data is written
         */
        void getColumn(size_t col, DataVectorSP& vec);

        /**
           * Sets a column of the DataMatrixSP to the values of a DataVectorSP vec.
         *
         * @param col The column which is to be overwritten
         * @param vec DataVectorSP containing the data of the column
         */
        void setColumn(size_t col, DataVectorSP& vec);


        /**
         * Adds the values from another DataMatrixSP to the current values.
         * Modifies the current values.
         *
         * @param matr The DataMatrixSP which is added to the current values
         */
        void add(DataMatrixSP& matr);

        /**
           * Subtracts the values from another DataMatrixSP of the current values.
           * Modifies the current values.
         *
         * @param matr The DataMatrixSP which is subtracted from the current values
         */
        void sub(const DataMatrixSP& matr);

        /**
         * Reduce the DataMatrix along the
         * columns by adding all entries in one row.
         *
         * @param reduction DataVector into which the reduce columns are stored
         */
        void addReduce(DataVectorSP& reduction);

        /**
         * expands a given DataVector into a
         * DataMatrix.
         *
         * @param expand DataVector that should be expanded
         */
        void expand(const DataVectorSP& expand);

        /**
         * Multiplies the current DataMatrixSP component-wise with another DataMatrixSP.
         * Modifies the current values.
         * Performs
         * @code
         * for i from 1 to this.getSize()
         *   this[i] *= matr[i]
         * @endcode
         *
         * @param matr the DataMatrixSP which is multiplied to current DataMatrixSP
         */
        void componentwise_mult(DataMatrixSP& matr);

        /**
         * Divides the current DataMatrixSP component-wise by another DataMatrixSP.
         * Modifies the current values.
         * Performs
         * @code
         * for i from 1 to this.getTotalSize()
         *   this[i] /= matr[i]
         * @endcode
         * Note: <b>No check for division by zero!</b>
         *
         * @param matr the DataMatrixSP which the current DataMatrixSP is divided by
         */
        void componentwise_div(DataMatrixSP& matr);

        /**
         * Multiplies all elements by a constant factor
         *
         * @param scalar the constant
         */
        void mult(float scalar);

        /**
         * Squares all elements of the DataMatrixSP
         */
        void sqr();

        /**
         * Takes the square root of all elements of the DataMatrixSP
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
        float sum();

        /**
         * Returns the minimum value of column col.
         *
         * @param col Number of the column
         *
         * @return Minimum value
         */
        float min(size_t col);

        /**
        * Returns the minimum over all entries.
        *
        * @return Minimal value of all entries
        */
        float min();

        /**
        * Returns the maximum value of column col.
        *
        * @param col Number of the column
        *
        * @return Maximum value
        */
        float max(size_t col);

        /**
        * Returns the maximum over all entries.
        *
        * @return Maximal value of all entries
        */
        float max();

        /**
        * Determines minimum and maximum of column col.
        *
        * @param col Number of the column
        * @param min Reference variable for the minimum
        * @param max Reference variable for the maximum
        */
        void minmax(size_t col, float* min, float* max);

        /**
         * Determines minimum and maximum over all entries.
         *
         * @param min Reference variable for the minimum
         * @param max Reference variable for the maximum
         */
        void minmax(float* min, float* max);


        /**
        * Returns pointer to float array containing underlying data.
        *
        * @return Pointer to data
        */
        float* getPointer();

        /**
         * Returns the total number of (used) elements, i.e., getNrows()*getNCols()
         *
         * @return Number of elements stored in the matrix
         */
        inline size_t getSize() {
          return ncols * nrows;
        };

        /**
         * Returns the number of unused rows.
         *
         * @return number of unused rows
         */
        inline size_t getUnused() {
          return unused;
        };

        /**
        * Determines the number of non-zero elements in the vector.
        *
        * @return The number of non-zero elements
        */
        size_t getNumberNonZero();

        /**
         * Returns the number of rows of the DataMatrixSP.
         *
         * @return Number of rows
         */
        inline size_t getNrows() const {
          return nrows;
        }

        /**
         * Returns the number of columns of the DataMatrixSP.
         *
         * @return Number of columns
         */
        inline size_t getNcols() const {
          return ncols;
        }

        /**
         * Get the current number of rows by which the DataMatrixSP is extended,
         * if appendRow() is called and no unused rows are left
         *
         * @return Row increment
         */
        inline size_t getInc() const {
          return inc_rows;
        }

        /**
         * Sets the current number of rows by which the DataMatrixSP is extended,
         * if appendRow() is called and no unused rows are left.
         * Defaults to 100.
         *
         * @param inc_rows Row increment
         */
        void setInc(size_t inc_rows) {
          this->inc_rows = inc_rows;
        }


        /**
         * Normalizes the d-th dimension (entries in the d-th column) to @f$[0,1]@f$.
         * Considers contents of DataMatrixSP as a d-dimensional dataset, one
         * data point per row.
         *
         * @param d The dimension (column) that should be normalized (starting with 0)
         */
        void normalizeDimension(size_t d);

        /**
           * Normalizes the d-th dimension (entries in the d-th column) to @f$[border,1-border]@f$.
           * Considers contents of DataMatrixSP as a d-dimensional dataset, one
           * data point per row.
           *
           * @param d The dimension (column) that should be normalized (starting with 0)
         * @param border Width of the border
         */
        void normalizeDimension(size_t d, float border);


        /**
        * Writes the data stored in the DataMatrixSP into a string
        *
        * @param text String to which the data is written
        */
        void toString(std::string& text);

        /**
        * Returns a description of the DataMatrixSP as a string.
        *
        * @returns string of the DataMatrixSP
        */
        std::string toString();

        /**
        * Destructor
        */
        virtual ~DataMatrixSP();

      private:
        /// Pointer to the data
        float* data;
        /// Number of rows of the data matrix
        size_t nrows;
        /// Number of columns of the data matrix
        size_t ncols;
        /// Number of additional rows for which memory has already been reserved
        size_t unused;
        /// Number of rows by which the reserved memory is increased, if adding a row would exceed the storage reserved so far.
        size_t inc_rows;
    };
  }
}

#endif /* DATAMATRIXSP_HPP */
