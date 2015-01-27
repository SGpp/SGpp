// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#ifndef GRID_HPP
#define GRID_HPP


#include <sgpp/base/operation/hash/OperationEval.hpp>
#include <sgpp/base/operation/hash/common/basis/Basis.hpp>
#include <sgpp/base/grid/generation/GridGenerator.hpp>

#include <map>

#include <sgpp/globaldef.hpp>


namespace SGPP {
  namespace base {

    /**
     * enum to address different gridtypes in a standardized way
     *
     * @todo: use this also for unserialization!
     */
    enum GridType {
      Linear,
      LinearStretched,
      LinearBoundary,
      LinearTruncatedBoundary,
      LinearStretchedTruncatedBoundary,
      LinearGeneralizedTruncatedBoundary,
      ModLinear,
      Poly,
      ModPoly,
      ModWavelet,
      ModBspline,
      Prewavelet,
      SquareRoot,
      Periodic
    };

    /**
     * structure that can be used by applications to cluster regular grid information
     */
    struct RegularGridConfiguration {
      /// Grid Type, see enum
      SGPP::base::GridType type_;
      /// number of dimensions
      size_t dim_;
      /// number of levels
      int level_;
    };

    /**
     * structure that can be used by application to define adaptivity strategies
     */
    struct AdpativityConfiguration {
      /// number of refinements
      size_t numRefinements_;
      /// refinement threshold for surpluses
      double threshold_;
      /// refinement type: false: classic, true: maxLevel
      bool maxLevelType_;
      /// max. number of points to be refined
      size_t noPoints_;
      /// max. percent of points to be refined
      double percent_;
    };

    /**
     * abstract base class for all types grids used in sgpp
     * the class gives pure virtual function definitions that
     * have to be implemented by all types of grids
     */
    class Grid {
      public:
        /**
         * creates a stencil for a linear grid without boundaries
         *
         * @param dim the grid's dimension
         */
        static Grid* createLinearGridStencil(size_t dim);

        /**
         * creates a stencil for a modified linear grid (without boundaries)
         *
         * @param dim the grid's dimension
         */
        static Grid* createModLinearGridStencil(size_t dim);

        /**
         * creates a linear grid without boundaries
         *
         * @param dim the grid's dimension
         */
        static Grid* createLinearGrid(size_t dim);

        /**
         * creates a linear stretched grid without boundaries
         *
         * @param dim the grid's dimension
         */
        static Grid* createLinearStretchedGrid(size_t dim);

        /**
         * creates a linear boundary grid
         *
         * @param dim the grid's dimension
         */
        static Grid* createLinearBoundaryGrid(size_t dim);

        /**
         * creates a linear truncated boundary grid
         *
         * @param dim the grid's dimension
         */
        static Grid* createLinearTruncatedBoundaryGrid(size_t dim);

        /**
         * creates a linearstretched truncated boundary grid
         *
         * @param dim the grid's dimension
         */
        static Grid* createLinearStretchedTruncatedBoundaryGrid(size_t dim);
        /**
         * creates a mod linear grid
         *
         * @param dim the grid's dimension
         */
        static Grid* createModLinearGrid(size_t dim);

        /**
         * creates a mod polynomial grid
         *
         * @param dim the grid's dimension
         * @param degree the polynom's max. degree
         */
        static Grid* createPolyGrid(size_t dim, size_t degree);

        /**
         * creates a poly grid
         *
         * @param dim the grid's dimension
         * @param degree the polynom's max. degree
         */
        static Grid* createModPolyGrid(size_t dim, size_t degree);

        /**
         * creates a mod wavelet grid
         *
         * @param dim the grid's dimension
         */
        static Grid* createModWaveletGrid(size_t dim);

        /**
         * creates a mod-Bspline grid
         *
         * @param dim the grid's dimension
         * @param degree the polynom's max. degree
         */
        static Grid* createModBsplineGrid(size_t dim, size_t degree);

        /**
         * creates a prewavelet grid
         *
         * @param dim the grid's dimension
         */
        static Grid* createPrewaveletGrid(size_t dim);

        /**
         * creates a square root grid(h-grid)
         * @param dim the grid's dimension
         * */
        static Grid* createSquareRootGrid(size_t dim);

        /**
         * creates a truncated boundary grid=contains all the gridpoints of the fullgrids which have \f$|l|<level and li>=l_user\f$
         *
         * @param dim the grid's dimension
         * */
        static Grid* createLinearGeneralizedTruncatedBoundaryGrid(size_t dim);

        /**
         * creates a periodic grid
         *
         * @param dim the grid's dimension
         */
        static Grid* createPeriodicGrid(size_t dim);

        /**
         * reads a grid out of a string
         *
         * @param istr string that contains the grid information
         */
        static Grid* unserialize(const std::string& istr);

        /**
         * reads a grid out of a stream
         * @todo check for empty istream - error message is not very meaningful
         * @param istr inputstream that contains the grid information
         */
        static Grid* unserialize(std::istream& istr);

      protected:
        /**
         * This constructor creates a new GridStorage out of the stream.
         * For derived classes create an own constructor wich takes a std::istream and calls
         * this function. Add your own static unserialize function and add it in typeMap().
         *
         * @param istr inputstream that contains the grid information
         */
        Grid(std::istream& istr);

        /**
         * Standard Constructor
         */
        Grid();

      public:
        /**
         * Desctructor
         */
        virtual ~Grid();

        /**
         * gets a pointer to the GridStorage object
         *
         * @return pointer to the GridStorage obeject
         */
        virtual GridStorage* getStorage();

        /**
         * gets a pointer to the GridStorage's BoundingsBox object
         *
         * @return pointer to the GridStorage's BoundingsBox object
         */
        virtual BoundingBox* getBoundingBox();

        /**
         * gets a pointer to the GridStorage's Stretching object
         *
         * @return pointer to the GridStorage's Stretching object
         */
        virtual Stretching* getStretching();

        /**
         * sets the GridStorage's BoundingsBox pointer to a BoundingBox object
         *
         * @return pointer to the GridStorage's BoundingsBox object
         */
        virtual void setBoundingBox(BoundingBox& bb);

        /**
         * sets the GridStorage's Stretching pointer to a Stretching object
         *
         * @return pointer to the GridStorage's Stretching object
         */
        virtual void setStretching(Stretching& bb);

        /**
         * gets a pointer to GridGenerator object
         *
         * @return pointer to the GrdGenerator object
         */
        virtual GridGenerator* createGridGenerator() = 0;

        /**
         * Returns a string that identifies the grid type uniquely
         *
         * @todo use enum instead!!
         *
         * @return string that identifies the grid type uniquely
         */
        virtual const char* getType() = 0;

        /**
         * Returns the Basis class associated with the grid
         *
         * @return Basis class associated with the grid
         */
        virtual const SBasis& getBasis() = 0;

        /**
         * Serializes grid to a string.
         * Needed for Python compatibility. Calls serialize(std::ostream&).
         *
         * @param ostr string into which the grid is written
         */
        void serialize(std::string& ostr);

        /**
         * Serializes the grid.
         * Override if additional information need to be saved.
         * Call base function before writing anything!
         *
         * @param ostr stream to which the grid is written
         */
        virtual void serialize(std::ostream& ostr);

        /**
         * Serializes grid to a string.
         * Needed for Java compatibility.
         *
         * @returns string into which the grid is written
         */
        std::string serialize();

        /**
         * Refine grid
         * Refine the given number of points on the grid according to the vector
         *
         * @param vector DataVector vector with errors for each basis function or alpha-vector
         * @param numOfPoints integer number of points to refine
         */
        void refine(DataVector* vector, int numOfPoints);

        /**
         * Evaluate the value of function in the point
         *
         * @param alpha DataVector alpha vector of the grid
         * @param point DataVector point where the function should be evaluated
         */
        double eval(DataVector& alpha, DataVector& point);

        /**
         * Insert one point to the grid
         *
         * @param dim dimension of the grid
         * @param levels array with levels of the point
         * @param indices array with indices of the point
         * @param isLeaf indicator whether the point is a leaf
         */
        void insertPoint(size_t dim, unsigned int levels[], unsigned int indices[],
                         bool isLeaf);

        /**
         * Returns the number of points on the grid
         * @return the number of points on the grid
         */
        size_t getSize();

        /**
         * returns the algorithmic dimensions (the dimensions in which the Up Down
         * operations should be applied)
         *
         * @return the algorithmic dimensions
         */
        std::vector<size_t> getAlgorithmicDimensions();

        /**
         * sets the algorithmic dimensions (the dimensions in which the Up Down
         * operations should be applied)
         *
         * @param newAlgoDims std::vector containing the algorithmic dimensions
         */
        void setAlgorithmicDimensions(std::vector<size_t> newAlgoDims);

      protected:
        /// pointer the GridStorage object of the grid
        GridStorage* storage;

        typedef Grid* (*Factory)(std::istream&);
        typedef std::map<std::string, Grid::Factory> factoryMap;

        static Grid* nullFactory(std::istream&);


      private:
        /**
         * This method returns a map with all available grid types for serialization
         *
         * @return a map with all available grid types for serialization
         */
        static factoryMap& typeMap();

        //pointer to the Operation Eval used in Grid.eval()
        static OperationEval* evalOp;
    };

  }
}

#endif /* GRID_HPP */
