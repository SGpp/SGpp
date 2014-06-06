/* ****************************************************************************
* Copyright (C) 2011 Technische Universitaet Muenchen                         *
* This file is part of the SG++ project. For conditions of distribution and   *
* use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
**************************************************************************** */
// @author Sarpkan Selcuk (Sarpkan.Selcuk@mytum.de)

#ifndef STRETCHING_HPP
#define STRETCHING_HPP

#define LOOKUPSIZE 2047
#define LOOKUPMAX 11

#include <cstddef>
#include <math.h>
#include <string>
#include <iostream>
#include <fstream>
#include <vector>
#include <time.h>
#include <algorithm>
#include "base/grid/common/BoundingBox.hpp"
namespace sg {
  namespace base {

    struct Stretching1D {
      std::string type;
      double x_0, xsi;
      //Index 0 : current position
      //Index 1: left neighbor position
      //Index 2: right neighbor position
      double lookup[LOOKUPSIZE][3];
    };

    /**
     *
     * Stretching can be done in four different ways:
     *
     * 1. No Stretching operation is done
     * 2. Stretching is done via log/exp calculation
     * 3. Sinh Stretching
     * 4. Discrete Stretching (No function used, grid points are directly taken)
     */

    class Stretching : public BoundingBox {
      private:

        Stretching1D* stretching1Ds;
        int* discreteVectorLevel;
        std::string* stretchingMode;


        /*
         * generates the lookup table up to level 11 for all dimensions, the rest is calculated via the
         * Stretching type given.
         */
        void generateLookupTable();

        /*
         * generates the positions of the left and right neighbors of one specific level and index point.
         * Written to the lookup table
         */
        void generateLeftRightArrays(Stretching1D& str1D, size_t dim);

        /*
         * The main Stretching function uses the formula
         *
         * \f$ x_i = f^{-1}(f(a) + \frac{i}{n}.(f(b)-f(a)) ), i=0,...,n \f$
         * where f is the transformation function and [a,b] are the intervals the transformation is defined.
         *
         * Calculates the streched coordinate with respect to the given parameters
         *Note to self: may take other than double, take care
         * @param level level of the node
         * @param index index of the node
         * @param dimension dimension of the node
         */
        double stretchingXform(int level, int index, size_t dimension);

        /*
         * calculates the logarithm transform in the given dimension
         *
         * \f$ f = log(x) \text{and} f^{-1}=e^x \f$
         *
         * @param dimension describes in which dimension the Stretching occurs.
         */
        void logXform(Stretching1D& str1D, size_t dimension);

        /*
         * overloaded
         * calculates the logarithm transform for an arbitrary input
         *
         * @param level level of the node
         * @param index index of the node
         * @param dimension dimension of the node
         */
        double logXform(int level, int index, size_t dimension);

        /*
         * calculates the leentvar transform in the given dimension
         *
         * \f$ f= sinh_{-1}(\xi(x-x_0) \text{and} f^{-1}=frac{1}{\xi}sinh(y)+x_0\f$
         *
         * @param dimension describes in which dimension the Stretching occurs.
         */
        void leentvaarXform (Stretching1D& str1D, size_t dimension);

        /*
         * overloaded
         * calculates the leentvar transform in the given dimension
         *
         * \f$ f= sinh_{-1}(\xi(x-x_0) \text{and} f^{-1}=frac{1}{\xi}sinh(y)+x_0\f$
         *
         * @param dimension describes in which dimension the Stretching occurs.
         */
        double leentvaarXform (int level, int index, size_t dimension);

        /*
         * makes no Stretching on the dimension, the points are equidistant
         *
         * @param dimension describes in which dimension the Stretching occurs.
         */
        void noXform(Stretching1D& str1D, size_t dimension);

        /*
         * overloaded
         * returns the position of the node determined by the input arguments, on the non-streched grid.
         *
         * @param level level of the node
         * @param index index of the node
         * @param dimension dimension of the node
         *
         */
        double noXform(int level, int index, size_t dimension);

        /*
         * calculates the lookup table index
         *
         * @param level level of the point
         * @param index index of the point
         */
        int calculateLookupIndex(int level, int index);

        /*
         * gets the discrete points that are stretched and creates a lookup table
         * @param vec double vector containing the discrete points, includes left and right boundary
         * @param stretch1d Stretch1D to write on the type of stretching and lookup table
         * @param dim dimension of the node
         * @param discreteVectorLevel how many levels the grid vector points form
         */

        void parseVectorToLookupTable(std::vector<double> vec, Stretching1D& stretch1d, size_t dim, int& discreteVectorLevel);


        /*
         * calculates the left and right neighbor index and levels
         * @param level current level of the node
         * @param index current index of the node
         * @param leftLevel left neighbor level
         * @param leftIndex left neighbor index
         * @param rightLevel right neighbor level
         * @param rightIndex right neighbor index
         */
        void calculateNeighborSpecs(int level, int index, int& leftLevel, int& leftIndex, int& rightLevel, int& rightIndex);

      public:

        /**
         * Constructor for Stretching
         *
         * initializes the Stretching using the boundaries with the input type array given
         * for each dimension
         *
         * @param dim number of the dimensions used with the grid
         * @param boundaries DimensionBoundary struct to get the boundary values for Stretching
         * @param stretching1ds array to define the stretching type
         */
        Stretching(size_t dim, DimensionBoundary* boundaries, Stretching1D* stretching1ds);


        Stretching(size_t dim, std::vector<DimensionBoundary> boundaries, std::vector<Stretching1D> t);

        /**
         * Constructor for Stretching
         *
         * initializes the Stretching using the coordinates given. (For Janos' request)
         *
         * @param dim number of the dimensions used with the grid
         * @param coordinates vector<double> array to get the boundaries, as well as the coordinates
         * of the specific level the vector defines for each dimension.
         */
        Stretching(size_t dim, std::vector<double>* coordinates);

        /**
         * Copy-Constructor
         *
         * initializes the Stretching with values of another Stretching instance
         *
         * @param copyStretching reference to a Stretching Object whose values are copied
         */
        Stretching(Stretching& copyStretching);

        /**
         * Destructor
         */
        ~Stretching();

        /*
         * Gets the node position defined in the input parameters
         *
         * @param level level of the node
         * @param index index of the node
         * @param dimension dimension of the node
         *
         */
        double getCoordinates(int level, int index, size_t dimension);

        /*
         * Returns the type of the Stretching of the dimension given
         *
         * @param dim dimension index of the type array
         *
         */
        Stretching1D getStretching1D(size_t dim);

        /*
         * Prints the lookup table generated. Used for debugging
         */
        void printLookupTable();

        /*
         * Given the level, index and dimension, returns the lookup entries for left, current and right positions.
         *
         * @param level level of the node
         * @param index index of the node
         * @param dimension dimension of the node
         * @param posc the point itself
         * @param posl left point
         * @param posr right point
         */
        void getAdjacentPositions(int level, int index, size_t dimension, double& posc, double& posl, double& posr);

        /*
         * returns the stretchingMode, can be "analytic" or "discrete"
         */
        std::string* getStretchingMode();

        /*
         * returns the discreteVector for the discrete grid structure.
         * used for serializing the grid
         */
        std::vector<double>* getDiscreteVector(bool bSort);

        /*
         * returns the discretevectorlevel array of size nDim
         */
        int* getDiscreteVectorLevel();

        void calculateNeighborLookup(int maxlevel);
    };

  }
}

#endif /* STRETCHING_HPP */
