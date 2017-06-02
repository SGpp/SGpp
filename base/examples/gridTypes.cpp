// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

/**
 * \page GridTypes List of different Grid Types
 *
 * This example is supposed to simply demonstrate the available grid, boundary and basis function
 * types.
 */

/**
 * First, we include the meta-header <tt>sgpp_base.hpp</tt>, which includes itself all headers
 * from the base module and set a few parameters.
 */

#include <sgpp_base.hpp>

int main() {
  size_t dim = 2;
  size_t level = 3;
  size_t polyDegree = 4;
  sgpp::base::Grid* grid;

  /**
   * LinearGrid
   * @image html "createLinearGrid_C2J-small.png" ""
   */

  grid = sgpp::base::Grid::createLinearGrid(dim);
  grid->getGenerator().regular(level);

  /**
   * LinearBoundaryGrid
   * @image html "createLinearBoundaryGrid_C2,_1J-small.png"   ""
   */
  grid = sgpp::base::Grid::createLinearBoundaryGrid(dim);
  grid->getGenerator().regular(level);

  /**
   * LinearStretchedGrid
   * @image html createLinearStretchedGrid_C2J-small.png ""
   */
  grid = sgpp::base::Grid::createLinearStretchedGrid(dim);
  grid->getGenerator().regular(level);

  /**
   * LinearStretchedBoundaryGrid
   * @image html createLinearStretchedBoundaryGrid_C2J-small.png ""
   */
  grid = sgpp::base::Grid::createLinearStretchedBoundaryGrid(dim);
  grid->getGenerator().regular(level);

  /**
   * WaveletGrid
   * @image html createWaveletGrid_C2J-small.png ""
   */
  grid = sgpp::base::Grid::createWaveletGrid(dim);
  grid->getGenerator().regular(level);

  /**
   * WaveletBoundaryGrid
   * @image html createWaveletBoundaryGrid_C2J-small.png ""
   */
  grid = sgpp::base::Grid::createWaveletBoundaryGrid(dim);
  grid->getGenerator().regular(level);

  /**
   * SquareRootGrid
   * @image html createSquareRootGrid_C2J-small.png ""
   */
  grid = sgpp::base::Grid::createSquareRootGrid(dim);
  grid->getGenerator().regular(level);

  /**
   * PrewaveletGrid
   * @image html createPrewaveletGrid_C2J-small.png ""
   */
  grid = sgpp::base::Grid::createPrewaveletGrid(dim);
  grid->getGenerator().regular(level);

  /**
   * PolyGrid
   * @image html "createPolyGrid_C2,_3J-small.png" ""
   */
  grid = sgpp::base::Grid::createPolyGrid(dim, polyDegree);
  grid->getGenerator().regular(level);

  /**
   * PolyBoundaryGrid
   * @image html "createPolyBoundaryGrid_C2,_3J-small.png" ""
   */
  grid = sgpp::base::Grid::createPolyBoundaryGrid(dim, polyDegree);
  grid->getGenerator().regular(level);

  /**
   * PeriodicGrid
   * @image html createPeriodicGrid_C2J-small.png ""
   */
  grid = sgpp::base::Grid::createPeriodicGrid(dim);
  grid->getGenerator().regular(level);

  /**
   * ModWaveletGrid
   * @image html createModWaveletGrid_C2J-small.png ""
   */
  grid = sgpp::base::Grid::createModWaveletGrid(dim);
  grid->getGenerator().regular(level);

  /**
   * ModPolyGrid
   * @image html "createModPolyGrid_C2,_3J-small.png" ""
   */
  grid = sgpp::base::Grid::createModPolyGrid(dim, polyDegree);
  grid->getGenerator().regular(level);

  /**
   * ModLinearGridStencil
   * @image html createModLinearGridStencil_C2J-small.png ""
   */
  grid = sgpp::base::Grid::createModLinearGridStencil(dim);
  grid->getGenerator().regular(level);

  /**
   * ModLinearGrid
   * @image html createModLinearGrid_C2J-small.png ""
   */
  grid = sgpp::base::Grid::createModLinearGrid(dim);
  grid->getGenerator().regular(level);

  /**
   * ModFundamentalSplineGrid
   * @image html "createModFundamentalSplineGrid_C2,_3J-small.png" ""
   */
  grid = sgpp::base::Grid::createModFundamentalSplineGrid(dim, polyDegree);
  grid->getGenerator().regular(level);

  /**
   * ModBsplineGrid
   * @image html "createModBsplineGrid_C2,_3J-small.png" ""
   */
  grid = sgpp::base::Grid::createModBsplineGrid(dim, polyDegree);
  grid->getGenerator().regular(level);

  /**
   * ModBsplineClenshawCurtisGrid
   * @image html "createModBsplineClenshawCurtisGrid_C2,_3J-small.png" ""
   */
  grid = sgpp::base::Grid::createModBsplineClenshawCurtisGrid(dim, polyDegree);
  grid->getGenerator().regular(level);

  /**
   * LinearTruncatedBoundaryGrid
   * @image html createLinearTruncatedBoundaryGrid_C2J-small.png ""
   */
  grid = sgpp::base::Grid::createLinearTruncatedBoundaryGrid(dim);
  grid->getGenerator().regular(level);

  /**
   * LinearGridStencil
   * @image html createLinearGridStencil_C2J-small.png ""
   */
  grid = sgpp::base::Grid::createLinearGridStencil(dim);
  grid->getGenerator().regular(level);

  /**
   * LinearClenshawCurtisGrid
   * @image html createLinearClenshawCurtisGrid_C2J-small.png ""
   */
  grid = sgpp::base::Grid::createLinearClenshawCurtisGrid(dim);
  grid->getGenerator().regular(level);

  /**
   * FundamentalSplineGrid
   * @image html "createFundamentalSplineGrid_C2,_3J-small.png" ""
   */
  grid = sgpp::base::Grid::createFundamentalSplineGrid(dim, polyDegree);
  grid->getGenerator().regular(level);

  /**
   * BsplineGrid
   * @image html "createBsplineGrid_C2,_3J-small.png" ""
   */
  grid = sgpp::base::Grid::createBsplineGrid(dim, polyDegree);
  grid->getGenerator().regular(level);

  /**
   * BsplineClenshawCurtisGrid
   * @image html "createBsplineClenshawCurtisGrid_C2,_3J-small.png" ""
   */
  grid = sgpp::base::Grid::createBsplineClenshawCurtisGrid(dim, polyDegree);
  grid->getGenerator().regular(level);

  /**
   * BsplineBoundaryGrid
   * @image html "createBsplineBoundaryGrid_C2,_3J-small.png" ""
   */
  grid = sgpp::base::Grid::createBsplineBoundaryGrid(dim, polyDegree);
  grid->getGenerator().regular(level);
}
