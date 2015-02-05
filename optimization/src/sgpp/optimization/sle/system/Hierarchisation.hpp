// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#ifndef SGPP_OPTIMIZATION_SLE_SYSTEM_HIERARCHISATION_HPP
#define SGPP_OPTIMIZATION_SLE_SYSTEM_HIERARCHISATION_HPP

#include <sgpp/globaldef.hpp>

#include <sgpp/optimization/sle/system/Cloneable.hpp>
#include <sgpp/base/grid/Grid.hpp>
#include <sgpp/base/grid/GridStorage.hpp>

#include <sgpp/base/operation/hash/common/basis/BsplineBasis.hpp>
#include <sgpp/base/operation/hash/common/basis/BsplineBoundaryBasis.hpp>
#include <sgpp/base/operation/hash/common/basis/BsplineClenshawCurtisBasis.hpp>
#include <sgpp/base/operation/hash/common/basis/BsplineModifiedBasis.hpp>
#include <sgpp/base/operation/hash/common/basis/LinearBasis.hpp>
#include <sgpp/base/operation/hash/common/basis/LinearBoundaryBasis.hpp>
#include <sgpp/base/operation/hash/common/basis/LinearClenshawCurtisBasis.hpp>
#include <sgpp/base/operation/hash/common/basis/LinearModifiedBasis.hpp>
#include <sgpp/base/operation/hash/common/basis/WaveletBasis.hpp>
#include <sgpp/base/operation/hash/common/basis/WaveletBoundaryBasis.hpp>
#include <sgpp/base/operation/hash/common/basis/WaveletModifiedBasis.hpp>

#include <sgpp/base/grid/type/LinearClenshawCurtisGrid.hpp>
#include <sgpp/base/grid/type/BsplineGrid.hpp>
#include <sgpp/base/grid/type/BsplineTruncatedBoundaryGrid.hpp>
#include <sgpp/base/grid/type/BsplineClenshawCurtisGrid.hpp>
#include <sgpp/base/grid/type/ModBsplineGrid.hpp>

#include <vector>
#include <cstddef>
#include <cstring>
#include <stdexcept>
#include <memory>

namespace SGPP {
  namespace optimization {
    namespace sle {
      namespace system {

        /**
         * Linear system of the hierarchisation in a sparse grid.
         */
        class Hierarchisation : public Cloneable {
          public:
            /**
             * Constructor.
             * Do not destruct the grid before this object!
             *
             * @param grid              sparse grid
             */
            Hierarchisation(base::Grid& grid)
              : Cloneable(),
                grid(grid),
                gridStorage(grid.getStorage()),
                basisType(INVALID) {
              initialize();
            }

            /**
             * Constructor.
             * Do not destruct the grid before this object!
             *
             * @param grid              sparse grid
             * @param gridStorage       custom grid storage (use basis function according to grid,
             *                          but use another set of grid points according to gridStorage)
             */
            Hierarchisation(base::Grid& grid, base::GridStorage* gridStorage)
              : Cloneable(),
                grid(grid),
                gridStorage(gridStorage),
                basisType(INVALID) {
              initialize();
            }

            /**
             * @param i     row index
             * @param j     column index
             * @return      whether the i-th grid point lies in the support of the j-th basis function
             */
            inline bool isMatrixEntryNonZero(size_t i, size_t j) {
              return (evalBasisFunctionAtGridPoint(j, i) != 0.0);
            }

            /**
             * @param i     row index
             * @param j     column index
             * @return      value of the j-th basis function at the i-th grid point
             */
            inline float_t getMatrixEntry(size_t i, size_t j) {
              return evalBasisFunctionAtGridPoint(j, i);
            }

            /**
             * @return          sparse grid
             */
            base::Grid& getGrid() {
              return grid;
            }

            /**
             * @return grid     sparse grid
             */
            void setGrid(base::Grid& grid) {
              this->grid = grid;
            }

            /**
             * @return              grid storage
             */
            base::GridStorage* getGridStorage() {
              return gridStorage;
            }

            /**
             * @param gridStorage   grid storage (do not destruct before this object!)
             */
            void setGridStorage(base::GridStorage* gridStorage) {
              this->gridStorage = gridStorage;
            }

            size_t getDimension() const {
              return gridStorage->size();
            }

            /**
             * @return smart pointer to cloned object
             */
            virtual Cloneable* clone() {
              return new Hierarchisation(grid, gridStorage);
            }

          protected:
            /// sparse grid
            base::Grid& grid;
            /// grid storage
            base::GridStorage* gridStorage;

            /// B-spline basis
            std::unique_ptr<base::SBsplineBase> bsplineBasis;
            /// B-spline boundary basis
            std::unique_ptr<base::SBsplineBoundaryBase> bsplineBoundaryBasis;
            /// B-spline Clenshaw-Curtis basis
            std::unique_ptr<base::SBsplineClenshawCurtisBase> bsplineClenshawCurtisBasis;
            /// modified B-spline basis
            std::unique_ptr<base::SBsplineModifiedBase> modBsplineBasis;
            /// linear basis
            std::unique_ptr<base::SLinearBase> linearBasis;
            /// linear boundary basis
            std::unique_ptr<base::SLinearBoundaryBase> linearBoundaryBasis;
            /// linear Clenshaw-Curtis basis
            std::unique_ptr<base::SLinearClenshawCurtisBase> linearClenshawCurtisBasis;
            /// modified linear basis
            std::unique_ptr<base::SLinearModifiedBase> modLinearBasis;
            /// wavelet basis
            std::unique_ptr<base::SWaveletBase> waveletBasis;
            /// wavelet boundary basis
            std::unique_ptr<base::SWaveletBoundaryBase> waveletBoundaryBasis;
            /// modified wavelet basis
            std::unique_ptr<base::SWaveletModifiedBase> modWaveletBasis;

            /// type of grid/basis functions
            enum {
              INVALID,
              BSPLINE,
              BSPLINE_BOUNDARY,
              BSPLINE_CLENSHAW_CURTIS,
              BSPLINE_MODIFIED,
              LINEAR,
              LINEAR_BOUNDARY,
              LINEAR_CLENSHAW_CURTIS,
              LINEAR_MODIFIED,
              WAVELET,
              WAVELET_BOUNDARY,
              WAVELET_MODIFIED
            } basisType;

            /**
             * Initialize the correct basis (according to the grid).
             */
            void initialize() {
              if (strcmp(grid.getType(), "bspline") == 0) {
                bsplineBasis = std::unique_ptr<base::SBsplineBase>(
                                 new base::SBsplineBase(
                                   dynamic_cast<base::BsplineGrid&>(grid).getDegree()));
                basisType = BSPLINE;
              } else if (strcmp(grid.getType(), "bsplineTruncatedBoundary") == 0) {
                bsplineBoundaryBasis =
                  std::unique_ptr<base::SBsplineBoundaryBase>(
                    new base::SBsplineBoundaryBase(
                      dynamic_cast<base::BsplineTruncatedBoundaryGrid&>(grid)
                      .getDegree()));
                basisType = BSPLINE_BOUNDARY;
              } else if (strcmp(grid.getType(), "bsplineClenshawCurtis") == 0) {
                bsplineClenshawCurtisBasis = std::unique_ptr <
                                             base::SBsplineClenshawCurtisBase > (
                                               new base::SBsplineClenshawCurtisBase(
                                                 dynamic_cast<base::BsplineClenshawCurtisGrid&>(grid).getDegree()));
                basisType = BSPLINE_CLENSHAW_CURTIS;
              } else if (strcmp(grid.getType(), "modBspline") == 0) {
                modBsplineBasis = std::unique_ptr<base::SBsplineModifiedBase>(
                                    new base::SBsplineModifiedBase(
                                      dynamic_cast<base::ModBsplineGrid&>(grid).getDegree()));
                basisType = BSPLINE_MODIFIED;
              } else if (strcmp(grid.getType(), "linear") == 0) {
                linearBasis = std::unique_ptr<base::SLinearBase>(
                                new base::SLinearBase());
                basisType = LINEAR;
              } else if (strcmp(grid.getType(), "linearTruncatedBoundary") == 0) {
                linearBoundaryBasis = std::unique_ptr<base::SLinearBoundaryBase>(
                                        new base::SLinearBoundaryBase());
                basisType = LINEAR_BOUNDARY;
              } else if (strcmp(grid.getType(), "linearClenshawCurtis") == 0) {
                linearClenshawCurtisBasis = std::unique_ptr <
                                            base::SLinearClenshawCurtisBase > (
                                              new base::SLinearClenshawCurtisBase());
                basisType = LINEAR_CLENSHAW_CURTIS;
              } else if (strcmp(grid.getType(), "modlinear") == 0) {
                modLinearBasis = std::unique_ptr<base::SLinearModifiedBase>(
                                   new base::SLinearModifiedBase());
                basisType = LINEAR_MODIFIED;
              } else if (strcmp(grid.getType(), "wavelet") == 0) {
                waveletBasis = std::unique_ptr<base::SWaveletBase>(
                                 new base::SWaveletBase());
                basisType = WAVELET;
              } else if (strcmp(grid.getType(), "waveletTruncatedBoundary") == 0) {
                waveletBoundaryBasis = std::unique_ptr<base::SWaveletBoundaryBase>(
                                         new base::SWaveletBoundaryBase());
                basisType = WAVELET_BOUNDARY;
              } else if (strcmp(grid.getType(), "modWavelet") == 0) {
                modWaveletBasis = std::unique_ptr<base::SWaveletModifiedBase>(
                                    new base::SWaveletModifiedBase());
                basisType = WAVELET_MODIFIED;
              } else {
                throw std::invalid_argument("Grid type not supported.");
              }
            }

            /**
             * @param basisI    basis function index
             * @param pointJ    grid point index
             * @return          value of the basisI-th basis function at the pointJ-th grid point
             */
            inline float_t evalBasisFunctionAtGridPoint(size_t basisI, size_t pointJ) {
              if (basisType == BSPLINE) {
                return evalBsplineFunctionAtGridPoint(basisI, pointJ);
              } else if (basisType == BSPLINE_BOUNDARY) {
                return evalBsplineBoundaryFunctionAtGridPoint(basisI, pointJ);
              } else if (basisType == BSPLINE_CLENSHAW_CURTIS) {
                return evalBsplineClenshawCurtisFunctionAtGridPoint(basisI, pointJ);
              } else if (basisType == BSPLINE_MODIFIED) {
                return evalBsplineModifiedFunctionAtGridPoint(basisI, pointJ);
              } else if (basisType == LINEAR) {
                return evalLinearFunctionAtGridPoint(basisI, pointJ);
              } else if (basisType == LINEAR_BOUNDARY) {
                return evalLinearBoundaryFunctionAtGridPoint(basisI, pointJ);
              } else if (basisType == LINEAR_CLENSHAW_CURTIS) {
                return evalLinearClenshawCurtisFunctionAtGridPoint(basisI, pointJ);
              } else if (basisType == LINEAR_MODIFIED) {
                return evalLinearModifiedFunctionAtGridPoint(basisI, pointJ);
              } else if (basisType == WAVELET) {
                return evalWaveletFunctionAtGridPoint(basisI, pointJ);
              } else if (basisType == WAVELET_BOUNDARY) {
                return evalWaveletBoundaryFunctionAtGridPoint(basisI, pointJ);
              } else if (basisType == WAVELET_MODIFIED) {
                return evalWaveletModifiedFunctionAtGridPoint(basisI, pointJ);
              } else {
                return 0.0;
              }
            }

            /**
             * @param basisI    basis function index
             * @param pointJ    grid point index
             * @return          value of the basisI-th B-spline basis function
             *                  at the pointJ-th grid point
             */
            inline float_t evalBsplineFunctionAtGridPoint(size_t basisI, size_t pointJ) {
              const base::GridIndex* gpBasis = gridStorage->get(basisI);
              const base::GridIndex* gpPoint = gridStorage->get(pointJ);
              float_t result = 1.0;

              for (size_t t = 0; t < gridStorage->dim(); t++) {
                float_t result1d = bsplineBasis->eval(gpBasis->getLevel(t),
                                                      gpBasis->getIndex(t),
                                                      gpPoint->getCoord(t));

                if (result1d == 0.0) {
                  return 0.0;
                }

                result *= result1d;
              }

              return result;
            }

            /**
             * @param basisI    basis function index
             * @param pointJ    grid point index
             * @return          value of the basisI-th B-spline boundary basis function
             *                  at the pointJ-th grid point
             */
            inline float_t evalBsplineBoundaryFunctionAtGridPoint(size_t basisI,
                size_t pointJ) {
              const base::GridIndex* gpBasis = gridStorage->get(basisI);
              const base::GridIndex* gpPoint = gridStorage->get(pointJ);
              float_t result = 1.0;

              for (size_t t = 0; t < gridStorage->dim(); t++) {
                float_t result1d = bsplineBoundaryBasis->eval(gpBasis->getLevel(t),
                                   gpBasis->getIndex(t),
                                   gpPoint->getCoord(t));

                if (result1d == 0.0) {
                  return 0.0;
                }

                result *= result1d;
              }

              return result;
            }

            /**
             * @param basisI    basis function index
             * @param pointJ    grid point index
             * @return          value of the basisI-th B-spline Clenshaw-Curtis basis function
             *                  at the pointJ-th grid point
             */
            inline float_t evalBsplineClenshawCurtisFunctionAtGridPoint(size_t basisI,
                size_t pointJ) {
              const base::GridIndex* gpBasis = gridStorage->get(basisI);
              const base::GridIndex* gpPoint = gridStorage->get(pointJ);
              float_t result = 1.0;

              for (size_t t = 0; t < gridStorage->dim(); t++) {
                float_t result1d = bsplineClenshawCurtisBasis->eval(
                                     gpBasis->getLevel(t), gpBasis->getIndex(t), gpPoint->getCoord(t));

                if (result1d == 0.0) {
                  return 0.0;
                }

                result *= result1d;
              }

              return result;
            }

            /**
             * @param basisI    basis function index
             * @param pointJ    grid point index
             * @return          value of the basisI-th modified B-spline basis function
             *                  at the pointJ-th grid point
             */
            inline float_t evalBsplineModifiedFunctionAtGridPoint(size_t basisI,
                size_t pointJ) {
              const base::GridIndex* gpBasis = gridStorage->get(basisI);
              const base::GridIndex* gpPoint = gridStorage->get(pointJ);
              float_t result = 1.0;

              for (size_t t = 0; t < gridStorage->dim(); t++) {
                float_t result1d = modBsplineBasis->eval(gpBasis->getLevel(t),
                                   gpBasis->getIndex(t),
                                   gpPoint->getCoord(t));

                if (result1d == 0.0) {
                  return 0.0;
                }

                result *= result1d;
              }

              return result;
            }

            /**
             * @param basisI    basis function index
             * @param pointJ    grid point index
             * @return          value of the basisI-th linear basis function
             *                  at the pointJ-th grid point
             */
            inline float_t evalLinearFunctionAtGridPoint(size_t basisI, size_t pointJ) {
              const base::GridIndex* gpBasis = gridStorage->get(basisI);
              const base::GridIndex* gpPoint = gridStorage->get(pointJ);
              float_t result = 1.0;

              for (size_t t = 0; t < gridStorage->dim(); t++) {
                float_t result1d = linearBasis->eval(gpBasis->getLevel(t),
                                                     gpBasis->getIndex(t),
                                                     gpPoint->getCoord(t));

                if (result1d == 0.0) {
                  return 0.0;
                }

                result *= result1d;
              }

              return result;
            }

            /**
             * @param basisI    basis function index
             * @param pointJ    grid point index
             * @return          value of the basisI-th linear boundary basis function
             *                  at the pointJ-th grid point
             */
            inline float_t evalLinearBoundaryFunctionAtGridPoint(size_t basisI,
                size_t pointJ) {
              const base::GridIndex* gpBasis = gridStorage->get(basisI);
              const base::GridIndex* gpPoint = gridStorage->get(pointJ);
              float_t result = 1.0;

              for (size_t t = 0; t < gridStorage->dim(); t++) {
                float_t result1d = linearBoundaryBasis->eval(gpBasis->getLevel(t),
                                   gpBasis->getIndex(t),
                                   gpPoint->getCoord(t));

                if (result1d == 0.0) {
                  return 0.0;
                }

                result *= result1d;
              }

              return result;
            }

            /**
             * @param basisI    basis function index
             * @param pointJ    grid point index
             * @return          value of the basisI-th linear Clenshaw-Curtis basis function
             *                  at the pointJ-th grid point
             */
            inline float_t evalLinearClenshawCurtisFunctionAtGridPoint(size_t basisI,
                size_t pointJ) {
              const base::GridIndex* gpBasis = gridStorage->get(basisI);
              const base::GridIndex* gpPoint = gridStorage->get(pointJ);
              float_t result = 1.0;

              for (size_t t = 0; t < gridStorage->dim(); t++) {
                float_t result1d = linearClenshawCurtisBasis->eval(
                                     gpBasis->getLevel(t), gpBasis->getIndex(t), gpPoint->getCoord(t));

                if (result1d == 0.0) {
                  return 0.0;
                }

                result *= result1d;
              }

              return result;
            }

            /**
             * @param basisI    basis function index
             * @param pointJ    grid point index
             * @return          value of the basisI-th modified linear basis function
             *                  at the pointJ-th grid point
             */
            inline float_t evalLinearModifiedFunctionAtGridPoint(size_t basisI,
                size_t pointJ) {
              const base::GridIndex* gpBasis = gridStorage->get(basisI);
              const base::GridIndex* gpPoint = gridStorage->get(pointJ);
              float_t result = 1.0;

              for (size_t t = 0; t < gridStorage->dim(); t++) {
                float_t result1d = modLinearBasis->eval(gpBasis->getLevel(t),
                                                        gpBasis->getIndex(t),
                                                        gpPoint->getCoord(t));

                if (result1d == 0.0) {
                  return 0.0;
                }

                result *= result1d;
              }

              return result;
            }

            /**
             * @param basisI    basis function index
             * @param pointJ    grid point index
             * @return          value of the basisI-th wavelet basis function
             *                  at the pointJ-th grid point
             */
            inline float_t evalWaveletFunctionAtGridPoint(size_t basisI, size_t pointJ) {
              const base::GridIndex* gpBasis = gridStorage->get(basisI);
              const base::GridIndex* gpPoint = gridStorage->get(pointJ);
              float_t result = 1.0;

              for (size_t t = 0; t < gridStorage->dim(); t++) {
                float_t result1d = waveletBasis->eval(gpBasis->getLevel(t),
                                                      gpBasis->getIndex(t),
                                                      gpPoint->getCoord(t));

                if (result1d == 0.0) {
                  return 0.0;
                }

                result *= result1d;
              }

              return result;
            }

            /**
             * @param basisI    basis function index
             * @param pointJ    grid point index
             * @return          value of the basisI-th wavelet boundary basis function
             *                  at the pointJ-th grid point
             */
            inline float_t evalWaveletBoundaryFunctionAtGridPoint(size_t basisI,
                size_t pointJ) {
              const base::GridIndex* gpBasis = gridStorage->get(basisI);
              const base::GridIndex* gpPoint = gridStorage->get(pointJ);
              float_t result = 1.0;

              for (size_t t = 0; t < gridStorage->dim(); t++) {
                float_t result1d = waveletBoundaryBasis->eval(gpBasis->getLevel(t),
                                   gpBasis->getIndex(t),
                                   gpPoint->getCoord(t));

                if (result1d == 0.0) {
                  return 0.0;
                }

                result *= result1d;
              }

              return result;
            }

            /**
             * @param basisI    basis function index
             * @param pointJ    grid point index
             * @return          value of the basisI-th modified wavelet basis function
             *                  at the pointJ-th grid point
             */
            inline float_t evalWaveletModifiedFunctionAtGridPoint(size_t basisI,
                size_t pointJ) {
              const base::GridIndex* gpBasis = gridStorage->get(basisI);
              const base::GridIndex* gpPoint = gridStorage->get(pointJ);
              float_t result = 1.0;

              for (size_t t = 0; t < gridStorage->dim(); t++) {
                float_t result1d = modWaveletBasis->eval(gpBasis->getLevel(t),
                                   gpBasis->getIndex(t),
                                   gpPoint->getCoord(t));

                if (result1d == 0.0) {
                  return 0.0;
                }

                result *= result1d;
              }

              return result;
            }
        };

      }
    }
  }
}

#endif
