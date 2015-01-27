// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at 
// sgpp.sparsegrids.org

#include <sgpp/base/datatypes/DataVector.hpp>
#include <sgpp/base/datatypes/DataMatrix.hpp>
#include <sgpp/base/operation/hash/OperationMatrix.hpp>
#include <sgpp/combigrid/fullgrid/CombiFullGrid.hpp>

namespace combigrid {

  /**
   * Explicit representation of the matrix \f$(\Phi_i,\Phi_j)_{LTwo}\f$ for a full grid
   */
  class OperationMatrixLTwoExplicitFullGrid: public SGPP::base::OperationMatrix {
    public:

      /**
       * Constructor that uses a external matrix pointer to construct the matrix,
       * i.e. matrix is NOT destroyed by the destructor of OperationMatrixLTwoExplicitFullGrid
       *
       * @param m pointer to datamatrix of size (number of grid point) x (number of grid points)
       * @param grid the full grid
       */
      OperationMatrixLTwoExplicitFullGrid(SGPP::base::DataMatrix* m,
                                          combigrid::FullGrid<double>* grid);

      /**
       * Constructor that creates an own matrix
       * i.e. matrix is destroyed by the destructor of OperationMatrixLTwoExplicitFullGrid
       *
       * @param grid the full grid
       */
      OperationMatrixLTwoExplicitFullGrid(combigrid::FullGrid<double>* grid);

      /**
       * Destructor
       */
      virtual ~OperationMatrixLTwoExplicitFullGrid();

      /**
       * Implementation of standard matrix multiplication
       *
       * @param alpha DataVector that is multiplied to the matrix
       * @param result DataVector into which the result of multiplication is stored
       */
      virtual void mult(SGPP::base::DataVector& alpha,
                        SGPP::base::DataVector& result);

    private:
      /**
       * This method is used by both constructors to build the matrix
       */
      void buildMatrix(combigrid::FullGrid<double>* grid);

      SGPP::base::DataMatrix* m_;
      bool ownsMatrix_;
  };

} /* namespace combigrid */