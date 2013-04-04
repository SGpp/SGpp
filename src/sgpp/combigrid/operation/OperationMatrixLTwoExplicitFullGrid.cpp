/* ****************************************************************************
* Copyright (C) 2011 Technische Universitaet Muenchen                         *
* This file is part of the SG++ project. For conditions of distribution and   *
* use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
**************************************************************************** */
// @author Janos Benk (benk@in.tum.de)

#include "OperationMatrixLTwoExplicitFullGrid.hpp"
#include "base/exception/data_exception.hpp"
#include <vector>
#include <cmath>
namespace combigrid {

  OperationMatrixLTwoExplicitFullGrid::OperationMatrixLTwoExplicitFullGrid(sg::base::DataMatrix* m,
      combigrid::FullGrid<double>* grid) :
    ownsMatrix_(false) {
    m_ = m;
    buildMatrix(grid);
  }

  OperationMatrixLTwoExplicitFullGrid::OperationMatrixLTwoExplicitFullGrid(
    combigrid::FullGrid<double>* grid) :
    ownsMatrix_(true) {
    m_ = new sg::base::DataMatrix(grid->getNrElements(), grid->getNrElements());
    buildMatrix(grid);
  }

  OperationMatrixLTwoExplicitFullGrid::~OperationMatrixLTwoExplicitFullGrid() {
    if (ownsMatrix_)
      delete m_;
  }

  void OperationMatrixLTwoExplicitFullGrid::buildMatrix(
    combigrid::FullGrid<double>* grid) {
    unsigned int gridSize = grid->getNrElements();
    unsigned int gridDim = grid->getDimension();

    std::vector<int> levels = grid->getLevels();

    for (unsigned int i = 0; i < gridSize; i++) {
      std::vector<int> i_indexes(gridDim);
      grid->getVectorIndex(i, i_indexes);

      for (unsigned int j = i; j < gridSize; j++) {
        std::vector<int> j_indexes(gridDim);
        grid->getVectorIndex(j, j_indexes);
        double res = 1;

        for (unsigned int k = 0; k < gridDim; k++) {
          if (i_indexes[k] == j_indexes[k]) {
            //Use formula for identical ansatz functions:
            res *= 2 / pow(2, levels[k]) / 3;
          } else if (//std::abs(i_indexes[k] - j_indexes[k]) == 1
            i_indexes[k] - j_indexes[k] == 1
            || i_indexes[k] - j_indexes[k] == -1 ) { //workaround as abs(int) is not supported in gcc 4.4
            //Ansatz functions overlap:
            res *= 1 / pow(2, levels[k]) / 6;
          } else {
            //Ansatz functions do not overlap => result is zero
            res = 0;
            break;
          }
        }

        m_->set(i, j, res);
        m_->set(j, i, res);
      }
    }
  }

  void OperationMatrixLTwoExplicitFullGrid::mult(sg::base::DataVector& alpha,
      sg::base::DataVector& result) {

    unsigned int nrows( static_cast<unsigned int>( m_->getNrows() ) );
    unsigned int ncols( static_cast<unsigned int>( m_->getNcols() ) );

    if (alpha.getSize() != ncols || result.getSize() != nrows) {
      throw sg::base::data_exception("Dimensions do not match!");
    }

    double* data = m_->getPointer();

    //Standard matrix multiplication:
    double temp = 0.;
    unsigned int acc = 0;

    for (unsigned int i = 0; i < nrows; i++) {
      for (unsigned int j = 0; j < ncols; j++) {
        temp += data[j + acc] * alpha[j];
      }

      result[i] = temp;
      temp = 0.;
      acc += ncols;
    }
  }

} /* namespace combigrid */
